"""Warpper of ABF fine-mapping method."""

import json
import logging
from typing import List

import numpy as np
import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet, combine_creds
from mafm.locus import Locus

logger = logging.getLogger("ABF")


def run_abf(locus: Locus, max_causal: int = 1, coverage: float = 0.95, var_prior: float = 0.2) -> CredibleSet:
    """
    Run ABF.

    calculate the approximate Bayes factor (ABF) from BETA and SE, using the
    formula:
    SNP_BF = sqrt(SE/(SE + W^2))EXP(W^2/(SE + W^2)*(BETA^2/SE^2)/2)
    where W is variance prior, usually set to 0.15 for quantitative traits
    and 0.2 for binary traits.
    the posterior probability of each variant being causal is calculated
    using the formula:
    PP(causal) = SNP_BF / sum(all_SNP_BFs)

    Reference: Asimit, J. L. et al. Eur J Hum Genet (2016)

    Parameters
    ----------
    locus : Locus
        Locus object.
    max_causal : int, optional
        Maximum number of causal variants, by default 1, only support 1.
    coverage : float, optional
        Coverage, by default 0.95
    var_prior : float, optional
        Variance prior, by default 0.2, usually set to 0.15 for quantitative traits
        and 0.2 for binary traits.


    Returns
    -------
    CredibleSet
        Credible set.
    """
    if max_causal > 1:
        logger.warning("ABF only support single causal variant. max_causal is set to 1.")
        max_causal = 1
    logger.info(f"Running ABF on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "var_prior": var_prior,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    df = locus.original_sumstats.copy()
    df["W2"] = var_prior**2
    df["SNP_BF"] = np.sqrt((df[ColName.SE] ** 2 / (df[ColName.SE] ** 2 + df["W2"]))) * np.exp(
        df["W2"] / (df[ColName.BETA] ** 2 + df["W2"]) * (df[ColName.BETA] ** 2 / df[ColName.SE] ** 2) / 2
    )
    df[ColName.PIP] = df["SNP_BF"] / df["SNP_BF"].sum()
    pips = pd.Series(data=df[ColName.PIP].values, index=df[ColName.SNPID].tolist(), name=ColName.ABF)
    ordering = np.argsort(pips.to_numpy())[::-1]
    idx = np.where(np.cumsum(pips.to_numpy()[ordering]) > coverage)[0][0]
    cs_snps = pips.index[ordering][: (idx + 1)].to_list()
    lead_snps = str(df.loc[df[df[ColName.SNPID].isin(cs_snps)][ColName.P].idxmin(), ColName.SNPID])
    logger.info(f"Fished ABF on {locus}")
    logger.info("N of credible set: 1")
    logger.info(f"Credible set size: {len(cs_snps)}")
    return CredibleSet(
        tool=Method.ABF,
        n_cs=1,
        coverage=coverage,
        lead_snps=[lead_snps],
        snps=[cs_snps],
        cs_sizes=[len(cs_snps)],
        pips=pips,
        parameters=parameters,
    )
