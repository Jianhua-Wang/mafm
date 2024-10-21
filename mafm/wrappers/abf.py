"""Warpper of ABF fine-mapping method."""

import logging

import numpy as np
import pandas as pd

from mafm.constants import ColName

logger = logging.getLogger("ABF")


def run_abf(
    sumstats: pd.DataFrame, var_prior: float = 0.2, max_causal: int = 1, **kwargs
) -> pd.Series:
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
    sumstats : pd.DataFrame
        Summary statistics.
    var_prior : float, optional
        Variance prior, by default 0.2, usually set to 0.15 for quantitative traits
        and 0.2 for binary traits.
    max_causal : int, optional
        Maximum number of causal variants, by default 1

    Returns
    -------
    pd.Series
        The result of ABF.
    """
    if max_causal > 1:
        logger.warning("ABF only support single causal variant.")
        max_causal = 1
    df = sumstats.copy()
    df["W2"] = var_prior**2
    df["SNP_BF"] = np.sqrt(
        (df[ColName.SE] ** 2 / (df[ColName.SE] ** 2 + df["W2"]))
    ) * np.exp(
        df["W2"]
        / (df[ColName.BETA] ** 2 + df["W2"])
        * (df[ColName.BETA] ** 2 / df[ColName.SE] ** 2)
        / 2
    )
    df[ColName.PP_ABF] = df["SNP_BF"] / df["SNP_BF"].sum()
    return pd.Series(data=df[ColName.PP_ABF].values, index=df[ColName.SNPID].tolist())
