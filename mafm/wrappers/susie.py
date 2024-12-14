"""Warpper of SuSiE."""

import json
import logging

import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
from mafm.locus import Locus, intersect_sumstat_ld
from mafm.wrappers.susie_rss import susie_rss

logger = logging.getLogger("SuSiE")


def run_susie(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    max_iter: int = 100,
    estimate_residual_variance: bool = False,
    min_abs_corr: float = 0.5,
    convergence_tol: float = 1e-3,
):
    """
    Run SuSiE with shotgun stochastic search.

    Parameters
    ----------
    locus : Locus
        Locus object.
    max_causal : int, optional
        Maximum number of causal variants, by default 1.
    coverage : float, optional
        Coverage of the credible set, by default 0.95.
    max_iter : int, optional
        Maximum number of iterations, by default 100.
    estimate_residual_variance : bool, optional
        Whether to estimate residual variance, by default False.
    min_abs_corr : float, optional
        Minimum absolute correlation, by default 0.5.
    convergence_tol : float, optional
        Convergence tolerance, by default 1e-3.

    Returns
    -------
    CredibleSet
        Credible set.
    """
    if not locus.is_matched:
        logger.warning("The sumstat and LD are not matched, will match them in same order.")
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running SuSiE on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "max_iter": max_iter,
        "estimate_residual_variance": estimate_residual_variance,
        "min_abs_corr": min_abs_corr,
        "convergence_tol": convergence_tol,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    s = susie_rss(
        bhat=locus.sumstats[ColName.BETA].to_numpy(),
        shat=locus.sumstats[ColName.SE].to_numpy(),
        n=locus.sample_size,
        R=locus.ld.r,
        L=max_causal,
        coverage=coverage,
        max_iter=max_iter,
        estimate_residual_variance=estimate_residual_variance,
        min_abs_corr=min_abs_corr,
        tol=convergence_tol,
    )
    pip = s["pip"]
    cs_idx = list(s["sets"]["cs"].values())
    n_cs = len(cs_idx)
    cs_sizes = [len(idx) for idx in cs_idx]
    cred_snps = [locus.sumstats[ColName.SNPID].iloc[idx].tolist() for idx in cs_idx]
    pips = pd.Series(data=pip, index=locus.sumstats[ColName.SNPID].tolist())
    lead_snps = [str(pips[pips.index.isin(cred_snps[i])].idxmax()) for i in range(n_cs)]
    logger.info(f"Fished SuSiE on {locus}")
    logger.info(f"N of credible set: {n_cs}")
    logger.info(f"Credible set size: {cs_sizes}")
    return CredibleSet(
        tool=Method.SUSIE,
        n_cs=n_cs,
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cred_snps,
        cs_sizes=cs_sizes,
        pips=pips,
        parameters=parameters,
    )
