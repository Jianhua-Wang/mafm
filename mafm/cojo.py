"""Wrapper for COJO."""

import logging

import pandas as pd
from cojopy.cojopy import COJO

from mafm.locus import Locus
from mafm.sumstats import ColName

logger = logging.getLogger("COJO")


def conditional_selection(
    locus: Locus,
    p_cutoff: float = 5e-8,
    collinear_cutoff: float = 0.9,
    window_size: int = 10000000,
    maf_cutoff: float = 0.01,
    diff_freq_cutoff: float = 0.2,
) -> pd.DataFrame:
    """
    Perform conditional selection on the locus.

    Parameters
    ----------
    locus : Locus
        The locus to perform conditional selection on.
    p_cutoff : float, optional
        The p-value cutoff for the conditional selection, by default 5e-8.
    collinear_cutoff : float, optional
        The collinearity cutoff for the conditional selection, by default 0.9.
    window_size : int, optional
        The window size for the conditional selection, by default 10000000.
    maf_cutoff : float, optional
        The MAF cutoff for the conditional selection, by default 0.01.
    diff_freq_cutoff : float, optional
        The difference in frequency cutoff for the conditional selection, by default 0.2.

    Returns
    -------
    pd.DataFrame
        The conditional selection results.
    """
    sumstats = locus.sumstats.copy()
    sumstats = sumstats[[ColName.SNPID, ColName.EA, ColName.NEA, ColName.BETA, ColName.SE, ColName.P, ColName.EAF]]
    sumstats.columns = ["SNP", "A1", "A2", "b", "se", "p", "freq"]
    sumstats["N"] = locus.sample_size
    if p_cutoff < 1e-5 and len(sumstats[sumstats["p"] < p_cutoff]) == 0:
        logger.warning("No SNPs passed the p-value cutoff, using p_cutoff=1e-5")
        p_cutoff = 1e-5

    ld_matrix = locus.ld.r.copy()
    ld_freq = locus.ld.map.copy()
    if "AF2" not in ld_freq.columns:
        logger.warning("AF2 is not in the LD matrix.")
        ld_freq = None
    else:
        ld_freq = ld_freq[["SNPID", "AF2"]]
        ld_freq.columns = ["SNP", "freq"]
        ld_freq["freq"] = 1 - ld_freq["freq"]
    c = COJO(
        p_cutoff=p_cutoff,
        collinear_cutoff=collinear_cutoff,
        window_size=window_size,
        maf_cutoff=maf_cutoff,
        diff_freq_cutoff=diff_freq_cutoff,
    )
    c.load_sumstats(sumstats=sumstats, ld_matrix=ld_matrix, ld_freq=ld_freq)  # type: ignore
    cojo_result = c.conditional_selection()
    return cojo_result
