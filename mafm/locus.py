"""Class for the input data of the fine-mapping analysis."""

import logging
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName
from mafm.ldmatrix import LDMatrix

logger = logging.getLogger("Locus")


class Locus:
    """
    Class for the input data of the fine-mapping analysis.

    Attributes
    ----------
    sumstats : pd.DataFrame
        Sumstats file.
    original_sumstats : pd.DataFrame
        Original sumstats file.
    ld : LDMatrix
        LD matrix.
    popu : str
        Population code.
    sample_size : int
        Sample size.
    """

    def __init__(
        self,
        popu: str,
        cohort: str,
        sample_size: int,
        sumstats: pd.DataFrame,
        ld: Optional[LDMatrix] = None,
        if_intersect: bool = False,
    ):
        """
        Initialize the Locus object.

        Parameters
        ----------
        popu : str
            Population code. e.g. "EUR". Choose from ["AFR", "AMR", "EAS", "EUR", "SAS"].
        cohort : str
            Cohort name.
        sample_size : int
            Sample size.
        sumstats : pd.DataFrame
            Sumstats file.
        ld : LDMatrix, optional
            LD matrix.
        if_intersect : bool, optional
            Whether to intersect the LD matrix and sumstats file, by default True.

        """
        self.sumstats = sumstats
        self._original_sumstats = self.sumstats.copy()
        self._popu = popu
        self._cohort = cohort
        self._sample_size = sample_size
        if ld:
            self.ld = ld
            if if_intersect:
                self = intersect_sumstat_ld(self)
        else:
            logger.warning("LD matrix and map file not found. Can only run ABF method.")
            self.ld = LDMatrix(pd.DataFrame(), np.array([]))

    @property
    def original_sumstats(self):
        """Get the original sumstats file."""
        return self._original_sumstats

    @property
    def popu(self):
        """Get the population code."""
        return self._popu

    @property
    def cohort(self):
        """Get the cohort name."""
        return self._cohort

    @property
    def sample_size(self):
        """Get the sample size."""
        return self._sample_size

    @property
    def chrom(self):
        """Get the chromosome."""
        return self.sumstats[ColName.CHR].iloc[0]

    @property
    def start(self):
        """Get the start position."""
        return self.sumstats[ColName.BP].min()

    @property
    def end(self):
        """Get the end position."""
        return self.sumstats[ColName.BP].max()

    @property
    def n_snps(self):
        """Get the number of SNPs."""
        return len(self.sumstats)

    @property
    def locus_id(self):
        """Get the locus ID."""
        return f"{self.popu}_{self.cohort}_chr{self.chrom}:{self.start}-{self.end}"

    @property
    def is_matched(self):
        """Check if the LD matrix and sumstats file are matched."""
        # check the order of SNPID in the LD matrix and the sumstats file are the exact same
        if self.ld is None:
            return False
        return self.ld.map[ColName.SNPID].equals(self.sumstats[ColName.SNPID])

    def __repr__(self):
        """Return a string representation of the Locus object."""
        return f"Locus(locus_id={self.locus_id}, popu={self.popu}, cohort={self.cohort}, sample_size={self.sample_size}, sumstats={self.sumstats.shape}, ld={self.ld.r.shape})"

    def copy(self):
        """Copy the Locus object."""
        return Locus(self.popu, self.cohort, self.sample_size, self.sumstats.copy(), self.ld.copy(), if_intersect=False)


class LocusCollection:
    """
    Class to store a collection of Locus objects.

    All the loci are located in the same chromosome and the same region.

    Attributes
    ----------
    loci : list[Locus]
        List of Locus objects.
    """

    def __init__(self, loci: list[Locus]):
        self.loci = loci

def intersect_sumstat_ld(locus: Locus) -> Locus:
    """
    Intersect the Variant IDs in the LD matrix and the sumstats file.

    Raises
    ------
    ValueError
        If no common Variant IDs found between the LD matrix and the sumstats file.

    Returns
    -------
    Locus
        Object containing the intersected LD matrix and sumstats file.
    """
    if locus.ld is None:
        raise ValueError("LD matrix not found.")
    if locus.is_matched:
        logger.info("The LD matrix and sumstats file are matched.")
    ldmap = locus.ld.map.copy()
    r = locus.ld.r.copy()
    sumstats = locus.sumstats.copy()
    sumstats = sumstats.sort_values([ColName.CHR, ColName.BP], ignore_index=True)
    intersec_sumstats = sumstats[sumstats[ColName.SNPID].isin(ldmap[ColName.SNPID])].copy()
    intersec_variants = intersec_sumstats[ColName.SNPID].to_numpy()
    if len(intersec_variants) == 0:
        raise ValueError("No common Variant IDs found between the LD matrix and the sumstats file.")
    elif len(intersec_variants) <= 10:
        logger.warning("Only a few common Variant IDs found between the LD matrix and the sumstats file(<= 10).")
    ldmap["idx"] = ldmap.index
    ldmap.set_index(ColName.SNPID, inplace=True, drop=False)
    ldmap = ldmap.loc[intersec_variants].copy()
    intersec_index = ldmap["idx"].to_numpy()
    r = r[intersec_index, :][:, intersec_index]
    intersec_sumstats.reset_index(drop=True, inplace=True)
    ldmap.drop("idx", axis=1, inplace=True)
    ldmap = ldmap.reset_index(drop=True)
    intersec_ld = LDMatrix(ldmap, r)
    logger.info(
        "Intersected the Variant IDs in the LD matrix and the sumstats file. "
        f"Number of common Variant IDs: {len(intersec_index)}"
    )
    return Locus(locus.popu, locus.cohort, locus.sample_size, intersec_sumstats, intersec_ld)
