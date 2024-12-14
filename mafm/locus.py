"""Class for the input data of the fine-mapping analysis."""

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName
from mafm.ldmatrix import LDMatrix, load_ld
from mafm.sumstats import load_sumstats

logger = logging.getLogger("Locus")


class Locus:
    """
    Locus class to represent a genomic locus with associated summary statistics and linkage disequilibrium (LD) matrix.

    Attributes
    ----------
    original_sumstats : pd.DataFrame
        The original summary statistics file.
        Population code.
    chrom : int
        Chromosome.
    start : int
        Start position of the locus.
    end : int
        End position of the locus.
    n_snps : int
        Number of SNPs in the locus.
    locus_id : str
        Unique identifier for the locus.
    is_matched : bool
        Whether the LD matrix and summary statistics file are matched.

    Methods
    -------
    __init__(popu: str, cohort: str, sample_size: int, sumstats: pd.DataFrame, ld: Optional[LDMatrix] = None, if_intersect: bool = False)
    __repr__()
        Return a string representation of the Locus object.
    copy()
        Copy the Locus object.
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
                inters = intersect_sumstat_ld(self)
                self.sumstats = inters.sumstats
                self.ld = inters.ld
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
        return f"Locus(popu={self.popu}, cohort={self.cohort}, sample_size={self.sample_size}, chr={self.chrom}, start={self.start}, end={self.end}, sumstats={self.sumstats.shape}, ld={self.ld.r.shape})"

    def copy(self):
        """Copy the Locus object."""
        return Locus(self.popu, self.cohort, self.sample_size, self.sumstats.copy(), self.ld.copy(), if_intersect=False)


class LocusSet:
    """
    LocusSet class to represent a set of genomic loci.

    Attributes
    ----------
    loci : list[Locus]
        List of Locus objects.
    n_loci : int
        Number of loci.
    chrom : int
        Chromosome number.
    start : int
        Start position of the locus.
    end : int
        End position of the locus.
    locus_id : str
        Unique identifier for the locus.

    Methods
    -------
    __init__(loci: list[Locus])
    __repr__()
        Return a string representation of the LocusSet object.
    copy()
        Copy the LocusSet object.
    """

    def __init__(self, loci: list[Locus]):
        """
        Initialize the LocusSet object.

        Parameters
        ----------
        loci : list[Locus]
            List of Locus objects.

        """
        self.loci = loci

    @property
    def n_loci(self):
        """Get the number of loci."""
        return len(self.loci)

    @property
    def chrom(self):
        """Get the chromosome."""
        chrom_list = [locus.chrom for locus in self.loci]
        if len(set(chrom_list)) > 1:
            raise ValueError("The chromosomes of the loci are not the same.")
        return chrom_list[0]

    @property
    def start(self):
        """Get the start position."""
        return min([locus.start for locus in self.loci])

    @property
    def end(self):
        """Get the end position."""
        return max([locus.end for locus in self.loci])

    @property
    def locus_id(self):
        """Get the locus ID."""
        return f"{self.chrom}:{self.start}-{self.end}"

    def __repr__(self):
        """Return a string representation of the LocusSet object."""
        return (
            f"LocusSet(\n n_loci={len(self.loci)}, chrom={self.chrom}, start={self.start}, end={self.end}, locus_id={self.locus_id} \n"
            + "\n".join([locus.__repr__() for locus in self.loci])
            + "\n"
            + ")"
        )

    def copy(self):
        """Copy the LocusSet object."""
        return LocusSet([locus.copy() for locus in self.loci])


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
        return locus
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


def intersect_loci(list_loci: list[Locus]) -> list[Locus]:
    """
    Intersect the Variant IDs in the LD matrices and the sumstats files of a list of Locus objects.

    Parameters
    ----------
    list_loci : list[Locus]
        List of Locus objects.

    Returns
    -------
    list[Locus]
        List of Locus objects containing the intersected LD matrices and sumstats files.
    """
    raise NotImplementedError("Intersect the Variant IDs in the LD matrices and the sumstats files of a list of Locus objects.")


def load_locus(prefix: str, popu: str, cohort: str, sample_size: int, if_intersect: bool = False, **kwargs) -> Locus:
    """
    Load the input data of the fine-mapping analysis.

    Parameters
    ----------
    prefix : str
        Prefix of the input files.
    popu : str
        Population of the input data.
    cohort : str
        Cohort of the input data.
    sample_size : int
        Sample size of the input data.
    if_intersect : bool, optional
        Whether to intersect the input data with the LD matrix, by default False.

    Returns
    -------
    Locus
        Object containing the input data.

    Raises
    ------
    ValueError
        If the input files are not found.
    """
    sumstats = load_sumstats(f"{prefix}.sumstat", if_sort_alleles=True, **kwargs)
    if os.path.exists(f"{prefix}.ld"):
        ld = load_ld(f"{prefix}.ld", f"{prefix}.ldmap", if_sort_alleles=True, **kwargs)
    elif os.path.exists(f"{prefix}.ld.npz"):
        ld = load_ld(f"{prefix}.ld.npz", f"{prefix}.ldmap", if_sort_alleles=True, **kwargs)
    else:
        raise ValueError("LD matrix file not found.")

    return Locus(popu, cohort, sample_size, sumstats=sumstats, ld=ld, if_intersect=if_intersect)


def load_locus_set(locus_info: pd.DataFrame, if_intersect: bool = False, **kwargs) -> LocusSet:
    """
    Load the input data of the fine-mapping analysis.

    Parameters
    ----------
    locus_info : pd.DataFrame
        Dataframe containing the locus information.
    if_intersect : bool, optional
        Whether to intersect the input data with the LD matrix, by default False.

    Returns
    -------
    LocusSet
        Object containing the input data.
    """
    required_cols = ["prefix", "popu", "cohort", "sample_size"]
    missing_cols = [col for col in required_cols if col not in locus_info.columns]
    if len(missing_cols) > 0:
        raise ValueError(f"The following columns are required: {missing_cols}")
    # TODO: make sure the combination of popu and cohort is unique
    loci = []
    for i, row in locus_info.iterrows():
        loci.append(load_locus(row["prefix"], row["popu"], row["cohort"], row["sample_size"], if_intersect, **kwargs))
    return LocusSet(loci)
