"""Main module."""

import logging

import numpy as np
import pandas as pd

from mafm.constants import ColName
from mafm.ldmatrix import LDMatrix, load_ld
from mafm.sumstats import load_sumstats

logger = logging.getLogger("MAFM")


class FmInput:
    """
    Class for the input data of the fine-mapping analysis.

    Attributes
    ----------
    ld_path : str
        Path to the input LD matrix.
    map_path : str
        Path to the input map file.
    sumstats_path : str
        Path to the input sumstats file.
    r : np.ndarray
        LD matrix.
    map : pd.DataFrame
        Map file.
    sumstats : pd.DataFrame
        Sumstats file.
    """

    def __init__(
        self,
        ld_path: str,
        map_path: str,
        sumstats_path: str,
        if_intersect: bool = True,
        **kwargs,
    ):
        """
        Initialize the FmInput object.

        Parameters
        ----------
        ld_path : str
            Path to the input LD matrix.
        map_path : str
            Path to the input map file.
        sumstats_path : str

        """
        self.ld_path = ld_path
        self.map_path = map_path
        self.sumstats_path = sumstats_path
        ld = load_ld(ld_path, map_path, if_sort_alleles=True, **kwargs)
        self.r = ld.r
        self.map = ld.map
        self.sumstats = load_sumstats(sumstats_path, if_sort_alleles=True, **kwargs)
        if if_intersect:
            self.__intersect()

    def __intersect(self):
        """
        Intersect the Variant IDs in the LD matrix and the sumstats file.

        Raises
        ------
        ValueError
            If no common Variant IDs found between the LD matrix and the sumstats file.
        """
        intersec_index = self.map[
            self.map[ColName.SNPID].isin(self.sumstats[ColName.SNPID])
        ].index
        if len(intersec_index) == 0:
            raise ValueError(
                "No common Variant IDs found between the LD matrix and the sumstats file."
            )
        elif len(intersec_index) <= 10:
            logger.warning(
                "Only a few common Variant IDs found between the LD matrix and the sumstats file(<= 10)."
            )
        self.map = self.map.loc[intersec_index]
        self.r = self.r[intersec_index, :][:, intersec_index]
        self.sumstats = self.sumstats.loc[
            self.sumstats[ColName.SNPID].isin(self.map[ColName.SNPID])
        ]
        self.sumstats = self.sumstats.reset_index(drop=True)
        logger.info(
            "Intersected the Variant IDs in the LD matrix and the sumstats file. "
            f"Number of common Variant IDs: {len(intersec_index)}"
        )
