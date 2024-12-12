"""Main module."""

import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import toml

from mafm.credibleset import CredibleSet

logger = logging.getLogger("MAFM")


class FmOutput:
    """Class for managing multiple fine-mapping results.

    Attributes
    ----------
    sumstat : pd.DataFrame
        Summary statistics dataframe.
    pips : pd.DataFrame
        Posterior inclusion probabilities dataframe.
    credible_sets : List[CredibleSet]
        List of credible sets.
    r : Optional[np.ndarray], optional
        Linkage disequilibrium matrix, by default None.
    map_df : Optional[pd.DataFrame], optional
        Mapping dataframe, by default None.
    n_res : int
        Number of credible sets.
    """

    def __init__(
        self,
        sumstat: pd.DataFrame,
        credible_sets: Dict[str, CredibleSet],
        r: Optional[np.ndarray] = None,
        map_df: Optional[pd.DataFrame] = None,
    ):
        """
        Initialize FmOutput with summary statistics, PIPs, credible sets, and optional LD matrix and map.

        Parameters
        ----------
        sumstat : pd.DataFrame
            Summary statistics dataframe.
        credible_sets : List[CredibleSet]
            List of credible sets.
        r : Optional[np.ndarray], optional
            Linkage disequilibrium matrix, by default None.
        map_df : Optional[pd.DataFrame], optional
            Mapping dataframe, by default None.
        """
        self.sumstat = sumstat
        self.r = r
        self.map = map_df
        self.cs = credible_sets
        self.n_res = len(credible_sets)
        self.tools = list(credible_sets.keys())

        # Create tool-specific attributes dynamically
        pips_df = []
        for v in credible_sets.values():
            # setattr(self, k, v)
            pips_df.append(v.pips)
        self.pips = pd.concat(pips_df, axis=1)

    def __getattr__(self, name: str):
        """Allow access to credible sets as attributes."""
        if name in self.cs:
            return self.cs[name]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    def save_results(self, prefix: str):
        """Save all results to files.

        Parameters
        ----------
        prefix : str
            Prefix for the output files.
        """
        outdir = Path(prefix)
        outdir.parent.mkdir(parents=True, exist_ok=True)

        # Save PIPs
        pips_file = f"{prefix}_pips.txt"
        self.pips.to_csv(pips_file, sep="\t", index=True)

        # Save credible sets info to TOML
        cs_data = {cs.tool: cs.to_dict() for cs in self.cs.values()}
        toml_file = f"{prefix}_cs.toml"
        with open(toml_file, "w") as f:
            toml.dump({"credible_sets": cs_data}, f)

        # Save sumstats
        if self.sumstat is not None:
            sumstat_file = f"{prefix}.munged.sumstats"
            self.sumstat.to_csv(sumstat_file, sep="\t", index=False)

        # Save map
        if self.map is not None:
            map_file = f"{prefix}.munged.ldmap"
            self.map.to_csv(map_file, sep="\t", index=False)

        # Save LD matrix
        if self.r is not None:
            ld_file = f"{prefix}.ld.npz"
            np.savez(ld_file, ld=self.r)

    @classmethod
    def load_results(cls, prefix: str) -> "FmOutput":
        """Load results from files.

        Parameters
        ----------
        prefix : str
            Prefix for the input files.

        Returns
        -------
        FmOutput
            An instance of FmOutput with loaded data.

        Raises
        ------
        FileNotFoundError
            If the sumstats file is not found.
        """
        # Load PIPs
        pips_file = f"{prefix}_pips.txt"
        pips = pd.read_csv(pips_file, sep="\t", index_col=0)

        # Load credible sets from TOML
        toml_file = f"{prefix}_cs.toml"
        with open(toml_file, "r") as f:
            cs_data = toml.load(f)["credible_sets"]

        # Create CredibleSet objects
        credible_sets = {}
        for cs_dict in cs_data.values():
            tool = cs_dict["tool"]
            tool_pips = pips[tool] if tool in pips.columns else pd.Series()
            cs = CredibleSet.from_dict(cs_dict, tool_pips)
            credible_sets[tool] = cs

        # Load optional files
        sumstat = None
        map_df = None
        r = None

        sumstat_file = f"{prefix}.munged.sumstats"
        if Path(sumstat_file).exists():
            sumstat = pd.read_csv(sumstat_file, sep="\t")
        else:
            raise FileNotFoundError(f"Sumstats file not found: {sumstat_file}")

        map_file = f"{prefix}.munged.ldmap"
        if Path(map_file).exists():
            map_df = pd.read_csv(map_file, sep="\t")

        ld_file = f"{prefix}.ld.npz"
        if Path(ld_file).exists():
            r = np.load(ld_file)["ld"]

        return cls(sumstat=sumstat, credible_sets=credible_sets, r=r, map_df=map_df)

    def summary(self) -> Dict:
        """
        Generate summary of all fine-mapping results.

        Returns
        -------
        Dict
            Summary dictionary with tool names as keys and summary statistics as values.
        """
        return {
            cs.tool: {"n_cs": cs.n_cs, "coverage": cs.coverage, "total_snps": sum(cs.cs_sizes)}
            for cs in self.cs.values()
        }
