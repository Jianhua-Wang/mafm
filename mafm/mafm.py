"""Main module."""

import inspect
import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import toml

from mafm.credibleset import CredibleSet, combine_creds
from mafm.locus import LocusSet
from mafm.wrappers import run_abf, run_carma, run_finemap, run_multisusie, run_rsparsepro, run_susie, run_susiex

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


def fine_map(
    locus_set: LocusSet,
    strategy: str = "single_input",
    tool: str = "susie",
    max_causal: int = 1,
    combine_cred: str = "union",
    combine_pip: str = "max",
    jaccard_threshold: float = 0.1,
    **kwargs,
) -> CredibleSet:
    """
    Perform fine-mapping on a locus set.

    Parameters
    ----------
    locus_set : LocusSet
        Locus set to fine-mapping.
    strategy : str
        Fine-mapping strategy. Choose from ["single_input", "multi_input", "post_hoc_combine"]
        single_input: use fine-mapping tools which take a single locus as input, these tools are:
            abf, carma, finemap, rsparsepro, susie
        multi_input: use fine-mapping tools which take multiple loci as input, these tools are:
            multisusie, susiex
        post_hoc_combine: use fine-mapping tools which take single loci as input (see single_input options), and then combine the results,
    tool : str
        Fine-mapping tool. Choose from ["abf", "carma", "finemap", "rsparsepro", "susie", "multisusie", "susiex"]
    combine_cred : str, optional
        Method to combine credible sets, by default "union".
        Options: "union", "intersection", "cluster".
        "union":        Union of all credible sets to form a merged credible set.
        "intersection": Frist merge the credible sets from the same tool,
                        then take the intersection of all merged credible sets.
                        no credible set will be returned if no common SNPs found.
        "cluster":      Merge credible sets with Jaccard index > 0.1.
    combine_pip : str, optional
        Method to combine PIPs, by default "max".
        Options: "max", "min", "mean", "meta".
        "meta": PIP_meta = 1 - prod(1 - PIP_i), where i is the index of tools,
                PIP_i = 0 when the SNP is not in the credible set of the tool.
        "max":  Maximum PIP value for each SNP across all tools.
        "min":  Minimum PIP value for each SNP across all tools.
        "mean": Mean PIP value for each SNP across all tools.
    jaccard_threshold : float, optional
        Jaccard index threshold for the "cluster" method, by default 0.1.
    max_causal : int, optional
        Maximum number of causal variants, by default 1.
    """
    tool_func_dict = {
        "abf": run_abf,
        "carma": run_carma,
        "finemap": run_finemap,
        "rsparsepro": run_rsparsepro,
        "susie": run_susie,
        "multisusie": run_multisusie,
        "susiex": run_susiex,
    }
    inspect_dict = {
        "abf": set(inspect.signature(run_abf).parameters),
        "carma": set(inspect.signature(run_carma).parameters),
        "finemap": set(inspect.signature(run_finemap).parameters),
        "rsparsepro": set(inspect.signature(run_rsparsepro).parameters),
        "susie": set(inspect.signature(run_susie).parameters),
        "multisusie": set(inspect.signature(run_multisusie).parameters),
        "susiex": set(inspect.signature(run_susiex).parameters),
    }
    params_dict = {}
    for t, args in inspect_dict.items():
        params_dict[t] = {k: v for k, v in kwargs.items() if k in args}
    if strategy == "single_input":
        if locus_set.n_loci > 1:
            raise ValueError("Locus set must contain only one locus for single-input strategy")
        if tool in ["abf", "carma", "finemap", "rsparsepro", "susie"]:
            return tool_func_dict[tool](locus_set.loci[0], **params_dict[tool])
        else:
            raise ValueError(f"Tool {tool} not supported for single-input strategy")
    elif strategy == "multi_input":
        if locus_set.n_loci < 2:
            raise ValueError("Locus set must contain at least two loci for multi-input strategy")
        if tool in ["multisusie", "susiex"]:
            return tool_func_dict[tool](locus_set, **params_dict[tool])
        else:
            raise ValueError(f"Tool {tool} not supported for multi-input strategy")
    elif strategy == "post_hoc_combine":
        if locus_set.n_loci < 2:
            raise ValueError("Locus set must contain at least two loci for post-hoc combine strategy")
        if tool in ["abf", "carma", "finemap", "rsparsepro", "susie"]:
            all_creds = []
            for locus in locus_set.loci:
                creds = tool_func_dict[tool](locus, max_causal=max_causal, **params_dict[tool])
                all_creds.append(creds)
            return combine_creds(all_creds, combine_cred=combine_cred, combine_pip=combine_pip, jaccard_threshold=jaccard_threshold)
        else:
            raise ValueError(f"Tool {tool} not supported for post-hoc combine strategy")
    else:
        raise ValueError(f"Strategy {strategy} not supported")
