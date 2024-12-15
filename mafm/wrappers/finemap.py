"""Wrapper for FINEMAP."""

import json
import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet, combine_creds
from mafm.locus import Locus, intersect_sumstat_ld
from mafm.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("FINEMAP")


@io_in_tempdir("./tmp/FINEMAP")
def run_finemap(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    n_iter: int = 100000,
    n_threads: int = 1,
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """
    Run FINEMAP with shotgun stochastic search.

    Parameters
    ----------
    locus : Locus
        Locus object.
    max_causal : int, optional
        Maximum number of causal variants, by default 1, only support 1.
    coverage : float, optional
        Coverage of the credible set, by default 0.95.
    n_iter : int, optional
        Number of iterations, by default 100000.
    n_threads : int, optional
        Number of threads, by default 1.
    temp_dir : Optional[str], optional
        Temporary directory, by default None.

    Returns
    -------
    CredibleSet
        Credible set.

    """
    logger.info(f"Running FINEMAP on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "n_iter": n_iter,
        "n_threads": n_threads,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    if not locus.is_matched:
        logger.warning("The sumstat and LD are not matched, will match them in same order.")
        locus = intersect_sumstat_ld(locus)
    # write z file
    if ColName.MAF not in locus.sumstats.columns:
        raise ValueError(f"{ColName.MAF} is required for FINEMAP.")
    finemap_input = locus.sumstats.copy()
    finemap_input[ColName.MAF] = finemap_input[ColName.MAF].replace(0, 0.00001)
    finemap_input = finemap_input[
        [
            ColName.SNPID,
            ColName.CHR,
            ColName.BP,
            ColName.EA,
            ColName.NEA,
            ColName.MAF,
            ColName.BETA,
            ColName.SE,
        ]
    ]
    finemap_input.rename(
        columns={
            ColName.SNPID: "rsid",
            ColName.CHR: "chromosome",
            ColName.BP: "position",
            ColName.MAF: "maf",
            ColName.BETA: "beta",
            ColName.SE: "se",
            ColName.EA: "allele1",
            ColName.NEA: "allele2",
        },
        inplace=True,
    )
    logger.info(f"Writing FINEMAP input to {temp_dir}/finemap.z")
    finemap_input.to_csv(f"{temp_dir}/finemap.z", sep=" ", index=False, float_format="%0.5f")

    # write ld file
    logger.info(f"Writing FINEMAP LD file to {temp_dir}/finemap.ld")
    np.savetxt(f"{temp_dir}/finemap.ld", locus.ld.r, delimiter=" ", fmt="%0.4f")
    # TODO: write ld file only once for multiple tools
    # TODO: use BCOR file for LD

    # write master file
    logger.info(f"Writing FINEMAP master file to {temp_dir}/finemap.master")
    with open(f"{temp_dir}/finemap.master", "w") as f:
        master_content = [
            f"{temp_dir}/finemap.z",
            f"{temp_dir}/finemap.ld",
            f"{temp_dir}/finemap.snp",
            f"{temp_dir}/finemap.config",
            f"{temp_dir}/finemap.cred",
            f"{temp_dir}/finemap.log",
            str(locus.sample_size),
        ]
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(";".join(master_content))

    # run finemap
    cmd = [
        "--sss",
        "--in-files",
        f"{temp_dir}/finemap.master",
        "--n-causal-snps",
        str(max_causal),
        "--n-iter",
        str(n_iter),
        "--n-threads",
        str(n_threads),
        "--prob-cred-set",
        str(coverage),
    ]
    required_output_files = [f"{temp_dir}/finemap.snp", f"{temp_dir}/finemap.config"]
    logger.info(f"Running FINEMAP with command: {' '.join(cmd)}.")
    tool_manager.run_tool("finemap", cmd, f"{temp_dir}/run.log", required_output_files)

    # get PIPs
    if os.path.getsize(f"{temp_dir}/finemap.snp") == 0:
        logger.warning("FINEMAP output is empty.")
        pip = pd.Series(index=finemap_input["rsid"].values.tolist())
    else:
        finemap_res = pd.read_csv(f"{temp_dir}/finemap.snp", sep=" ", usecols=["rsid", "prob"])
        finemap_res = pd.Series(finemap_res["prob"].values, index=finemap_res["rsid"].values)  # type: ignore
        pip = finemap_res

    # get credible set
    if os.path.getsize(f"{temp_dir}/finemap.config") == 0:
        logger.warning("FINEMAP output is empty.")
        no_cred = True
    else:
        finemap_config = pd.read_csv(f"{temp_dir}/finemap.config", sep=" ", usecols=["config", "prob"])
        finemap_config = finemap_config.sort_values("prob", ascending=False)
        # TODO: limit the number of causal SNPs when there are too many SNPs with very low PIPs
        finemap_config = finemap_config[finemap_config["prob"].shift().fillna(0).cumsum() <= coverage]
        cs_snps = list(set(finemap_config["config"].str.cat(sep=",").split(",")))
        lead_snps = str(
            locus.sumstats.loc[
                locus.sumstats[locus.sumstats[ColName.SNPID].isin(cs_snps)][ColName.P].idxmin(), ColName.SNPID
            ]
        )
        no_cred = False

    # output
    logger.info(f"Fished FINEMAP on {locus}")
    logger.warning('FINEMAP outputs configuration file, not credible set. Concatenate the configurations to one credible set.')
    logger.info("N of credible set: 1")
    logger.info(f"Credible set size: {len(cs_snps)}")
    return CredibleSet(
        tool=Method.FINEMAP,
        n_cs=1 if not no_cred else 0,
        coverage=coverage,
        lead_snps=[lead_snps] if not no_cred else [],
        snps=[cs_snps] if not no_cred else [],
        cs_sizes=[len(cs_snps)] if not no_cred else [],
        pips=pip,
        parameters={"max_causal": max_causal, "method": "sss", "n_iter": n_iter},
    )


def run_finemap_multi(
    inputs: list[Locus],
    max_causal: int = 1,
    coverage: float = 0.95,
    n_iter: int = 100000,
    n_threads: int = 1,
    combine_cred: str = "union",
    combine_pip: str = "max",
    jaccard_threshold: float = 0.1,
) -> CredibleSet:
    """
    Run FINEMAP for multiple datasets.

    Run FINEMAP for each dataset in the input list. Then combine the results.

    Parameters
    ----------
    inputs : list[FmInput]
        List of input data.
    max_causal : int, optional
        Maximum number of causal variants, by default 1, only support 1.
    coverage : float, optional
        Coverage, by default 0.95.
    n_iter : int, optional
        Number of iterations, by default 100000.
    n_threads : int, optional
        Number of threads, by default 1.
    combine_cred : str, optional
        Method to combine credible sets, by default "union".
    combine_pip : str, optional
        Method to combine PIPs, by default "max".
    jaccard_threshold : float, optional
        Jaccard index threshold for merging credible sets, by default 0.1.

    Returns
    -------
    CredibleSet
        Combined credible set.

    """
    logger.info("Running FINEMAP for multiple datasets.")
    logger.info(f"max_causal={max_causal}, coverage={coverage}, n_iter={n_iter}, n_threads={n_threads}")
    logger.info(f"combine_cred={combine_cred}, combine_pip={combine_pip}, jaccard_threshold={jaccard_threshold}")

    # run FINEMAP for each dataset
    cs_list = []
    for input in inputs:
        cs = run_finemap(input, max_causal, coverage, n_iter, n_threads)
        cs_list.append(cs)

    # combine credible sets
    combined_cs = combine_creds(cs_list, combine_cred, combine_pip, jaccard_threshold)
    logger.info("Finished FINEMAP for multiple datasets.")
    return combined_cs
