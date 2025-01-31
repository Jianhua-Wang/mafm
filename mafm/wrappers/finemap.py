"""Wrapper for FINEMAP."""

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
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
    # change maf to 0.00001 if maf is 0
    finemap_input.loc[finemap_input["maf"] <= 0.00001, "maf"] = 0.00001
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

    # get credible set
    cred_file_list = Path(f"{temp_dir}/").glob("finemap.cred*")
    cred_prob = {}
    pip = pd.Series(index=finemap_input["rsid"].values.tolist(), data=0.0)
    n_cs = 0
    cs_snps = []
    lead_snps = []
    cs_sizes = []
    for cred_file in cred_file_list:
        with open(cred_file, "r") as f:
            n_causal = int(cred_file.name[-1])
            first_line = f.readline()
            cred_prob[n_causal] = float(first_line.split()[-1])
    if len(cred_prob) == 0:
        logger.warning("FINEMAP output is empty.")
    else:
        # get the credible set with the highest posterior probability
        n_cs = max(cred_prob, key=lambda k: float(cred_prob[k]))
        logger.info(f"FINEMAP found {n_cs} causal SNPs with the post-prob {cred_prob[n_cs]}.")
        cred_set = pd.read_csv(f"{temp_dir}/finemap.cred{n_cs}", sep=" ", comment="#")
        for cred_idx in range(1, n_cs + 1):
            cred_df = cred_set[[f"cred{cred_idx}", f"prob{cred_idx}"]].copy()
            cred_df.rename(columns={f"cred{cred_idx}": "snp", f"prob{cred_idx}": "pip"}, inplace=True)
            cred_df.dropna(inplace=True)
            cs_snps.append(cred_df["snp"].values.tolist())
            cs_sizes.append(len(cred_df["snp"].values.tolist()))
            pip[cred_df["snp"].values.tolist()] = cred_df["pip"].values.tolist()
            lead_snps.append(
                str(
                    locus.sumstats.loc[
                        locus.sumstats[locus.sumstats[ColName.SNPID].isin(cred_df["snp"].values.tolist())][
                            ColName.P
                        ].idxmin(),
                        ColName.SNPID,
                    ]
                )
            )

    # output
    logger.info(f"Fished FINEMAP on {locus}")
    logger.info(f"N of credible set: {n_cs}")
    logger.info(f"Credible set size: {cs_sizes}")
    return CredibleSet(
        tool=Method.FINEMAP,
        n_cs=n_cs,
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=cs_sizes,
        pips=pip,
        parameters=parameters,
    )
