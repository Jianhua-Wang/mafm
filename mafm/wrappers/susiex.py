"""Wrapper for SuSiEx."""

import json
import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
from mafm.locus import Locus, LocusSet, intersect_sumstat_ld
from mafm.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("SuSiEx")


@io_in_tempdir("./tmp/SuSiEx")
def run_susiex(
    locus_set: LocusSet,
    max_causal: int = 1,
    coverage: float = 0.95,
    pval_thresh: float = 1e-5,
    maf_thresh: float = 0.005,
    mult_step: bool = False,
    keep_ambig: bool = True,
    n_threads: int = 1,
    min_purity: float = 0.5,
    max_iter: int = 100,
    tol: float = 1e-3,
    temp_dir: Optional[str] = None,
):
    """
    Run SuSiEx on a LocusSet.

    Parameters
    ----------
    locus_set : LocusSet
        The LocusSet to run SuSiEx on.
    max_causal : int, optional
        The maximum number of causal SNPs, default is 1.
    coverage : float, optional
        The coverage of the credible set, default is 0.95.
    pval_thresh : float, optional
        The p-value threshold for SuSiEx, default is 1e-5.
    maf_thresh : float, optional
        The MAF threshold for SuSiEx, default is 0.005.
    mult_step : bool, optional
        Whether to use multiple steps in SuSiEx, default is False.
    keep_ambig : bool, optional
        Whether to keep ambiguous SNPs in SuSiEx, default is True.
    n_threads : int, optional
        The number of threads to use in SuSiEx, default is 1.
    min_purity : float, optional
        The minimum purity for SuSiEx, default is 0.5.
    max_iter : int, optional
        The maximum number of iterations for SuSiEx, default is 100.
    tol : float, optional
        The tolerance for SuSiEx, default is 1e-3.
    temp_dir : str, optional
        The temporary directory to use in SuSiEx, default is None.

    Returns
    -------
    CredibleSet
        The credible set.
    """
    logger.info(f"Running SuSiEx on {locus_set}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "pval_thresh": pval_thresh,
        "maf_thresh": maf_thresh,
        "mult_step": mult_step,
        "keep_ambig": keep_ambig,
        "n_threads": n_threads,
        "min_purity": min_purity,
        "max_iter": max_iter,
        "tol": tol,
    }
    logger.info(f"Parameters: {parameters}")

    input_prefix_list = []
    for locus in locus_set.loci:
        locus = intersect_sumstat_ld(locus)
        input_prefix = f"{temp_dir}/{locus.popu}.{locus.cohort}"
        logger.debug(f"Writing {input_prefix}.sumstats")
        locus.sumstats.to_csv(f"{input_prefix}.sumstats", sep="\t", index=False)
        ldmap = locus.ld.map.copy()
        ldmap["cm"] = 0
        logger.debug(f"Writing {input_prefix}_ref.bim")
        ldmap[[ColName.CHR, ColName.SNPID, "cm", ColName.BP, ColName.A1, ColName.A2]].to_csv(
            f"{input_prefix}_ref.bim", sep="\t", index=False, header=False
        )
        ldmap["MAF"] = ldmap["AF2"].apply(lambda x: min(x, 1 - x))
        ldmap["NCHROBS"] = locus.sample_size * 2
        ldmap.rename(columns={"SNPID": "SNP"}, inplace=True)
        logger.debug(f"Writing {input_prefix}_frq.frq")
        ldmap[["CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"]].to_csv(
            f"{input_prefix}_frq.frq", sep="\t", index=False
        )
        logger.debug(f"Writing {input_prefix}.ld.bin")
        ld = locus.ld.r**2
        ld.astype(np.float32).tofile(f"{input_prefix}.ld.bin")
        input_prefix_list.append(input_prefix)

    sst_file = ",".join([i + ".sumstats" for i in input_prefix_list])
    n_gwas = ",".join([str(locus.sample_size) for locus in locus_set.loci])
    ld_file = ",".join(input_prefix_list)
    chrom, start, end = locus_set.chrom, locus_set.start, locus_set.end
    cmd = [
        f"--sst_file={sst_file}",
        f"--n_gwas={n_gwas}",
        f"--ld_file={ld_file}",
        f"--out_dir={temp_dir}",
        f"--out_name=chr{chrom}_{start}_{end}",
        f"--level={coverage}",
        f"--pval_thresh={pval_thresh}",
        f"--maf={maf_thresh}",
        f"--chr={chrom}",
        f"--bp={start},{end}",
        f'--snp_col={",".join(["1"]*locus_set.n_loci)}',
        f'--chr_col={",".join(["2"]*locus_set.n_loci)}',
        f'--bp_col={",".join(["3"]*locus_set.n_loci)}',
        f'--a1_col={",".join(["5"]*locus_set.n_loci)}',
        f'--a2_col={",".join(["6"]*locus_set.n_loci)}',
        f'--eff_col={",".join(["9"]*locus_set.n_loci)}',
        f'--se_col={",".join(["10"]*locus_set.n_loci)}',
        f'--pval_col={",".join(["11"]*locus_set.n_loci)}',
        "--plink=../utilities/plink",
        f"--n_sig={max_causal}",
        f"--mult-step={mult_step}",
        f"--keep-ambig={keep_ambig}",
        f"--threads={n_threads}",
        f"--min_purity={min_purity}",
        f"--max_iter={max_iter}",
        f"--tol={tol}",
    ]
    required_output_files = [f"{temp_dir}/chr{chrom}_{start}_{end}.snp", f"{temp_dir}/chr{chrom}_{start}_{end}.cs"]
    logger.info(f"Running SuSiEx with command: {' '.join(cmd)}.")
    tool_manager.run_tool("SuSiEx", cmd, f"{temp_dir}/run.log", required_output_files)

    pip_df = pd.read_csv(f"{temp_dir}/chr{chrom}_{start}_{end}.snp", sep="\t")
    cs_snp = []
    if len(pip_df.columns) == 2:
        logger.warning("No credible set found, please try other parameters.")
        pip = pd.Series(index=pip_df["SNP"].values.tolist())
    else:
        cs_df = pd.read_csv(f"{temp_dir}/chr{chrom}_{start}_{end}.cs", sep="\t")
        for _, sub_df in cs_df.groupby("CS_ID"):
            cs_snp.append(sub_df["SNP"].values.tolist())
        pip_cols = [col for col in pip_df.columns if col.startswith("PIP")]
        pip_df = pip_df[pip_cols + ["SNP"]].copy()
        pip_df.set_index("SNP", inplace=True)
        pip_df["PIP"] = pip_df[pip_cols].max(axis=1)
        pip = pd.Series(index=pip_df.index.values.tolist(), data=pip_df["PIP"].values.tolist())

    logger.info(f"Fished SuSiEx on {locus_set}")
    logger.info(f"N of credible set: {len(cs_snp)}")
    logger.info(f"Credible set size: {[len(i) for i in cs_snp]}")
    return CredibleSet(
        tool=Method.SUSIEX,
        n_cs=len(cs_snp),
        coverage=coverage,
        lead_snps=[str(pip[pip.index.isin(i)].idxmax()) for i in cs_snp],
        snps=cs_snp,
        cs_sizes=[len(i) for i in cs_snp],
        pips=pip,
        parameters=parameters,
    )
