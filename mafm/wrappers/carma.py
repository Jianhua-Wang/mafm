"""Wrapper for CARMA fine-mapping."""

import io
import json
import logging
from contextlib import redirect_stdout
from typing import Optional

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
from mafm.locus import Locus, intersect_sumstat_ld
from mafm.utils import check_r_package, io_in_tempdir

logger = logging.getLogger("CARMA")


@io_in_tempdir("./tmp/CARMA")
def run_carma(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    effect_size_prior: str = "Spike-slab",
    input_alpha: float = 0.0,
    y_var: float = 1.0,
    bf_threshold: float = 10.0,
    outlier_bf_threshold: float = 1 / 3.2,
    outlier_switch: bool = True,
    max_model_dim: int = 200000,
    all_inner_iter: int = 10,
    all_iter: int = 3,
    tau: float = 0.04,
    epsilon_threshold: float = 1e-5,
    printing_log: bool = False,
    em_dist: str = "logistic",
    temp_dir: Optional[str] = None,
) -> CredibleSet:
    """Run CARMA fine-mapping using R through rpy2.

    Parameters
    ----------
    locus : Locus
        Locus object containing summary statistics and LD matrix
    max_causal : int, optional
        Maximum number of causal variants assumed per locus, by default 1
    coverage : float, optional
        Coverage probability for credible sets, by default 0.95
    effect_size_prior : str, optional
        Prior distribution for effect sizes ('Cauchy' or 'Spike-slab'), by default "Spike-slab"
    input_alpha : float, optional
        Elastic net mixing parameter (0 ≤ input_alpha ≤ 1), by default 0.0
    y_var : float, optional
        Variance of the summary statistics, by default 1.0
    bf_threshold : float, optional
        Bayes factor threshold for credible models, by default 10.0
    outlier_bf_threshold : float, optional
        Bayes factor threshold for outlier detection, by default 1/3.2
    outlier_switch : bool, optional
        Whether to perform outlier detection, by default True
    max_model_dim : int, optional
        Maximum number of top candidate models, by default 200000
    all_inner_iter : int, optional
        Maximum iterations for Shotgun algorithm within EM, by default 10
    all_iter : int, optional
        Maximum iterations for EM algorithm, by default 3
    tau : float, optional
        Prior precision parameter of effect size, by default 0.04
    epsilon_threshold : float, optional
        Convergence threshold for Bayes factors, by default 1e-5
    printing_log : bool, optional
        Whether to print CARMA running log, by default False
    em_dist : str, optional
        Distribution for modeling prior probability ('logistic'), by default "logistic"
    temp_dir : Optional[str], optional
        Temporary directory, by default None.

    Returns
    -------
    CredibleSet
        Credible set

    Notes
    -----
    The function interfaces with the R package CARMA through rpy2 to perform
    Bayesian fine-mapping of GWAS loci.
    """
    if not check_r_package("CARMA"):
        raise RuntimeError("CARMA is not installed or R version is not supported.")
    if not locus.is_matched:
        logger.warning("The sumstat and LD are not matched, will match them in same order.")
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running CARMA on {locus}")

    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "effect_size_prior": effect_size_prior,
        "input_alpha": input_alpha,
        "y_var": y_var,
        "bf_threshold": bf_threshold,
        "outlier_bf_threshold": outlier_bf_threshold,
        "outlier_switch": outlier_switch,
        "max_model_dim": max_model_dim,
        "all_inner_iter": all_inner_iter,
        "all_iter": all_iter,
        "tau": tau,
        "epsilon_threshold": epsilon_threshold,
        "em_dist": em_dist,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")
    sumstats = locus.sumstats.copy()
    ld = locus.ld.r.copy()
    sumstats[ColName.Z] = sumstats[ColName.BETA] / sumstats[ColName.SE]

    # Import required R packages
    carma = importr("CARMA")

    # Create R lists for input
    z_list = ro.ListVector({"1": ro.FloatVector(sumstats[ColName.Z].values)})
    ld_matrix = ro.r.matrix(ro.FloatVector(ld.flatten()), nrow=ld.shape[0], ncol=ld.shape[1])  # type: ignore
    ld_list = ro.ListVector({"1": ld_matrix})
    lambda_list = ro.ListVector({"1": ro.FloatVector([1.0])})

    # Run CARMA with all parameters
    f = io.StringIO()
    with redirect_stdout(f):
        carma_results = carma.CARMA(
            z_list,
            ld_list,
            lambda_list=lambda_list,
            effect_size_prior=effect_size_prior,
            input_alpha=input_alpha,
            y_var=y_var,
            rho_index=coverage,
            BF_index=bf_threshold,
            outlier_BF_index=outlier_bf_threshold,
            outlier_switch=outlier_switch,
            num_causal=max_causal,
            Max_Model_Dim=max_model_dim,
            all_inner_iter=all_inner_iter,
            all_iter=all_iter,
            tau=tau,
            epsilon_threshold=epsilon_threshold,
            printing_log=printing_log,
            EM_dist=em_dist,
            output_labels=temp_dir,
        )
    logger.debug(f.getvalue())

    # Extract PIPs
    pips = np.array(carma_results[0].rx2("PIPs"))

    # Extract credible sets
    cs = np.zeros(len(sumstats))
    credible_sets = carma_results[0].rx2("Credible set")[1]
    if credible_sets:
        for i, cs_indices in enumerate(credible_sets, 1):
            cs[np.array(cs_indices, dtype=int) - 1] = i  # R uses 1-based indexing

    # Add results to summary statistics
    result_df = sumstats.copy()
    result_df["PIP"] = pips
    result_df["CS"] = cs.astype(int)

    pips = pd.Series(data=result_df["PIP"].to_numpy(), index=result_df["SNPID"].to_numpy())

    cs_snps = []
    lead_snps = []
    for cs_i, sub_df in result_df.groupby("CS"):
        if cs_i == 0:
            continue
        cs_snps.append(sub_df["SNPID"].values.tolist())
        lead_snps.append(pips[pips.index.isin(sub_df["SNPID"].values)].idxmax())
    return CredibleSet(
        tool=Method.CARMA,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=[len(i) for i in cs_snps],
        pips=pips,
        parameters=parameters,
    )
