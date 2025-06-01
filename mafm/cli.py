"""Console script for mafm."""

import json
import logging
import os
from enum import Enum
from multiprocessing import Pool

import numpy as np
import pandas as pd
import typer
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeRemainingColumn

from mafm import __version__
from mafm.locus import load_locus_set
from mafm.mafm import fine_map, pipeline
from mafm.meta import meta_loci
from mafm.qc import loci_qc

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)


class MetaMethod(str, Enum):
    """The method to perform meta-analysis."""

    meta_all = "meta_all"
    meta_by_population = "meta_by_population"
    no_meta = "no_meta"


class Strategy(str, Enum):
    """The strategy to perform fine-mapping."""

    single_input = "single_input"
    multi_input = "multi_input"
    post_hoc_combine = "post_hoc_combine"


class Tool(str, Enum):
    """The tool to perform fine-mapping."""

    abf = "abf"
    carma = "carma"
    finemap = "finemap"
    rsparsepro = "rsparsepro"
    susie = "susie"
    multisusie = "multisusie"
    susiex = "susiex"


@app.callback(invoke_without_command=True, no_args_is_help=True)
def main(
    version: bool = typer.Option(False, "--version", "-V", help="Show version."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show verbose info."),
):
    """MAFM: Multi-ancestry fine-mapping pipeline."""
    console = Console()
    console.rule("[bold blue]MAFM[/bold blue]")
    console.print(f"Version: {__version__}", justify="center")
    console.print("Author: Jianhua Wang", justify="center")
    console.print("Email: jianhua.mert@gmail.com", justify="center")
    if version:
        typer.echo(f"MAFM version: {__version__}")
        raise typer.Exit()
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.info("Verbose mode is on.")
    else:
        for name in [
            "MAFM",
            "FINEMAP",
            "RSparsePro",
            "COJO",
            "SuSiE",
            "MULTISUSIE",
            "SUSIEX",
            "CARMA",
            "ABF",
            "Locus",
            "LDMatrix",
            "QC",
            "Sumstats",
            "Utils",
        ]:
            logging.getLogger(name).setLevel(logging.INFO)
        # logging.getLogger().setLevel(logging.INFO)
        # from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

        # rpy2_logger.setLevel(logging.ERROR)


@app.command(
    name="meta",
    help="Meta-analysis of summary statistics and LD matrices.",
)
def run_meta(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
    meta_method: MetaMethod = typer.Option(MetaMethod.meta_all, "--meta-method", "-m", help="Meta-analysis method."),
):
    """Meta-analysis of summary statistics and LD matrices."""
    meta_loci(inputs, outdir, threads, meta_method)


@app.command(
    name="qc",
    help="Quality control of summary statistics and LD matrices.",
)
def run_qc(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads."),
):
    """Quality control of summary statistics and LD matrices."""
    loci_qc(inputs, outdir, threads)


@app.command(
    name="finemap",
    help="Perform fine-mapping on three strategies.",
)
def run_fine_map(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    strategy: Strategy = typer.Option(Strategy.single_input, "--strategy", "-s", help="Fine-mapping strategy."),
    tool: Tool = typer.Option(Tool.susie, "--tool", "-t", help="Fine-mapping tool."),
    max_causal: int = typer.Option(1, "--max-causal", "-c", help="Maximum number of causal SNPs."),
    set_L_by_cojo: bool = typer.Option(True, "--set-L-by-cojo", "-sl", help="Set L by COJO."),
    p_cutoff: float = typer.Option(5e-8, "--p-cutoff", "-pc", help="P-value cutoff for COJO."),
    collinear_cutoff: float = typer.Option(0.9, "--collinear-cutoff", "-cc", help="Collinearity cutoff for COJO."),
    window_size: int = typer.Option(10000000, "--window-size", "-ws", help="Window size for COJO."),
    maf_cutoff: float = typer.Option(0.01, "--maf-cutoff", "-mc", help="MAF cutoff for COJO."),
    diff_freq_cutoff: float = typer.Option(
        0.2, "--diff-freq-cutoff", "-dfc", help="Difference in frequency cutoff for COJO."
    ),
    coverage: float = typer.Option(0.95, "--coverage", "-cv", help="Coverage of the credible set."),
    combine_cred: str = typer.Option("union", "--combine-cred", "-cc", help="Method to combine credible sets."),
    combine_pip: str = typer.Option("max", "--combine-pip", "-cp", help="Method to combine PIPs."),
    jaccard_threshold: float = typer.Option(
        0.1, "--jaccard-threshold", "-j", help="Jaccard threshold for combining credible sets."
    ),
    # susie parameters
    max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations."),
    estimate_residual_variance: bool = typer.Option(
        False, "--estimate-residual-variance", "-er", help="Estimate residual variance."
    ),
    min_abs_corr: float = typer.Option(0.5, "--min-abs-corr", "-mc", help="Minimum absolute correlation."),
    convergence_tol: float = typer.Option(1e-3, "--convergence-tol", "-ct", help="Convergence tolerance."),
):
    """Perform fine-mapping on three strategies."""
    loci_info = pd.read_csv(inputs, sep="\t")
    # Create progress bar
    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
    )

    # Get total number of loci
    locus_groups = list(loci_info.groupby("locus_id"))
    total_loci = len(locus_groups)

    with progress:
        task = progress.add_task("[cyan]Fine-mapping loci...", total=total_loci)

        for locus_id, locus_info in locus_groups:
            locus_set = load_locus_set(locus_info)
            creds = fine_map(
                locus_set,
                strategy=strategy,
                tool=tool,
                max_causal=max_causal,
                set_L_by_cojo=set_L_by_cojo,
                p_cutoff=p_cutoff,
                collinear_cutoff=collinear_cutoff,
                window_size=window_size,
                maf_cutoff=maf_cutoff,
                diff_freq_cutoff=diff_freq_cutoff,
                coverage=coverage,
                combine_cred=combine_cred,
                combine_pip=combine_pip,
                jaccard_threshold=jaccard_threshold,
                # susie parameters
                max_iter=max_iter,
                estimate_residual_variance=estimate_residual_variance,
                min_abs_corr=min_abs_corr,
                convergence_tol=convergence_tol,
            )
            out_dir = f"{outdir}/{locus_id}"
            os.makedirs(out_dir, exist_ok=True)
            creds.pips.to_csv(f"{out_dir}/pips.txt", sep="\t", header=False, index=True)
            with open(f"{out_dir}/creds.json", "w") as f:
                json.dump(creds.to_dict(), f, indent=4)

            progress.advance(task)


@app.command(
    name="pipeline",
    help="Run whole fine-mapping pipeline on a list of loci.",
)
def run_pipeline(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
    meta_method: MetaMethod = typer.Option(MetaMethod.meta_all, "--meta-method", "-m", help="Meta-analysis method."),
    skip_qc: bool = typer.Option(False, "--skip-qc", "-q", help="Skip quality control."),
    strategy: Strategy = typer.Option(Strategy.single_input, "--strategy", "-s", help="Fine-mapping strategy."),
    tool: Tool = typer.Option(Tool.susie, "--tool", "-t", help="Fine-mapping tool."),
    max_causal: int = typer.Option(1, "--max-causal", "-c", help="Maximum number of causal SNPs."),
    coverage: float = typer.Option(0.95, "--coverage", "-cv", help="Coverage of the credible set."),
    combine_cred: str = typer.Option("union", "--combine-cred", "-cc", help="Method to combine credible sets."),
    combine_pip: str = typer.Option("max", "--combine-pip", "-cp", help="Method to combine PIPs."),
    jaccard_threshold: float = typer.Option(
        0.1, "--jaccard-threshold", "-j", help="Jaccard threshold for combining credible sets."
    ),
    # ABF parameters
    var_prior: float = typer.Option(
        0.2,
        "--var-prior",
        "-vp",
        help="Variance prior, by default 0.2, usually set to 0.15 for quantitative traits and 0.2 for binary traits.",
        rich_help_panel="ABF",
    ),
    # FINEMAP parameters
    n_iter: int = typer.Option(100000, "--n-iter", "-ni", help="Number of iterations.", rich_help_panel="FINEMAP"),
    n_threads: int = typer.Option(1, "--n-threads", "-nt", help="Number of threads.", rich_help_panel="FINEMAP"),
    # susie parameters
    max_iter: int = typer.Option(
        100, "--max-iter", "-i", help="Maximum number of iterations.", rich_help_panel="SuSie"
    ),
    estimate_residual_variance: bool = typer.Option(
        False, "--estimate-residual-variance", "-er", help="Estimate residual variance.", rich_help_panel="SuSie"
    ),
    min_abs_corr: float = typer.Option(
        0.5, "--min-abs-corr", "-mc", help="Minimum absolute correlation.", rich_help_panel="SuSie"
    ),
    convergence_tol: float = typer.Option(
        1e-3, "--convergence-tol", "-ct", help="Convergence tolerance.", rich_help_panel="SuSie"
    ),
    # RSparsePro parameters
    eps: float = typer.Option(1e-5, "--eps", "-e", help="Convergence criterion.", rich_help_panel="RSparsePro"),
    ubound: int = typer.Option(
        100000, "--ubound", "-ub", help="Upper bound for convergence.", rich_help_panel="RSparsePro"
    ),
    cthres: float = typer.Option(0.7, "--cthres", "-ct", help="Threshold for coverage.", rich_help_panel="RSparsePro"),
    eincre: float = typer.Option(
        1.5, "--eincre", "-ei", help="Adjustment for error parameter.", rich_help_panel="RSparsePro"
    ),
    minldthres: float = typer.Option(
        0.7, "--minldthres", "-mlt", help="Threshold for minimum LD within effect groups.", rich_help_panel="RSparsePro"
    ),
    maxldthres: float = typer.Option(
        0.2, "--maxldthres", "-mlt", help="Threshold for maximum LD across effect groups.", rich_help_panel="RSparsePro"
    ),
    varemax: float = typer.Option(
        100.0, "--varemax", "-vm", help="Maximum error parameter.", rich_help_panel="RSparsePro"
    ),
    varemin: float = typer.Option(
        1e-3, "--varemin", "-vm", help="Minimum error parameter.", rich_help_panel="RSparsePro"
    ),
    # CARMA parameters
    effect_size_prior: str = typer.Option(
        "Spike-slab",
        "--effect-size-prior",
        "-es",
        help="Prior distribution for effect sizes ('Cauchy' or 'Spike-slab'), by default Spike-slab.",
        rich_help_panel="CARMA",
    ),
    input_alpha: float = typer.Option(
        0.0, "--input-alpha", "-ia", help="Elastic net mixing parameter (0 ≤ input_alpha ≤ 1).", rich_help_panel="CARMA"
    ),
    y_var: float = typer.Option(
        1.0, "--y-var", "-yv", help="Variance of the summary statistics.", rich_help_panel="CARMA"
    ),
    bf_threshold: float = typer.Option(
        10.0, "--bf-threshold", "-bf", help="Bayes factor threshold for credible models.", rich_help_panel="CARMA"
    ),
    outlier_bf_threshold: float = typer.Option(
        1 / 3.2,
        "--outlier-bf-threshold",
        "-obf",
        help="Bayes factor threshold for outlier detection.",
        rich_help_panel="CARMA",
    ),
    max_model_dim: int = typer.Option(
        200000, "--max-model-dim", "-mmd", help="Maximum number of top candidate models.", rich_help_panel="CARMA"
    ),
    all_inner_iter: int = typer.Option(
        10,
        "--all-inner-iter",
        "-aie",
        help="Maximum iterations for Shotgun algorithm within EM.",
        rich_help_panel="CARMA",
    ),
    all_iter: int = typer.Option(
        3, "--all-iter", "-ai", help="Maximum iterations for EM algorithm.", rich_help_panel="CARMA"
    ),
    tau: float = typer.Option(
        0.04, "--tau", "-t", help="Prior precision parameter of effect size.", rich_help_panel="CARMA"
    ),
    epsilon_threshold: float = typer.Option(
        1e-5, "--epsilon-threshold", "-et", help="Convergence threshold for Bayes factors.", rich_help_panel="CARMA"
    ),
    em_dist: str = typer.Option(
        "logistic", "--em-dist", "-ed", help="Distribution for modeling prior probability.", rich_help_panel="CARMA"
    ),
    # SuSiEx parameters
    # pval_thresh: float = typer.Option(1e-5, "--pval-thresh", "-pt", help="P-value threshold for SuSiEx.", rich_help_panel="SuSiEx"),
    # maf_thresh: float = typer.Option(0.005, "--maf-thresh", "-mt", help="MAF threshold for SuSiEx.", rich_help_panel="SuSiEx"),
    mult_step: bool = typer.Option(
        False, "--mult-step", "-ms", help="Whether to use multiple steps in SuSiEx.", rich_help_panel="SuSiEx"
    ),
    keep_ambig: bool = typer.Option(
        True, "--keep-ambig", "-ka", help="Whether to keep ambiguous SNPs in SuSiEx.", rich_help_panel="SuSiEx"
    ),
    # n_threads: int = typer.Option(1, "--n-threads", "-nt", help="Number of threads.", rich_help_panel="SuSiEx"),
    min_purity: float = typer.Option(
        0.5, "--min-purity", "-mp", help="Minimum purity for SuSiEx.", rich_help_panel="SuSiEx"
    ),
    # max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations.", rich_help_panel="SuSiEx"),
    tol: float = typer.Option(1e-3, "--tol", "-t", help="Convergence tolerance.", rich_help_panel="SuSiEx"),
    # MULTISUSIE parameters
    rho: float = typer.Option(
        0.75, "--rho", "-r", help="The prior correlation between causal variants.", rich_help_panel="MULTISUSIE"
    ),
    scaled_prior_variance: float = typer.Option(
        0.2, "--scaled-prior-variance", "-spv", help="The scaled prior variance.", rich_help_panel="MULTISUSIE"
    ),
    standardize: bool = typer.Option(
        False,
        "--standardize",
        "-s",
        help="Whether to standardize the summary statistics.",
        rich_help_panel="MULTISUSIE",
    ),
    pop_spec_standardization: bool = typer.Option(
        True,
        "--pop-spec-standardization",
        "-pss",
        help="Whether to use population-specific standardization.",
        rich_help_panel="MULTISUSIE",
    ),
    # estimate_residual_variance: bool = typer.Option(True, "--estimate-residual-variance", "-er", help="Estimate residual variance.", rich_help_panel="MULTISUSIE"),
    estimate_prior_variance: bool = typer.Option(
        True, "--estimate-prior-variance", "-epv", help="Estimate prior variance.", rich_help_panel="MULTISUSIE"
    ),
    estimate_prior_method: str = typer.Option(
        "early_EM",
        "--estimate-prior-method",
        "-epm",
        help="Method to estimate prior variance.",
        rich_help_panel="MULTISUSIE",
    ),
    pop_spec_effect_priors: bool = typer.Option(
        True,
        "--pop-spec-effect-priors",
        "-pesp",
        help="Whether to use population-specific effect priors.",
        rich_help_panel="MULTISUSIE",
    ),
    iter_before_zeroing_effects: int = typer.Option(
        5,
        "--iter-before-zeroing-effects",
        "-ibe",
        help="Number of iterations before zeroing out effects.",
        rich_help_panel="MULTISUSIE",
    ),
    prior_tol: float = typer.Option(
        1e-9, "--prior-tol", "-pt", help="Tolerance for prior variance.", rich_help_panel="MULTISUSIE"
    ),
    # min_abs_corr: float = typer.Option(0, "--min-abs-corr", "-mc", help="Minimum absolute correlation.", rich_help_panel="MULTISUSIE"),
    # max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations.", rich_help_panel="MULTISUSIE"),
    # tol: float = typer.Option(1e-3, "--tol", "-t", help="Convergence tolerance.", rich_help_panel="MULTISUSIE"),
):
    """Run whole fine-mapping pipeline on a list of loci."""
    loci_info = pd.read_csv(inputs, sep="\t")
    for locus_id, locus_info in loci_info.groupby("locus_id"):
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)
        pipeline(
            locus_info,
            outdir=out_dir,
            meta_method=meta_method,
            skip_qc=skip_qc,
            strategy=strategy,
            tool=tool,
            max_causal=max_causal,
            coverage=coverage,
            combine_cred=combine_cred,
            combine_pip=combine_pip,
            jaccard_threshold=jaccard_threshold,
            # susie parameters
            max_iter=max_iter,
            estimate_residual_variance=estimate_residual_variance,
            min_abs_corr=min_abs_corr,
            convergence_tol=convergence_tol,
            # ABF parameters
            var_prior=var_prior,
            # FINEMAP parameters
            n_iter=n_iter,
            n_threads=n_threads,
            # RSparsePro parameters
            eps=eps,
            ubound=ubound,
            cthres=cthres,
            eincre=eincre,
            minldthres=minldthres,
            maxldthres=maxldthres,
            varemax=varemax,
            varemin=varemin,
            # CARMA parameters
            effect_size_prior=effect_size_prior,
            input_alpha=input_alpha,
            y_var=y_var,
            bf_threshold=bf_threshold,
            outlier_bf_threshold=outlier_bf_threshold,
            max_model_dim=max_model_dim,
            all_inner_iter=all_inner_iter,
            all_iter=all_iter,
            tau=tau,
            epsilon_threshold=epsilon_threshold,
            em_dist=em_dist,
            # SuSiEx parameters
            mult_step=mult_step,
            keep_ambig=keep_ambig,
            min_purity=min_purity,
            tol=tol,
            # MULTISUSIE parameters
            rho=rho,
            scaled_prior_variance=scaled_prior_variance,
            standardize=standardize,
            pop_spec_standardization=pop_spec_standardization,
            estimate_prior_variance=estimate_prior_variance,
            estimate_prior_method=estimate_prior_method,
            pop_spec_effect_priors=pop_spec_effect_priors,
            iter_before_zeroing_effects=iter_before_zeroing_effects,
            prior_tol=prior_tol,
        )


if __name__ == "__main__":
    app(main)
