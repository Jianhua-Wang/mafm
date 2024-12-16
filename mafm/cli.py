"""Console script for mafm."""

import json
import logging
import os
from enum import Enum

import numpy as np
import pandas as pd
import typer
from rich.console import Console

from mafm import __version__
from mafm.locus import load_locus_set
from mafm.mafm import fine_map, pipeline
from mafm.meta import meta
from mafm.qc import locus_qc

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
        for name in ["MAFM", "FINEMAP", "RSparsePro", "SUSIE", "MULTISUSIE", "SUSIEX", "CARMA", "ABF", "Locus"]:
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
    meta_method: MetaMethod = typer.Option(MetaMethod.meta_all, "--meta-method", "-m", help="Meta-analysis method."),
):
    """Meta-analysis of summary statistics and LD matrices."""
    loci_info = pd.read_csv(inputs, sep="\t")
    for locus_id, locus_info in loci_info.groupby("locus_id"):
        locus_set = load_locus_set(locus_info)
        locus_set = meta(locus_set, meta_method)
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)
        for locus in locus_set.loci:
            out_prefix = f"{out_dir}/{locus.prefix}"
            locus.sumstats.to_csv(f"{out_prefix}.sumstat", sep="\t", index=False)
            np.savez_compressed(f"{out_prefix}.ld.npz", ld=locus.ld.r.astype(np.float16))
            locus.ld.map.to_csv(f"{out_prefix}.ldmap", sep="\t", index=False)


@app.command(
    name="qc",
    help="Quality control of summary statistics and LD matrices.",
)
def run_qc(
    inputs: str = typer.Argument(..., help="Input files."),
    outdir: str = typer.Argument(..., help="Output directory."),
):
    """Quality control of summary statistics and LD matrices."""
    loci_info = pd.read_csv(inputs, sep="\t")
    for locus_id, locus_info in loci_info.groupby("locus_id"):
        locus_set = load_locus_set(locus_info)
        qc_metrics = locus_qc(locus_set)
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)
        for k, v in qc_metrics.items():
            v.to_csv(f"{out_dir}/{k}.txt", sep="\t", index=False)


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
    for locus_id, locus_info in loci_info.groupby("locus_id"):
        locus_set = load_locus_set(locus_info)
        creds = fine_map(
            locus_set,
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
        )
        out_dir = f"{outdir}/{locus_id}"
        os.makedirs(out_dir, exist_ok=True)
        creds.pips.to_csv(f"{out_dir}/pips.txt", sep="\t", header=False, index=True)
        with open(f"{out_dir}/creds.json", "w") as f:
            json.dump(creds.to_dict(), f, indent=4)


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
    # susie parameters
    max_iter: int = typer.Option(100, "--max-iter", "-i", help="Maximum number of iterations."),
    estimate_residual_variance: bool = typer.Option(
        False, "--estimate-residual-variance", "-er", help="Estimate residual variance."
    ),
    min_abs_corr: float = typer.Option(0.5, "--min-abs-corr", "-mc", help="Minimum absolute correlation."),
    convergence_tol: float = typer.Option(1e-3, "--convergence-tol", "-ct", help="Convergence tolerance."),
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
            max_iter=max_iter,
            estimate_residual_variance=estimate_residual_variance,
            min_abs_corr=min_abs_corr,
            convergence_tol=convergence_tol,
        )


if __name__ == "__main__":
    app(main)
