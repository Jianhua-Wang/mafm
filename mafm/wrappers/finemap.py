"""Wrapper for FINEMAP."""

import logging
import os
import subprocess
from typing import List, Optional, Union

from mafm.constants import ColName
from mafm.mafm import FmInput
from mafm.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("FINEMAP")


def run_finemap(
    mfinput: FmInput,
    subprogram: str = "--sss",
    cond_pvalue: float = 5e-8,
    corr_config: float = 0.95,
    force_n_samples: bool = False,
    n_causal_snps: int = 5,
    n_configs_top: int = 50000,
    n_conv_sss: int = 100,
    n_iter: int = 100000,
    n_threads: int = 1,
    prior_k: bool = False,
    prior_k0: float = 0.0,
    prob_conv_sss_tol: float = 0.001,
    prob_cred_set: float = 0.95,
    pvalue_snps: float = 1.0,
    rsids: Optional[List[str]] = None,
    std_effects: bool = False,
) -> None:
    """
    Run FINEMAP with specified parameters.

    Parameters
    ----------
    executable_path : str
        Path to the FINEMAP executable.
    in_files : str
        Path to the master file.
    subprogram : str, optional
        Subprogram to run ('--cond', '--config', or '--sss') (default: '--sss').
    cond_pvalue : float, optional
        P-value threshold for declaring genome-wide significance (default: 5e-8).
    corr_config : float, optional
        Correlation threshold for causal configuration (default: 0.95).
    force_n_samples : bool, optional
        Allow correlations with different sample size than GWAS (default: False).
    n_causal_snps : int, optional
        Maximum number of allowed causal SNPs (default: 5).
    n_configs_top : int, optional
        Number of top causal configurations to save (default: 50000).
    n_conv_sss : int, optional
        Number of iterations for convergence in SSS (default: 100).
    n_iter : int, optional
        Maximum number of iterations for SSS (default: 100000).
    n_threads : int, optional
        Number of parallel threads (default: 1).
    prior_k : bool, optional
        Use prior probabilities for number of causal SNPs (default: False).
    prior_k0 : float, optional
        Prior probability of no causal SNP (default: 0.0).
    prob_conv_sss_tol : float, optional
        Tolerance for SSS convergence (default: 0.001).
    prob_cred_set : float, optional
        Probability for credible interval (default: 0.95).
    pvalue_snps : float, optional
        P-value threshold for including SNPs (default: 1.0).
    rsids : List[str], optional
        List of SNP identifiers (for --config subprogram).
    std_effects : bool, optional
        Print mean and std dev of posterior effect size distribution (default: False).

    Returns
    -------
    subprocess.CompletedProcess
        Result of the subprocess run.

    Examples
    --------
    >>> result = run_finemap(
    ...     executable_path="/path/to/finemap",
    ...     in_files="master_file.txt",
    ...     subprogram="--sss",
    ...     n_threads=4,
    ...     n_causal_snps=10
    ... )
    >>> print(f"FINEMAP exited with return code {result.returncode}")
    >>> print(f"FINEMAP output: {result.stdout}")
    """
    # command = [executable_path, subprogram, "--in-files", in_files]

    # Add other parameters
    # command.extend(
    #     [
    #         "--cond-pvalue",
    #         str(cond_pvalue),
    #         "--corr-config",
    #         str(corr_config),
    #         "--n-causal-snps",
    #         str(n_causal_snps),
    #         "--n-configs-top",
    #         str(n_configs_top),
    #         "--n-conv-sss",
    #         str(n_conv_sss),
    #         "--n-iter",
    #         str(n_iter),
    #         "--n-threads",
    #         str(n_threads),
    #         "--prior-k0",
    #         str(prior_k0),
    #         "--prob-conv-sss-tol",
    #         str(prob_conv_sss_tol),
    #         "--prob-cred-set",
    #         str(prob_cred_set),
    #         "--pvalue-snps",
    #         str(pvalue_snps),
    #     ]
    # )

    # if force_n_samples:
    #     command.append("--force-n-samples")

    # if prior_k:
    #     command.append("--prior-k")
    # if prior_snps:
    #     command.append("--prior-snps")
    # if isinstance(prior_std, list):
    #     command.extend(["--prior-std", ",".join(map(str, prior_std))])
    # else:
    #     command.extend(["--prior-std", str(prior_std)])
    # if rsids:
    #     command.extend(["--rsids", ",".join(rsids)])
    # if std_effects:
    #     command.append("--std-effects")

    # Run FINEMAP
    # result = subprocess.run(command, capture_output=True, text=True, check=True)
    # return result
    pass
