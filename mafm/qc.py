"""Quality control functions for MAFM data."""

import logging
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn.mixture import GaussianMixture

from mafm.constants import ColName
from mafm.Locus import Locus, intersect

logger = logging.getLogger("QC")


def get_eigen(ldmatrix: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Compute eigenvalues and eigenvectors of R.

    TODO: accelerate with joblib

    Parameters
    ----------
    R : np.ndarray
        A p by p symmetric, positive semidefinite correlation matrix.

    Returns
    -------
    dict
        Dictionary containing eigenvalues and eigenvectors.
    """
    ldmatrix = ldmatrix.astype(np.float32)
    eigvals, eigvecs = np.linalg.eigh(ldmatrix)
    return {"eigvals": eigvals, "eigvecs": eigvecs}


def estimate_s_rss(
    locus: Locus, r_tol: float = 1e-8, method: str = "null-mle", eigvens: Optional[Dict[str, np.ndarray]] = None
) -> float:
    """
    Estimate s in the susie_rss Model Using Regularized LD.

    This function estimates the parameter s, which provides information about the consistency between z-scores
    and the LD matrix. A larger s indicates a strong inconsistency between z-scores and the LD matrix.

    Parameters
    ----------
    locus : Locus
        Locus object.
    r_tol : float, default=1e-8
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.
    method : str, default="null-mle"
        Method to estimate s. Options are "null-mle", "null-partialmle", or "null-pseudomle".

    Returns
    -------
    float
        Estimated s value between 0 and 1 (or potentially > 1 for "null-partialmle").
    """
    # make sure the LD matrix and sumstats file are matched
    input_locus = locus.copy()
    input_locus = intersect(input_locus)
    z = (input_locus.sumstats[ColName.BETA] / input_locus.sumstats[ColName.SE]).to_numpy()
    n = input_locus.sample_size
    # Check and process input arguments z, R
    z = np.where(np.isnan(z), 0, z)
    if eigvens is not None:
        eigvals = eigvens["eigvals"]
        eigvecs = eigvens["eigvecs"]
    else:
        eigens = get_eigen(input_locus.ld.r)
        eigvals = eigens["eigvals"]
        eigvecs = eigens["eigvecs"]

    if np.any(eigvals < -r_tol):
        logger.warning("The LD matrix is not positive semidefinite. Negative eigenvalues are set to zero")
    eigvals[eigvals < r_tol] = 0

    if n <= 1:
        raise ValueError("n must be greater than 1")

    sigma2 = (n - 1) / (z**2 + n - 2)
    z = np.sqrt(sigma2) * z

    if method == "null-mle":

        def negloglikelihood(s, ztv, d):
            denom = (1 - s) * d + s
            term1 = 0.5 * np.sum(np.log(denom))
            term2 = 0.5 * np.sum((ztv / denom) * ztv)
            return term1 + term2

        ztv = eigvecs.T @ z
        result = minimize_scalar(
            negloglikelihood,
            bounds=(0, 1),
            method="bounded",
            args=(ztv, eigvals),
            options={"xatol": np.sqrt(np.finfo(float).eps)},
        )
        s = result.x  # type: ignore

    elif method == "null-partialmle":
        colspace = np.where(eigvals > 0)[0]
        if len(colspace) == len(z):
            s = 0
        else:
            znull = eigvecs[:, ~np.isin(np.arange(len(z)), colspace)].T @ z
            s = np.sum(znull**2) / len(znull)

    elif method == "null-pseudomle":

        def pseudolikelihood(s: float, z: np.ndarray, eigvals: np.ndarray, eigvecs: np.ndarray) -> float:
            precision = eigvecs @ (eigvecs.T / ((1 - s) * eigvals + s))
            postmean = np.zeros_like(z)
            postvar = np.zeros_like(z)
            for i in range(len(z)):
                postmean[i] = -(1 / precision[i, i]) * precision[i, :].dot(z) + z[i]
                postvar[i] = 1 / precision[i, i]
            return -np.sum(stats.norm.logpdf(z, loc=postmean, scale=np.sqrt(postvar)))

        result = minimize_scalar(
            pseudolikelihood,
            bounds=(0, 1),
            method="bounded",
            args=(z, eigvals, eigvecs),
        )
        s = result.x  # type: ignore

    else:
        raise ValueError("The method is not implemented")

    return s  # type: ignore


def kriging_rss(
    locus: Locus,
    r_tol: float = 1e-8,
    s: Optional[float] = None,
    eigvens: Optional[Dict[str, np.ndarray]] = None,
) -> pd.DataFrame:
    """
    Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue.

    Under the null, the rss model with regularized LD matrix is z|R,s ~ N(0, (1-s)R + s I)).
    We use a mixture of normals to model the conditional distribution of z_j given other z scores.

    Parameters
    ----------
    locus : Locus
        Locus object.
    r_tol : float = 1e-8
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.
    s : Optional[float] = None
        An estimated s from estimate_s_rss function.
    eigvens : Optional[Dict[str, np.ndarray]] = None
        A dictionary containing eigenvalues and eigenvectors of R.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the results of the kriging RSS test.
    """
    # Check and process input arguments z, R
    input_locus = locus.copy()
    input_locus = intersect(input_locus)
    z = (input_locus.sumstats[ColName.BETA] / input_locus.sumstats[ColName.SE]).to_numpy()
    n = input_locus.sample_size
    z = np.where(np.isnan(z), 0, z)

    # Compute eigenvalues and eigenvectors
    if eigvens is not None:
        eigvals = eigvens["eigvals"]
        eigvecs = eigvens["eigvecs"]
    else:
        eigens = get_eigen(input_locus.ld.r)
        eigvals = eigens["eigvals"]
        eigvecs = eigens["eigvecs"]
    if s is None:
        s = estimate_s_rss(locus, eigvens={"eigvals": eigvals, "eigvecs": eigvecs})
    eigvals = eigvals[::-1]
    eigvecs = eigvecs[:, ::-1]

    eigvals[eigvals < r_tol] = 0

    if n <= 1:
        raise ValueError("n must be greater than 1")

    sigma2 = (n - 1) / (z**2 + n - 2)
    z = np.sqrt(sigma2) * z

    dinv = 1 / ((1 - s) * eigvals + s)
    dinv[np.isinf(dinv)] = 0
    precision = eigvecs @ (eigvecs * dinv).T
    condmean = np.zeros_like(z)
    condvar = np.zeros_like(z)
    for i in range(len(z)):
        condmean[i] = -(1 / precision[i, i]) * precision[i, :i].dot(z[:i]) - (1 / precision[i, i]) * precision[
            i, i + 1 :
        ].dot(z[i + 1 :])
        condvar[i] = 1 / precision[i, i]
    z_std_diff = (z - condmean) / np.sqrt(condvar)

    # Obtain grid
    a_min = 0.8
    a_max = 2 if np.max(z_std_diff**2) < 1 else 2 * np.sqrt(np.max(z_std_diff**2))
    npoint = int(np.ceil(np.log2(a_max / a_min) / np.log2(1.05)))
    a_grid = 1.05 ** np.arange(-npoint, 1) * a_max

    # Compute likelihood
    sd_mtx = np.outer(np.sqrt(condvar), a_grid)
    matrix_llik = stats.norm.logpdf(z[:, np.newaxis] - condmean[:, np.newaxis], scale=sd_mtx)
    lfactors = np.max(matrix_llik, axis=1)
    matrix_llik = matrix_llik - lfactors[:, np.newaxis]

    # Estimate weight using Gaussian Mixture Model
    gmm = GaussianMixture(n_components=len(a_grid), covariance_type="diag", max_iter=1000)
    gmm.fit(matrix_llik)
    w = gmm.weights_

    # Compute denominators in likelihood ratios
    logl0mix = np.log(np.sum(np.exp(matrix_llik) * (w + 1e-15), axis=1)) + lfactors  # type: ignore

    # Compute numerators in likelihood ratios
    matrix_llik = stats.norm.logpdf(z[:, np.newaxis] + condmean[:, np.newaxis], scale=sd_mtx)
    lfactors = np.max(matrix_llik, axis=1)
    matrix_llik = matrix_llik - lfactors[:, np.newaxis]
    logl1mix = np.log(np.sum(np.exp(matrix_llik) * (w + 1e-15), axis=1)) + lfactors  # type: ignore

    # Compute (log) likelihood ratios
    logLRmix = logl1mix - logl0mix

    res = pd.DataFrame(
        {
            "z": z,
            "condmean": condmean,
            "condvar": condvar,
            "z_std_diff": z_std_diff,
            "logLR": logLRmix,
        },
        index=input_locus.sumstats[ColName.SNPID].to_numpy(),
    )

    # plt.figure(figsize=(5, 5))
    # plt.scatter(condmean, z)
    # plt.xlabel("Expected value")
    # plt.ylabel("Observed z scores")
    # plt.plot([min(condmean), max(condmean)], [min(condmean), max(condmean)], "r--")

    # idx = (logLRmix > 2) & (np.abs(z) > 2)
    # if np.any(idx):
    #     plt.scatter(condmean[idx], z[idx], color="red")

    # plt.title("Observed vs Expected z-scores")
    # plt.tight_layout()

    # return {"plot": plt.gcf(), "conditional_dist": res}

    return res


def compute_dentist_s(locus: Locus) -> pd.DataFrame:
    """
    Compute Dentist-S statistic and p-value.

    Reference: https://github.com/mkanai/slalom/blob/854976f8e19e6fad2db3123eb9249e07ba0e1c1b/slalom.py#L254

    Parameters
    ----------
    locus : Locus
        Locus object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the results of the Dentist-S test.
    """
    input_locus = locus.copy()
    input_locus = intersect(input_locus)
    df = input_locus.sumstats.copy()
    df["Z"] = df[ColName.BETA] / df[ColName.SE]
    lead_idx = df[ColName.P].idxmin()
    # TODO: use abf to select lead variant, although in most cases the lead variant is the one with the smallest p-value
    lead_z = df.loc[lead_idx, ColName.Z]
    df["r"] = input_locus.ld.r[lead_idx]
    df['r'] = df['r'].astype(np.float32)

    df["t_dentist_s"] = (df.Z - df.r * lead_z) ** 2 / (1 - df.r**2)  # type: ignore
    df["t_dentist_s"] = np.where(df["t_dentist_s"] < 0, np.inf, df["t_dentist_s"])
    df.at[lead_idx, "t_dentist_s"] = np.nan
    df["p_dentist_s"] = stats.chi2.logsf(df["t_dentist_s"], df=1)

    df = df[[ColName.SNPID, "t_dentist_s", "p_dentist_s"]]
    df.set_index(ColName.SNPID, inplace=True)
    df.index.name = None
    return df


# TODO: implement HEELS to evaluate the local heritability

def locus_qc(
    locus: Locus,
    r_tol: float = 1e-3,
    method: str = "null-mle",
) -> Dict[str, Any]:
    """
    Quality control for a locus.

    Parameters
    ----------
    locus : Locus
        Locus object.
    r_tol : float, default=1e-3
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.
    method : str, default="null-mle"
        Method to estimate s. Options are "null-mle", "null-partialmle", or "null-pseudomle".

    Returns
    -------
    dict
        Dictionary of quality control results.
    """
    pass
