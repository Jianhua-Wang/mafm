"""Quality control functions for MAFM data."""

import logging
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn.mixture import GaussianMixture

logger = logging.getLogger("QC")


def estimate_s_rss(z: np.ndarray, R: np.ndarray, n: int, r_tol: float = 1e-8, method: str = "null-mle") -> float:
    """
    Estimate s in the susie_rss Model Using Regularized LD.

    This function estimates the parameter s, which provides information about the consistency between z-scores
    and the LD matrix. A larger s indicates a strong inconsistency between z-scores and the LD matrix.

    Parameters
    ----------
    z : np.ndarray
        A p-vector of z-scores.
    R : np.ndarray
        A p by p symmetric, positive semidefinite correlation matrix.
    n : int
        The sample size.
    r_tol : float, default=1e-8
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.
    method : str, default="null-mle"
        Method to estimate s. Options are "null-mle", "null-partialmle", or "null-pseudomle".

    Returns
    -------
    float
        Estimated s value between 0 and 1 (or potentially > 1 for "null-partialmle").

    Examples
    --------
    >>> np.random.seed(1)
    >>> n, p = 500, 1000
    >>> beta = np.zeros(p)
    >>> beta[:4] = 0.01
    >>> X = np.random.randn(n, p)
    >>> X = (X - X.mean(axis=0)) / X.std(axis=0)
    >>> y = X @ beta + np.random.randn(n)
    >>> ss = univariate_regression(X, y)
    >>> R = np.corrcoef(X.T)
    >>> zhat = ss['betahat'] / ss['sebetahat']
    >>> s1 = estimate_s_rss(zhat, R, n)
    """
    # Check and process input arguments z, R
    z = np.where(np.isnan(z), 0, z)
    R = R.astype(np.float32)
    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(R)

    if np.any(eigvals < -r_tol):
        logger.warning("The matrix R is not positive semidefinite. Negative eigenvalues are set to zero")
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


def kriging_rss(z: np.ndarray, R: np.ndarray, n: int, s: float, r_tol: float = 1e-8) -> Dict[str, Any]:
    """
    Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue.

    Under the null, the rss model with regularized LD matrix is z|R,s ~ N(0, (1-s)R + s I)).
    We use a mixture of normals to model the conditional distribution of z_j given other z scores.

    Parameters
    ----------
    z : np.ndarray
        A p-vector of z scores.
    R : np.ndarray
        A p by p symmetric, positive semidefinite correlation matrix.
    n : int
        The sample size.
    s : float
        An estimated s from estimate_s_rss function.
    r_tol : float, default=1e-8
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.

    Returns
    -------
    dict
        A dictionary containing a matplotlib plot object and a pandas DataFrame.
        The plot compares observed z score vs the expected value.
        The DataFrame summarizes the conditional distribution for each variant and the likelihood ratio test.
    """
    # Check and process input arguments z, R
    z = np.where(np.isnan(z), 0, z)

    # Compute eigenvalues and eigenvectors
    R = R.astype(np.float32)
    eigvals, eigvecs = np.linalg.eigh(R)
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
        }
    )

    plt.figure(figsize=(5, 5))
    plt.scatter(condmean, z)
    plt.xlabel("Expected value")
    plt.ylabel("Observed z scores")
    plt.plot([min(condmean), max(condmean)], [min(condmean), max(condmean)], "r--")

    idx = (logLRmix > 2) & (np.abs(z) > 2)
    if np.any(idx):
        plt.scatter(condmean[idx], z[idx], color="red")

    plt.title("Observed vs Expected z-scores")
    plt.tight_layout()

    return {"plot": plt.gcf(), "conditional_dist": res}
