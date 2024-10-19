"""Functions for reading and converting lower triangle matrices."""

import logging
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from sklearn.mixture import GaussianMixture

from mafm.constants import ColName
from mafm.sumstats import make_SNPID_unique, munge_bp, munge_chr

logger = logging.getLogger("LDMatrix")


class LDMatrix:
    """
    Class to store the LD matrix and the corresponding Variant IDs.

    Attributes
    ----------
    map : pd.DataFrame
        DataFrame containing the Variant IDs.
    r : np.ndarray
        LD matrix.
    """

    def __init__(self, map_df: pd.DataFrame, r: np.ndarray):
        """
        Initialize the LDMatrix object.

        Parameters
        ----------
        map_df : pd.DataFrame
            DataFrame containing the Variant IDs.
        r : np.ndarray
            LD matrix.
        """
        self.map = map_df
        self.r = r


def read_lower_triangle(file_path: str, delimiter: str = "\t") -> np.ndarray:
    """
    Read a lower triangle matrix from a file.

    Parameters
    ----------
    file_path : str
        Path to the input text file containing the lower triangle matrix.
    delimiter : str, optional
        Delimiter used in the input file (default is tab).

    Returns
    -------
    np.ndarray
        Lower triangle matrix.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain a valid lower triangle matrix.
    FileNotFoundError
        If the specified file does not exist.
    """
    try:
        with open(file_path, "r") as file:
            rows = [
                list(map(float, line.strip().split(delimiter)))
                for line in file
                if line.strip()
            ]
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    if not rows:
        raise ValueError("The input file is empty.")

    n = len(rows)
    lower_triangle = np.zeros((n, n))

    for i, row in enumerate(rows):
        if len(row) != i + 1:
            raise ValueError(
                f"Invalid number of elements in row {i + 1}. Expected {i + 1}, got {len(row)}."
            )
        lower_triangle[i, : len(row)] = row

    return lower_triangle


def load_ld_matrix(file_path: str, delimiter: str = "\t") -> np.ndarray:
    """
    Convert a lower triangle matrix from a file to a symmetric square matrix.

    Parameters
    ----------
    file_path : str
        Path to the input text file containing the lower triangle matrix.
    delimiter : str, optional
        Delimiter used in the input file (default is tab).

    Returns
    -------
    np.ndarray
        Symmetric square matrix with diagonal filled with 1.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain a valid lower triangle matrix.
    FileNotFoundError
        If the specified file does not exist.

    Notes
    -----
    This function assumes that the input file contains a valid lower triangle matrix
    with each row on a new line and elements separated by the specified delimiter.

    Examples
    --------
    >>> lower_triangle_to_symmetric('lower_triangle.txt')
    array([[1.  , 0.1 , 0.2 , 0.3 ],
            [0.1 , 1.  , 0.4 , 0.5 ],
            [0.2 , 0.4 , 1.  , 0.6 ],
            [0.3 , 0.5 , 0.6 , 1.  ]])
    """
    lower_triangle = read_lower_triangle(file_path, delimiter)

    # Create the symmetric matrix
    symmetric_matrix = lower_triangle + lower_triangle.T

    # Fill the diagonal with 1
    np.fill_diagonal(symmetric_matrix, 1)

    # convert to float16
    symmetric_matrix = symmetric_matrix.astype(np.float16)
    return symmetric_matrix


def load_ld_map(map_path: str, delimiter: str = "\t") -> pd.DataFrame:
    r"""
    Read Variant IDs from a file.

    Parameters
    ----------
    map_path : str
        Path to the input text file containing the Variant IDs.
    delimiter : str, optional
        Delimiter used in the input file (default is tab).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the Variant IDs.

    Raises
    ------
    ValueError
        If the input file is empty or does not contain the required columns.

    Notes
    -----
    This function assumes that the input file contains the required columns:
    - Chromosome (CHR)
    - Base pair position (BP)
    - Allele 1 (A1)
    - Allele 2 (A2)

    Examples
    --------
    >>> contents = "CHR\tBP\tA1\tA2\n1\t1000\tA\tG\n1\t2000\tC\tT\n2\t3000\tT\tC"
    >>> open('map.txt', 'w') as file:
    >>>     file.write(contents)
    >>> load_ld_map('map.txt')
        SNPID   CHR        BP A1 A2
    0   1-1000-A-G  1  1000  A  G
    1   1-2000-C-T  1  2000  C  T
    2   2-3000-C-T  2  3000  T  C
    """
    map_df = pd.read_csv(map_path, sep=delimiter)
    missing_cols = [col for col in ColName.map_cols if col not in map_df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in the input file: {missing_cols}")
    outdf = munge_chr(map_df)
    outdf = munge_bp(outdf)
    for col in [ColName.A1, ColName.A2]:
        pre_n = outdf.shape[0]
        outdf = outdf[outdf[col].notnull()]
        outdf[col] = outdf[col].astype(str).str.upper()
        outdf = outdf[outdf[col].str.match(r"^[ACGT]+$")]
        after_n = outdf.shape[0]
        logger.debug(f"Remove {pre_n - after_n} rows because of invalid {col}.")
    outdf = outdf[outdf[ColName.A1] != outdf[ColName.A2]]
    outdf = make_SNPID_unique(
        outdf, col_ea=ColName.A1, col_nea=ColName.A2, remove_duplicates=False
    )
    outdf.reset_index(drop=True, inplace=True)
    return outdf


def sort_alleles(ld: LDMatrix) -> LDMatrix:
    """
    Sort alleles in the LD map in alphabetical order. Change the sign of the LD matrix if the alleles are swapped.

    Parameters
    ----------
    ld : LDMatrix
        Dictionary containing the Variant IDs and the LD matrix.

    Returns
    -------
    LDMatrix
        Dictionary containing the Variant IDs and the LD matrix with alleles sorted.

    Examples
    --------
    >>> ld = {
    ...     'map': pd.DataFrame({
    ...         'SNPID': ['1-1000-A-G', '1-2000-C-T'],
    ...         'CHR': [1, 1],
    ...         'BP': [1000, 2000],
    ...         'A1': ['A', 'T'],
    ...         'A2': ['G', 'C']
    ...     }),
    ...     'r': np.array([[1. , 0.1],
    ...                    [0.1, 1. ]])
    ... }
    >>> ld = LDMatrix(**ld)
    >>> sort_alleles(ld)
    LDMatrix(map=   SNPID  CHR    BP A1 A2
    0  1-1000-A-G    1  1000  A  G
    1  1-2000-C-T    1  2000  C  T, r=array([[ 1. , -0.1],
            [-0.1,  1. ]]))
    """
    ld_df = ld.r.copy()
    ld_map = ld.map.copy()
    ld_map[["sort_a1", "sort_a2"]] = np.sort(ld_map[[ColName.A1, ColName.A2]], axis=1)
    swapped_index = ld_map[ld_map[ColName.A1] != ld_map["sort_a1"]].index
    # Change the sign of the rows and columns the LD matrix if the alleles are swapped
    ld_df[swapped_index] *= -1
    ld_df[:, swapped_index] *= -1
    np.fill_diagonal(ld_df, 1)

    ld_map[ColName.A1] = ld_map["sort_a1"]
    ld_map[ColName.A2] = ld_map["sort_a2"]
    ld_map.drop(columns=["sort_a1", "sort_a2"], inplace=True)
    return LDMatrix(ld_map, ld_df)


def estimate_s_rss(
    z: np.ndarray, R: np.ndarray, n: int, r_tol: float = 1e-8, method: str = "null-mle"
) -> float:
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

    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(R)

    if np.any(eigvals < -r_tol):
        logger.warning(
            "The matrix R is not positive semidefinite. Negative eigenvalues are set to zero"
        )
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

        def pseudolikelihood(
            s: float, z: np.ndarray, eigvals: np.ndarray, eigvecs: np.ndarray
        ) -> float:
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
    z: np.ndarray, R: np.ndarray, n: int, s: float, r_tol: float = 1e-8
) -> Dict[str, Any]:
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
        condmean[i] = -(1 / precision[i, i]) * precision[i, :i].dot(z[:i]) - (
            1 / precision[i, i]
        ) * precision[i, i + 1 :].dot(z[i + 1 :])
        condvar[i] = 1 / precision[i, i]
    z_std_diff = (z - condmean) / np.sqrt(condvar)

    # Obtain grid
    a_min = 0.8
    a_max = 2 if np.max(z_std_diff**2) < 1 else 2 * np.sqrt(np.max(z_std_diff**2))
    npoint = int(np.ceil(np.log2(a_max / a_min) / np.log2(1.05)))
    a_grid = 1.05 ** np.arange(-npoint, 1) * a_max

    # Compute likelihood
    sd_mtx = np.outer(np.sqrt(condvar), a_grid)
    matrix_llik = stats.norm.logpdf(
        z[:, np.newaxis] - condmean[:, np.newaxis], scale=sd_mtx
    )
    lfactors = np.max(matrix_llik, axis=1)
    matrix_llik = matrix_llik - lfactors[:, np.newaxis]

    # Estimate weight using Gaussian Mixture Model
    gmm = GaussianMixture(
        n_components=len(a_grid), covariance_type="diag", max_iter=1000
    )
    gmm.fit(matrix_llik)
    w = gmm.weights_

    # Compute denominators in likelihood ratios
    logl0mix = np.log(np.sum(np.exp(matrix_llik) * (w + 1e-15), axis=1)) + lfactors  # type: ignore

    # Compute numerators in likelihood ratios
    matrix_llik = stats.norm.logpdf(
        z[:, np.newaxis] + condmean[:, np.newaxis], scale=sd_mtx
    )
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


def load_ld(
    ld_path: str, map_path: str, delimiter: str = "\t", if_sort_alleles: bool = True
) -> LDMatrix:
    """
    Read LD matrices and Variant IDs from files. Pair each matrix with its corresponding Variant IDs.

    Parameters
    ----------
    ld_path : str
        Path to the input text file containing the lower triangle matrix.
    map_path : str
        Path to the input text file containing the Variant IDs.
    delimiter : str, optional
        Delimiter used in the input file (default is tab).
    if_sort_alleles : bool, optional
        Sort alleles in the LD map in alphabetical order and change the sign of the LD matrix if the alleles are swapped
        (default is True).

    Returns
    -------
    LDMatrix
        Object containing the LD matrix and the Variant IDs.

    Raises
    ------
    ValueError
        If the number of variants in the map file does not match the number of rows in the LD matrix.
    """
    ld_df = load_ld_matrix(ld_path, delimiter)
    logger.info(f"Loaded LD matrix with shape {ld_df.shape} from '{ld_path}'.")
    map_df = load_ld_map(map_path, delimiter)
    logger.info(f"Loaded map file with shape {map_df.shape} from '{map_path}'.")
    if ld_df.shape[0] != map_df.shape[0]:
        raise ValueError(
            "The number of variants in the map file does not match the number of rows in the LD matrix."
            f"Number of variants in the map file: {map_df.shape[0]}, number of rows in the LD matrix: {ld_df.shape[0]}"
        )
    ld = LDMatrix(map_df, ld_df)
    if if_sort_alleles:
        ld = sort_alleles(ld)

    return ld
