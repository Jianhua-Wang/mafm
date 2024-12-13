"""Quality control functions for MAFM data."""

import logging
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit, minimize_scalar
from sklearn.mixture import GaussianMixture

from mafm.constants import ColName
from mafm.locus import Locus, LocusSet, intersect_sumstat_ld

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
    # ldmatrix = ldmatrix.astype(np.float32)
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
    input_locus = intersect_sumstat_ld(input_locus)
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

    # if np.any(eigvals < -r_tol):
    #     logger.warning("The LD matrix is not positive semidefinite. Negative eigenvalues are set to zero")
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
    input_locus = intersect_sumstat_ld(input_locus)
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
            "SNPID": input_locus.sumstats[ColName.SNPID].to_numpy(),
            "z": z,
            "condmean": condmean,
            "condvar": condvar,
            "z_std_diff": z_std_diff,
            "logLR": logLRmix,
        },
        # index=input_locus.sumstats[ColName.SNPID].to_numpy(),
    )
    # TODO: remove variants with logLR > 2 and abs(z) > 2

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
    input_locus = intersect_sumstat_ld(input_locus)
    df = input_locus.sumstats.copy()
    df["Z"] = df[ColName.BETA] / df[ColName.SE]
    lead_idx = df[ColName.P].idxmin()
    # TODO: use abf to select lead variant, although in most cases the lead variant is the one with the smallest p-value
    lead_z = df.loc[lead_idx, ColName.Z]
    df["r"] = input_locus.ld.r[lead_idx]

    df["t_dentist_s"] = (df.Z - df.r * lead_z) ** 2 / (1 - df.r**2)  # type: ignore
    df["t_dentist_s"] = np.where(df["t_dentist_s"] < 0, np.inf, df["t_dentist_s"])
    df.at[lead_idx, "t_dentist_s"] = np.nan
    df["p_dentist_s"] = stats.chi2.logsf(df["t_dentist_s"], df=1)

    df = df[[ColName.SNPID, "t_dentist_s", "p_dentist_s"]].copy()
    # df.set_index(ColName.SNPID, inplace=True)
    # df.index.name = None
    return df


def compare_maf(locus: Locus) -> pd.DataFrame:
    """
    Compare the allele frequency in the sumstats and the allele frequency in the LD reference.

    Parameters
    ----------
    locus : Locus
        Locus object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the results of the comparison.
    """
    input_locus = locus.copy()
    input_locus = intersect_sumstat_ld(input_locus)
    df = input_locus.sumstats[[ColName.SNPID, ColName.MAF]].copy()
    df.rename(columns={ColName.MAF: "MAF_sumstats"}, inplace=True)
    df.set_index(ColName.SNPID, inplace=True)
    af_ld = pd.Series(index=input_locus.ld.map[ColName.SNPID], data=input_locus.ld.map["AF2"])
    maf_ld = np.minimum(af_ld, 1 - af_ld)
    df["MAF_ld"] = maf_ld
    df[ColName.SNPID] = df.index
    df.reset_index(drop=True, inplace=True)
    return df


def snp_missingness(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the missingness of each cohort.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the missingness of each cohort.
    """
    missingness_df = []
    for locus in locus_set.loci:
        loc = intersect_sumstat_ld(locus)
        loc = loc.sumstats[[ColName.SNPID]].copy()
        loc[f"{locus.popu}_{locus.cohort}"] = 1
        loc.set_index(ColName.SNPID, inplace=True)
        missingness_df.append(loc)
    missingness_df = pd.concat(missingness_df, axis=1)
    missingness_df.fillna(0, inplace=True)
    # log warning if missing rate > 0.1
    for col in missingness_df.columns:
        missing_rate = float(round(1 - missingness_df[col].sum() / missingness_df.shape[0], 3))
        if missing_rate > 0.1:
            logger.warning(f"The missing rate of {col} is {missing_rate}")
        else:
            logger.info(f"The missing rate of {col} is {missing_rate}")

    return missingness_df


def ld_4th_moment(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the 4th moment of the LD matrix.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the 4th moment of the LD matrix.
    """
    ld_4th_res = []
    # intersect between loci
    overlap_snps = set(locus_set.loci[0].sumstats[ColName.SNPID])
    for locus in locus_set.loci[1:]:
        overlap_snps = overlap_snps.intersection(set(locus.sumstats[ColName.SNPID]))
    for locus in locus_set.loci:
        locus = locus.copy()
        locus.sumstats = locus.sumstats[locus.sumstats[ColName.SNPID].isin(overlap_snps)]
        locus = intersect_sumstat_ld(locus)
        r_4th = pd.Series(index=locus.ld.map[ColName.SNPID], data=np.power(locus.ld.r, 4).sum(axis=0))
        r_4th = r_4th - 1
        r_4th.name = f"{locus.popu}_{locus.cohort}"
        ld_4th_res.append(r_4th)
    return pd.concat(ld_4th_res, axis=1)


def ld_decay(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the decay of the LD matrix.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the decay of the LD matrix.
    """

    def fit_exp(x, a, b):
        with np.errstate(over="ignore"):
            return a * np.exp(-b * x)

    binsize = 1000
    decay_res = []
    for locus in locus_set.loci:
        ldmap = locus.ld.map.copy()
        r = locus.ld.r.copy()
        distance_mat = np.array([ldmap["BP"] - ldmap["BP"].values[i] for i in range(len(ldmap))])
        distance_mat = distance_mat[np.tril_indices_from(distance_mat, k=-1)].flatten()
        distance_mat = np.abs(distance_mat)
        r = r[np.tril_indices_from(r, k=-1)].flatten()
        r = np.square(r)
        bins = np.arange(0, ldmap["BP"].max() - ldmap["BP"].min() + binsize, binsize)

        r_sum, _ = np.histogram(distance_mat, bins=bins, weights=r)
        count, _ = np.histogram(distance_mat, bins=bins)

        with np.errstate(divide="ignore", invalid="ignore"):
            r2_avg = np.where(count > 0, r_sum / count, 0)
        popt, _ = curve_fit(fit_exp, bins[1:] / binsize, r2_avg)
        res = pd.DataFrame(
            {
                "distance_kb": bins[1:] / binsize,
                "r2_avg": r2_avg,
                "decay_rate": popt[0],
                "cohort": f"{locus.popu}_{locus.cohort}",
            }
        )
        decay_res.append(res)
    return pd.concat(decay_res, axis=0)


def cochran_q(locus_set: LocusSet) -> pd.DataFrame:
    """
    Compute the Cochran-Q statistic.

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the Cochran-Q statistic.
    """
    merged_df = locus_set.loci[0].original_sumstats[[ColName.SNPID]].copy()
    for i, df in enumerate(locus_set.loci):
        df = df.sumstats[[ColName.SNPID, ColName.BETA, ColName.SE, ColName.EAF]].copy()
        df.rename(columns={ColName.BETA: f"BETA_{i}", ColName.SE: f"SE_{i}", ColName.EAF: f"EAF_{i}"}, inplace=True)
        merged_df = pd.merge(merged_df, df, on=ColName.SNPID, how="inner", suffixes=("", f"_{i}"))

    k = len(locus_set.loci)
    weights = []
    effects = []
    for i in range(k):
        weights.append((1 / (merged_df[f"SE_{i}"] ** 2)))
        effects.append(merged_df[f"BETA_{i}"])

    # Calculate weighted mean effect size
    weighted_mean = np.sum([w * e for w, e in zip(weights, effects)], axis=0) / np.sum(weights, axis=0)

    # Calculate Q statistic
    Q = np.sum([w * (e - weighted_mean) ** 2 for w, e in zip(weights, effects)], axis=0)

    # Calculate degrees of freedom
    df = k - 1

    # Calculate P-value
    p_value = stats.chi2.sf(Q, df)

    # Calculate I^2
    with np.errstate(invalid="ignore"):
        I_squared = np.maximum(0, (Q - df) / Q * 100)

    # Create output dataframe
    output_df = pd.DataFrame({"SNPID": merged_df["SNPID"], "Q": Q, "Q_pvalue": p_value, "I_squared": I_squared})
    return output_df.set_index(ColName.SNPID)


def locus_qc(
    locus_set: LocusSet,
    r_tol: float = 1e-3,
    method: str = "null-mle",
):
    """
    Quality control for a locus.

    TODO: add LAVA

    Parameters
    ----------
    locus_set : LocusSet
        LocusSet object.
    r_tol : float, default=1e-3
        Tolerance level for eigenvalue check of positive semidefinite matrix of R.
    method : str, default="null-mle"
        Method to estimate s. Options are "null-mle", "null-partialmle", or "null-pseudomle".

    Returns
    -------
    dict
        Dictionary of quality control results.
    """
    qc_metrics = {}
    all_expected_z = []
    all_dentist_s = []
    all_compare_maf = []
    for locus in locus_set.loci:
        lo = intersect_sumstat_ld(locus)
        eigens = get_eigen(lo.ld.r)
        lambda_s = estimate_s_rss(locus, r_tol, method, eigens)
        expected_z = kriging_rss(locus, r_tol, lambda_s, eigens)
        expected_z["lambda_s"] = lambda_s
        expected_z["cohort"] = f"{locus.popu}_{locus.cohort}"
        dentist_s = compute_dentist_s(locus)
        dentist_s["cohort"] = f"{locus.popu}_{locus.cohort}"
        compare_maf_res = compare_maf(locus)
        compare_maf_res["cohort"] = f"{locus.popu}_{locus.cohort}"
        all_expected_z.append(expected_z)
        all_dentist_s.append(dentist_s)
        all_compare_maf.append(compare_maf_res)
    all_expected_z = pd.concat(all_expected_z, axis=0)
    all_dentist_s = pd.concat(all_dentist_s, axis=0)
    all_compare_maf = pd.concat(all_compare_maf, axis=0)
    qc_metrics["expected_z"] = all_expected_z
    qc_metrics["dentist_s"] = all_dentist_s
    qc_metrics["compare_maf"] = all_compare_maf

    qc_metrics["ld_4th_moment"] = ld_4th_moment(locus_set)
    qc_metrics["ld_decay"] = ld_decay(locus_set)

    if len(locus_set.loci) > 1:
        qc_metrics["cochran_q"] = cochran_q(locus_set)
        qc_metrics["snp_missingness"] = snp_missingness(locus_set)

    return qc_metrics
