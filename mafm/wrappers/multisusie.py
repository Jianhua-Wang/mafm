"""Wrapper for MultiSuSiE."""

import json
import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
from mafm.locus import LocusSet, intersect_sumstat_ld
from mafm.wrappers.multisusie_rss import multisusie_rss

logger = logging.getLogger("MULTISUSIE")


def run_multisusie(
    locus_set: LocusSet,
    max_causal: int = 1,
    coverage: float = 0.95,
    rho: float = 0.75,
    scaled_prior_variance: float = 0.2,
    standardize: bool = False,
    pop_spec_standardization: bool = True,
    estimate_residual_variance: bool = True,
    estimate_prior_variance: bool = True,
    estimate_prior_method: str = "early_EM",
    pop_spec_effect_priors: bool = True,
    iter_before_zeroing_effects: int = 5,
    prior_tol: float = 1e-9,
    max_iter: int = 100,
    tol: float = 1e-3,
    min_abs_corr: float = 0,
):
    """
    Run MultiSuSiE.

    Parameters
    ----------
    locus_set : LocusSet
        The LocusSet to run MultiSuSiE on.
    max_causal : int, optional
        Maximum number of causal variants, by default 1.
    coverage : float, optional
        Coverage of the credible set, by default 0.95.
    rho : float, optional
        The prior correlation between causal variants, by default 0.75.
    scaled_prior_variance : float, optional
        The scaled prior variance, by default 0.2.
    standardize : bool, optional
        Whether to standardize the data, by default False.
    pop_spec_standardization : bool, optional
        Whether to perform population-specific standardization, by default True.
    estimate_residual_variance : bool, optional
        Whether to estimate the residual variance, by default True.
    estimate_prior_variance : bool, optional
        Whether to estimate the prior variance, by default True.
    estimate_prior_method : str, optional
        The method to estimate the prior variance, by default "early_EM".
    pop_spec_effect_priors : bool, optional
        Whether to use population-specific effect priors, by default True.
    iter_before_zeroing_effects : int, optional
        The number of iterations before zeroing out effects, by default 5.
    prior_tol : float, optional
        The tolerance for the prior, by default 1e-9.
    max_iter : int, optional
        The maximum number of iterations, by default 100.
    tol : float, optional
        The tolerance for convergence, by default 1e-3.
    min_abs_corr : float, optional
        The minimum absolute correlation, by default 0.

    Returns
    -------
    CredibleSet
        The credible set.
    """
    logger.info(f"Running MultiSuSiE on {locus_set}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "rho": rho,
        "scaled_prior_variance": scaled_prior_variance,
        "standardize": standardize,
        "pop_spec_standardization": pop_spec_standardization,
        "estimate_residual_variance": estimate_residual_variance,
        "estimate_prior_variance": estimate_prior_variance,
        "estimate_prior_method": estimate_prior_method,
        "pop_spec_effect_priors": pop_spec_effect_priors,
        "iter_before_zeroing_effects": iter_before_zeroing_effects,
        "prior_tol": prior_tol,
        "max_iter": max_iter,
        "tol": tol,
        "min_abs_corr": min_abs_corr,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")

    all_variants = []
    for locus in locus_set.loci:
        locus = intersect_sumstat_ld(locus)
        all_variants.append(locus.sumstats[[ColName.SNPID, ColName.CHR, ColName.BP, ColName.EA, ColName.NEA]])
    all_variants = pd.concat(all_variants, axis=0)
    all_variants.drop_duplicates(subset=[ColName.SNPID], inplace=True)
    all_variants.sort_values([ColName.CHR, ColName.BP], inplace=True)
    variant_to_index = {variant: i for i, variant in enumerate(all_variants["SNPID"])}
    # TODO: make concat of all_variants as a function of locus_set
    # TODO: make intersect of loci as a function of all_variants
    # TODO: add a switch of either using all_variants or the joint_variants
    # joint_variants = []
    # for i, locus in enumerate(locus_set.loci):
    #     locus = intersect_sumstat_ld(locus)
    #     sumstat = locus.sumstats.copy()
    #     if i == 0:
    #         joint_variants = set(sumstat["SNPID"].values)
    #     else:
    #         joint_variants = joint_variants.intersection(set(sumstat["SNPID"].values))
    # joint_variants = sumstat[sumstat["SNPID"].isin(joint_variants)][["SNPID", "CHR", "BP", "EA", "NEA"]].copy()
    # joint_variants.sort_values(["CHR", "BP"], inplace=True)
    # all_variants = joint_variants.copy()
    # variant_to_index = {variant: i for i, variant in enumerate(all_variants["SNPID"])}

    z_list = []
    R_list = []
    for locus in locus_set.loci:
        sumstat = locus.sumstats.copy()
        ldmap = locus.ld.map
        ld = locus.ld.r
        sumstat[ColName.Z] = sumstat[ColName.BETA] / sumstat[ColName.SE]
        sumstat.set_index(ColName.SNPID, inplace=True)
        z = all_variants[ColName.SNPID].map(sumstat[ColName.Z]).values
        z_list.append(z)
        expand_ld = np.zeros((all_variants.shape[0], all_variants.shape[0]))
        intersec_index = ldmap[ldmap["SNPID"].isin(all_variants["SNPID"])].index
        ldmap = ldmap.loc[intersec_index]
        ld = ld[intersec_index, :][:, intersec_index]
        study_indices = np.array([variant_to_index[snp] for snp in ldmap["SNPID"]])
        idx_i, idx_j = np.meshgrid(study_indices, study_indices)
        expand_ld[idx_i, idx_j] += ld.astype(np.float32)
        np.fill_diagonal(expand_ld, 1)
        R_list.append(expand_ld)

    rho_array = np.full((len(locus_set.loci), len(locus_set.loci)), rho)
    np.fill_diagonal(rho_array, 1)
    ss_fit = multisusie_rss(
        z_list=z_list,
        R_list=R_list,
        population_sizes=[locus.sample_size for locus in locus_set.loci],
        rho=rho_array,  # type: ignore
        L=max_causal,
        coverage=coverage,
        scaled_prior_variance=scaled_prior_variance,
        max_iter=max_iter,
        tol=tol,
        pop_spec_standardization=pop_spec_standardization,
        estimate_residual_variance=estimate_residual_variance,
        estimate_prior_variance=estimate_prior_variance,
        estimate_prior_method=estimate_prior_method,
        pop_spec_effect_priors=pop_spec_effect_priors,
        iter_before_zeroing_effects=iter_before_zeroing_effects,
        prior_tol=prior_tol,
        min_abs_corr=min_abs_corr,
        float_type=np.float32,
        low_memory_mode=False,
        recover_R=False,
        single_population_mac_thresh=20,
        mac_list=None,
        multi_population_maf_thresh=0,
        maf_list=None,
    )
    pip = pd.Series(index=all_variants[ColName.SNPID].tolist(), data=ss_fit.pip)
    cs_snp = []
    for cs_snp_idx in ss_fit.sets[0]:
        if len(cs_snp_idx) > 0 and len(cs_snp_idx) < len(pip):
            cs_snp.append(all_variants[ColName.SNPID].to_numpy()[cs_snp_idx])
    cs_sizes = [len(snpids) for snpids in cs_snp]
    lead_snps = [str(pip[pip.index.isin(snpids)].idxmax()) for snpids in cs_snp]

    logger.info(f"Finished MultiSuSiE on {locus_set}")
    logger.info(f"N of credible set: {len(cs_snp)}")
    logger.info(f"Credible set size: {[len(i) for i in cs_snp]}")

    return CredibleSet(
        tool=Method.MULTISUSIE,
        n_cs=len(cs_snp),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snp,
        cs_sizes=cs_sizes,
        pips=pip,
        parameters=parameters,
    )
