"""
RSparsePro wrapper for multi-ancestry fine-mapping.

Original code from https://github.com/zhwm/RSparsePro_LD
"""

import json
import logging

import numpy as np
import pandas as pd
from scipy.special import softmax

from mafm.constants import ColName, Method
from mafm.credibleset import CredibleSet
from mafm.locus import Locus, intersect_sumstat_ld

logger = logging.getLogger("RSparsePro")


class RSparsePro(object):
    """
    RSparsePro for robust fine-mapping in the presence of LD mismatch.

    This class implements the RSparsePro algorithm for robust fine-mapping in the presence of LD mismatch.

    Attributes
    ----------
        p : int
            The number of variants.
        k : int
            The number of causal signals.
        vare : float
            The error parameter.
        mat : numpy.ndarray
            The matrix used in the algorithm.
        beta_mu : numpy.ndarray
            The posterior mean of the effect sizes.
        gamma : numpy.ndarray
            The posterior inclusion probabilities.
        tilde_b : numpy.ndarray
            The estimated effect sizes.

    Methods
    -------
        infer_q_beta(R): Infer the posterior mean of the effect sizes.
        infer_tilde_b(bhat): Infer the posterior mean of the effect sizes.
        train(bhat, R, maxite, eps, ubound): Train the RSparsePro model.
        get_PIP(): Get the posterior inclusion probabilities.
        get_effect(cthres): Get the effect sizes.
        get_ztilde(): Get the estimated effect sizes.
        get_eff_maxld(eff, ld): Get the maximum LD within effect groups.
        get_eff_minld(eff, ld): Get the minimum LD within effect groups.
        get_ordered(eff_mu): Check if the effect sizes are ordered.
        adaptive_train(zscore, ld, K, maxite, eps, ubound, cthres, minldthres, maxldthres, eincre, varemax, varemin): Train the RSparsePro model.
        parse_args(): Parse command line arguments.
    """

    def __init__(self, P, K, R, vare):
        """
        Initialize the RSparsePro model.

        Parameters
        ----------
            P : int
                The number of variants.
            K : int
                The number of causal signals.
            R : numpy.ndarray
                The LD matrix.
            vare : float
                The error parameter.
        """
        self.p = P
        self.k = K
        self.vare = vare
        if vare != 0:
            self.mat = np.dot(R, np.linalg.inv(np.eye(self.p) + 1 / vare * R))
        self.beta_mu = np.zeros([self.p, self.k])
        self.gamma = np.zeros([self.p, self.k])
        self.tilde_b = np.zeros((self.p,))

    def infer_q_beta(self, R):
        """
        Infer the posterior mean of the effect sizes.

        Parameters
        ----------
            R : numpy.ndarray
                The LD matrix.
        """
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:, idxall] * self.beta_mu[:, idxall]).sum(axis=1)
            res_beta = self.tilde_b - np.dot(R, beta_all_k)
            self.beta_mu[:, k] = res_beta
            u = 0.5 * self.beta_mu[:, k] ** 2
            self.gamma[:, k] = softmax(u)

    def infer_tilde_b(self, bhat):
        """
        Infer the posterior mean of the effect sizes.

        Parameters
        ----------
            bhat : numpy.ndarray
                The summary statistics.
        """
        if self.vare == 0:
            self.tilde_b = bhat
        else:
            beta_all = (self.gamma * self.beta_mu).sum(axis=1)
            self.tilde_b = np.dot(self.mat, (1 / self.vare * bhat + beta_all))

    def train(self, bhat, R, maxite, eps, ubound):
        """
        Train the RSparsePro model.

        Parameters
        ----------
            bhat : numpy.ndarray
                The summary statistics.
            R : numpy.ndarray
                The LD matrix.
            maxite : int
                The maximum number of iterations.
            eps : float
                The convergence criterion.
            ubound : int
                The upper bound for convergence.
        """
        for ite in range(maxite):
            old_gamma = self.gamma.copy()
            old_beta = self.beta_mu.copy()
            old_tilde = self.tilde_b.copy()
            self.infer_tilde_b(bhat)
            self.infer_q_beta(R)
            diff_gamma = np.linalg.norm(self.gamma - old_gamma)
            diff_beta = np.linalg.norm(self.beta_mu - old_beta)
            diff_b = np.linalg.norm(self.tilde_b - old_tilde)
            all_diff = diff_gamma + diff_beta + diff_b
            logger.info(
                "Iteration-->{} . Diff_b: {:.1f} . Diff_s: {:.1f} . Diff_mu: {:.1f} . ALL: {:.1f}".format(
                    ite, diff_b, diff_gamma, diff_beta, all_diff
                )
            )
            if all_diff < eps:
                logger.info("The RSparsePro algorithm has converged.")
                converged = True
                break
            if ite == (maxite - 1) or abs(all_diff) > ubound:
                logger.info("The RSparsePro algorithm didn't converge.")
                converged = False
                break
        return converged

    def get_PIP(self):
        """Get the posterior inclusion probabilities."""
        return np.max((self.gamma), axis=1).round(4)

    def get_effect(self, cthres):
        """
        Get the effect sizes.

        Parameters
        ----------
            cthres : float
                The threshold for the coverage.
        """
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        mat_eff[mat_eff < 1 / (self.p + 1)] = 0
        csum = mat_eff.sum(axis=0).round(2)
        logger.info("Attainable coverage for effect groups: {}".format(csum))
        eff = {}
        eff_gamma = {}
        eff_mu = {}
        for k in range(self.k):
            if csum[k] >= cthres:
                p = 0
                while np.sum(mat_eff[matidx[0:p, k], k]) < cthres * csum[k]:
                    p = p + 1
                cidx = matidx[0:p, k].tolist()
                eff[k] = cidx
                eff_gamma[k] = mat_eff[cidx, k].round(4)
                eff_mu[k] = self.beta_mu[cidx, k].round(4)
        return eff, eff_gamma, eff_mu

    def get_ztilde(self):
        """Get the estimated effect sizes."""
        return self.tilde_b.round(4)

    # def get_resz(self, bhat, ld, eff):
    #    idx = [i[0] for i in eff.values()]
    #    realmu = np.zeros(len(bhat))
    #    realmu[idx] = np.dot(np.linalg.inv(ld[np.ix_(idx, idx)]), bhat[idx])
    #    estz = np.dot(ld, realmu)
    #    resz = bhat - estz
    #    return resz.round(4)


def get_eff_maxld(eff, ld):
    """
    Get the maximum LD within effect groups.

    Parameters
    ----------
        eff : dict
            The effect sizes.
        ld : numpy.ndarray
            The LD matrix.
    """
    idx = [i[0] for i in eff.values()]
    if len(eff) > 1:
        maxld = np.abs(np.tril(ld[np.ix_(idx, idx)], -1)).max()
    else:
        maxld = 0.0
    return maxld


def get_eff_minld(eff, ld):
    """
    Get the minimum LD within effect groups.

    Parameters
    ----------
        eff : dict
            The effect sizes.
        ld : numpy.ndarray
            The LD matrix.
    """
    if len(eff) == 0:
        minld = 1.0
    else:
        minld = min([abs(ld[np.ix_(v, v)]).min() for _, v in eff.items()])
    return minld


def get_ordered(eff_mu):
    """
    Check if the effect sizes are ordered.

    Parameters
    ----------
        eff_mu : dict
            The effect sizes.
    """
    if len(eff_mu) > 1:
        val_mu = [round(-abs(i[0])) for _, i in eff_mu.items()]
        ordered = list(eff_mu.keys())[-1] == len(eff_mu) - 1  # and (sorted(val_mu) == val_mu)
    else:
        ordered = True
    return ordered


def adaptive_train(zscore, ld, K, maxite, eps, ubound, cthres, minldthres, maxldthres, eincre, varemax, varemin):
    """
    Train the RSparsePro model.

    Parameters
    ----------
        zscore : numpy.ndarray
            The summary statistics.
        ld : numpy.ndarray
            The LD matrix.
        K : int
            The number of causal signals.
    """
    vare = 0
    mc = False
    eff = {}
    eff_mu = {}
    minld = 1.0
    maxld = 0.0
    while (not mc) or (not get_ordered(eff_mu)) or (minld < minldthres) or (maxld > maxldthres):
        model = RSparsePro(len(zscore), K, ld, vare)
        mc = model.train(zscore, ld, maxite, eps, ubound)
        eff, eff_gamma, eff_mu = model.get_effect(cthres)
        maxld = get_eff_maxld(eff, ld)
        minld = get_eff_minld(eff, ld)
        logging.info("Max ld across effect groups: {}.".format(maxld))
        logging.info("Min ld within effect groups: {}.".format(minld))
        logging.info("vare = {}".format(round(vare, 4)))
        if vare > varemax or (len(eff) < 2 and get_ordered(eff_mu)):
            # logging.info("Algorithm didn't converge at the max vare. Setting K to 1.")
            model = RSparsePro(len(zscore), 1, ld, 0)
            mc = model.train(zscore, ld, maxite, eps, ubound)
            eff, eff_gamma, eff_mu = model.get_effect(cthres)
            break
        elif vare == 0:
            vare = varemin
        else:
            vare *= eincre
    ztilde = model.get_ztilde()
    # resz = model.get_resz(zscore, ld, eff)
    PIP = model.get_PIP()
    return eff, eff_gamma, eff_mu, PIP, ztilde  # resz


def rsparsepro_main(
    zfile: pd.DataFrame,
    ld: np.ndarray,
    K: int = 10,
    maxite: int = 100,
    eps: float = 1e-5,
    ubound: int = 100000,
    cthres: float = 0.95,
    eincre: float = 1.5,
    minldthres: float = 0.7,
    maxldthres: float = 0.2,
    varemax: float = 100.0,
    varemin: float = 1e-3,
) -> pd.DataFrame:
    """
    Run RSparsePro.

    Parameters
    ----------
        zfile : pandas.DataFrame
            The summary statistics.
        ld : numpy.ndarray
            The LD matrix.
        K : int
            The number of causal signals.
        maxite : int
            The maximum number of iterations.
        eps : float
            The convergence criterion.
        ubound : int
            The upper bound for convergence.
        cthres : float
            The threshold for the coverage.
        eincre : float
            The adjustment for the error parameter.
        minldthres : float
            The threshold for the minimum LD within effect groups.
        maxldthres : float
            The threshold for the maximum LD across effect groups.
        varemax : float
            The maximum error parameter.
        varemin : float
            The minimum error parameter.

    Returns
    -------
        eff : dict
            The effect sizes.
        eff_gamma : dict
            The posterior inclusion probabilities.
        eff_mu : dict
            The estimated effect sizes.
        PIP : numpy.ndarray
            The posterior inclusion probabilities.
        ztilde : numpy.ndarray
            The estimated effect sizes.
    """
    eff, eff_gamma, eff_mu, PIP, ztilde = adaptive_train(
        zfile["Z"],
        ld,
        K,
        maxite,
        eps,
        ubound,
        cthres,
        minldthres,
        maxldthres,
        eincre,
        varemax,
        varemin,
    )
    zfile["PIP"] = PIP
    zfile["z_estimated"] = ztilde
    zfile["cs"] = 0
    for e in eff:
        mcs_idx = [zfile["RSID"][j] for j in eff[e]]
        logger.info(f"The {e}-th effect group contains effective variants:")
        logger.info(f"causal variants: {mcs_idx}")
        logger.info(f"variant probabilities for this effect group: {eff_gamma[e]}")
        logger.info(f"zscore for this effect group: {eff_mu[e]}\n")
        zfile.loc[list(eff[e]), "cs"] = e + 1
    # zfile.to_csv("{}.rsparsepro.txt".format(save), sep="\t", header=True, index=False)
    return zfile


def run_rsparsepro(
    locus: Locus,
    max_causal: int = 1,
    coverage: float = 0.95,
    maxite: int = 100,
    eps: float = 1e-5,
    ubound: int = 100000,
    eincre: float = 1.5,
    minldthres: float = 0.7,
    maxldthres: float = 0.2,
    varemax: float = 100.0,
    varemin: float = 1e-3,
) -> CredibleSet:
    """
    Run RSparsePro.

    Parameters
    ----------
        locus : Locus
            The locus.
        max_causal : int
            The maximum number of causal signals.
        coverage : float
            The coverage.
        maxite : int
            The maximum number of iterations.
        ubound : int
            The upper bound for convergence.
        cthres : float
            The threshold for the coverage.
        eincre : float
            The adjustment for the error parameter.
        minldthres : float
            The threshold for the minimum LD within effect groups.
        maxldthres : float
            The threshold for the maximum LD across effect groups.
        varemax : float
            The maximum error parameter.
        varemin : float
            The minimum error parameter.

    Returns
    -------
        CredibleSet
            The credible set.
    """
    if not locus.is_matched:
        logger.warning("The sumstat and LD are not matched, will match them in same order.")
        locus = intersect_sumstat_ld(locus)
    logger.info(f"Running RSparsePro on {locus}")
    parameters = {
        "max_causal": max_causal,
        "coverage": coverage,
        "maxite": maxite,
        "eps": eps,
        "ubound": ubound,
        "eincre": eincre,
        "minldthres": minldthres,
        "maxldthres": maxldthres,
        "varemax": varemax,
        "varemin": varemin,
    }
    logger.info(f"Parameters: {json.dumps(parameters, indent=4)}")

    sumstats = locus.sumstats.copy()
    ld = locus.ld.r.copy()
    sumstats["RSID"] = sumstats[ColName.SNPID]
    sumstats[ColName.Z] = sumstats[ColName.BETA] / sumstats[ColName.SE]
    zfile = rsparsepro_main(
        zfile=sumstats,
        ld=ld,
        K=max_causal,
        maxite=maxite,
        eps=eps,
        ubound=ubound,
        cthres=coverage,
        eincre=eincre,
        minldthres=minldthres,
        maxldthres=maxldthres,
        varemax=varemax,
        varemin=varemin,
    )

    pips = pd.Series(data=zfile["PIP"].to_numpy(), index=zfile["SNPID"].to_numpy())
    cs_snps = []
    lead_snps = []
    for cs_i, sub_df in zfile.groupby("cs"):
        if cs_i == 0:
            continue
        cs_snps.append(sub_df["SNPID"].values.tolist())
        lead_snps.append(pips[pips.index.isin(sub_df["SNPID"].values)].idxmax())
    return CredibleSet(
        tool=Method.RSparsePro,
        n_cs=len(cs_snps),
        coverage=coverage,
        lead_snps=lead_snps,
        snps=cs_snps,
        cs_sizes=[len(i) for i in cs_snps],
        pips=pips,
        parameters=parameters,
    )
