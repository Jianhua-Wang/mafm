"""Meta analysis of multi-ancestry gwas data."""

import logging
import os
from multiprocessing import Pool
import numpy as np
import pandas as pd
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeRemainingColumn
from scipy import stats

from mafm.constants import ColName
from mafm.ldmatrix import LDMatrix
from mafm.locus import Locus, LocusSet, intersect_sumstat_ld, load_locus_set
from mafm.sumstats import munge

logger = logging.getLogger("META")


def meta_sumstats(
    inputs: LocusSet,
) -> pd.DataFrame:
    """
    Perform fixed effect meta-analysis.

    Parameters
    ----------
    inputs : LocusSet
        List of input data.

    Returns
    -------
    pd.DataFrame
        Meta-analysis summary statistics.
    """
    # Merge all dataframes on SNPID
    merged_df = inputs.loci[0].original_sumstats[[ColName.SNPID]].copy()
    n_sum = sum([input.sample_size for input in inputs.loci])
    eaf_weights = [input.sample_size / n_sum for input in inputs.loci]
    for i, df in enumerate(inputs.loci):
        df = df.sumstats[[ColName.SNPID, ColName.BETA, ColName.SE, ColName.EAF]].copy()
        df.rename(columns={ColName.BETA: f"BETA_{i}", ColName.SE: f"SE_{i}", ColName.EAF: f"EAF_{i}"}, inplace=True)
        merged_df = pd.merge(merged_df, df, on=ColName.SNPID, how="outer", suffixes=("", f"_{i}"))

    # Calculate weights (inverse of variance)
    for i in range(len(inputs.loci)):
        merged_df[f"weight_{i}"] = 1 / (merged_df[f"SE_{i}"] ** 2)
        # merged_df[f"EAF_{i}"] = merged_df[f"EAF_{i}"] * eaf_weights[i]

    merged_df.fillna(0, inplace=True)

    # Calculate meta-analysis beta
    beta_numerator = sum(merged_df[f"BETA_{i}"] * merged_df[f"weight_{i}"] for i in range(len(inputs.loci)))
    weight_sum = sum(merged_df[f"weight_{i}"] for i in range(len(inputs.loci)))
    meta_beta = beta_numerator / weight_sum

    # Calculate meta-analysis SE
    meta_se = np.sqrt(1 / weight_sum)

    # Calculate meta-analysis Z-score and p-value
    meta_z = meta_beta / meta_se
    meta_p = 2 * stats.norm.sf(abs(meta_z))

    # Calculate meta-analysis EAF
    meta_eaf = sum(merged_df[f"EAF_{i}"] * eaf_weights[i] for i in range(len(inputs.loci)))

    # Create output dataframe
    output_df = pd.DataFrame(
        {
            ColName.SNPID: merged_df[ColName.SNPID],
            ColName.BETA: meta_beta,
            ColName.SE: meta_se,
            ColName.P: meta_p,
            ColName.EAF: meta_eaf,
        }
    )
    output_df[[ColName.CHR, ColName.BP, ColName.EA, ColName.NEA]] = merged_df[ColName.SNPID].str.split(
        "-", expand=True
    )[[0, 1, 2, 3]]
    return munge(output_df)


def meta_lds(
    inputs: LocusSet,
) -> LDMatrix:
    """
    Perform meta-analysis of LD matrices.

    Parameters
    ----------
    inputs : LocusSet
        List of input data.

    Returns
    -------
    LDMatrix
        Meta-analysis LD matrix.
    """
    # Get unique variants across all studies
    variant_dfs = [input.ld.map for input in inputs.loci]
    ld_matrices = [input.ld.r for input in inputs.loci]
    sample_sizes = [input.sample_size for input in inputs.loci]

    # Concatenate all variants
    merged_variants = pd.concat(variant_dfs, ignore_index=True)
    merged_variants.drop_duplicates(subset=[ColName.SNPID], inplace=True)
    merged_variants.sort_values([ColName.CHR, ColName.BP], inplace=True)
    merged_variants.reset_index(drop=True, inplace=True)
    # meta allele frequency of LD reference, if exists
    if all("AF2" in variant_df.columns for variant_df in variant_dfs):
        n_sum = sum([input.sample_size for input in inputs.loci])
        weights = [input.sample_size / n_sum for input in inputs.loci]
        af_df = merged_variants[[ColName.SNPID]].copy()
        af_df.set_index(ColName.SNPID, inplace=True)
        for i, variant_df in enumerate(variant_dfs):
            df = variant_df.copy()
            df.set_index(ColName.SNPID, inplace=True)
            af_df[f"AF2_{i}"] = df["AF2"] * weights[i]
        af_df.fillna(0, inplace=True)
        af_df["AF2_meta"] = af_df.sum(axis=1)
        merged_variants["AF2"] = merged_variants[ColName.SNPID].map(af_df["AF2_meta"])
    all_variants = merged_variants[ColName.SNPID].values
    variant_to_index = {snp: idx for idx, snp in enumerate(all_variants)}
    n_variants = len(all_variants)

    # Initialize arrays using numpy operations
    merged_ld = np.zeros((n_variants, n_variants))
    weight_matrix = np.zeros((n_variants, n_variants))

    # Process each study
    for ld_mat, variants_df, sample_size in zip(ld_matrices, variant_dfs, sample_sizes):
        # coverte float16 to float32, to avoid overflow
        # ld_mat = ld_mat.astype(np.float32)

        # Get indices in the master matrix
        study_snps = variants_df["SNPID"].values
        study_indices = np.array([variant_to_index[snp] for snp in study_snps])

        # Create index meshgrid for faster indexing
        idx_i, idx_j = np.meshgrid(study_indices, study_indices)

        # Update matrices using vectorized operations
        merged_ld[idx_i, idx_j] += ld_mat * sample_size
        weight_matrix[idx_i, idx_j] += sample_size

    # Compute weighted average
    mask = weight_matrix != 0
    merged_ld[mask] /= weight_matrix[mask]

    return LDMatrix(merged_variants, merged_ld.astype(np.float32))


def meta_all(
    inputs: LocusSet,
) -> Locus:
    """
    Perform meta-analysis of summary statistics and LD matrices.

    Parameters
    ----------
    inputs : LocusSet
        List of input data.

    Returns
    -------
    Locus
        Meta-analysis result.

    """
    meta_sumstat = meta_sumstats(inputs)
    meta_ld = meta_lds(inputs)
    sample_size = sum([input.sample_size for input in inputs.loci])
    popu = set()
    for input in inputs.loci:
        for pop in input.popu.split(","):
            popu.add(pop)
    popu = "+".join(sorted(popu))
    cohort = set()
    for input in inputs.loci:
        for cohort_name in input.cohort.split(","):
            cohort.add(cohort_name)
    cohort = "+".join(sorted(cohort))

    return Locus(popu, cohort, sample_size, sumstats=meta_sumstat, ld=meta_ld, if_intersect=True)


def meta_by_population(
    inputs: LocusSet,
) -> dict[str, Locus]:
    """
    Perform meta-analysis of summary statistics and LD matrices within each population.

    Parameters
    ----------
    inputs : LocusSet
        List of input data.

    Returns
    -------
    Locus
        Meta-analysis result.

    """
    meta_popu = {}
    for input in inputs.loci:
        popu = input.popu
        if popu not in meta_popu:
            meta_popu[popu] = [input]
        else:
            meta_popu[popu].append(input)

    for popu in meta_popu:
        if len(meta_popu[popu]) > 1:
            meta_popu[popu] = meta_all(LocusSet(meta_popu[popu]))
        else:
            meta_popu[popu] = intersect_sumstat_ld(meta_popu[popu][0])
    return meta_popu


def meta(
    inputs: LocusSet,
    meta_method: str = "meta_all",
) -> LocusSet:
    """
    Perform meta-analysis of summary statistics and LD matrices.

    Parameters
    ----------
    inputs : LocusSet
        List of input data.
    meta_method : str, optional
        Meta-analysis method, by default "meta_all"
        Options: "meta_all", "meta_by_population", "no_meta".

    Returns
    -------
    LocusSet
        Meta-analysis result.

    """
    if meta_method == "meta_all":
        return LocusSet([meta_all(inputs)])
    elif meta_method == "meta_by_population":
        res = meta_by_population(inputs)
        return LocusSet([res[popu] for popu in res])
    elif meta_method == "no_meta":
        return LocusSet([intersect_sumstat_ld(i) for i in inputs.loci])
    else:
        raise ValueError(f"Unsupported meta-analysis method: {meta_method}")


def meta_locus(args):
    """
    Process a single locus for meta-analysis.

    Args
    ----
        args: A tuple containing (locus_id, locus_info, outdir, meta_method).
            locus_id: The ID of the locus.
            locus_info: DataFrame containing locus information.
            outdir: Output directory path.
            meta_method: Method for meta-analysis.

    Returns
    -------
        list: A list of results containing processed locus information.
    """
    locus_id, locus_info, outdir, meta_method = args
    results = []
    locus_set = load_locus_set(locus_info)
    locus_set = meta(locus_set, meta_method)
    out_dir = f"{outdir}/{locus_id}"
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    for locus in locus_set.loci:
        out_prefix = f"{out_dir}/{locus.prefix}"
        locus.sumstats.to_csv(f"{out_prefix}.sumstats.gz", sep="\t", index=False, compression="gzip")
        np.savez_compressed(f"{out_prefix}.ld.npz", ld=locus.ld.r.astype(np.float16))
        locus.ld.map.to_csv(f"{out_prefix}.ldmap.gz", sep="\t", index=False, compression="gzip")
        chrom, start, end = locus.chrom, locus.start, locus.end
        results.append(
            [
                chrom,
                start,
                end,
                locus.popu,
                locus.sample_size,
                locus.cohort,
                out_prefix,
                f'chr{chrom}_{start}_{end}',
            ]
        )
    return results


def meta_loci(
    inputs: str,
    outdir: str,
    threads: int = 1,
    meta_method: str = "meta_all",
):
    """
    Meta-analysis of summary statistics and LD matrices.

    Parameters
    ----------
    inputs : str
        Input files.
    outdir : str
        Output directory.
    threads : int, optional
        Number of threads, by default 1.
    meta_method : str, optional
        Meta-analysis method, by default "meta_all".

    Returns
    -------
    None
    """
    loci_info = pd.read_csv(inputs, sep="\t")
    new_loci_info = pd.DataFrame(columns=loci_info.columns)

    # Group loci by locus_id
    grouped_loci = list(loci_info.groupby("locus_id"))
    total_loci = len(grouped_loci)
    os.makedirs(outdir, exist_ok=True)

    # Create process pool and process loci in parallel with progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Meta-analysing...", total=total_loci)

        with Pool(threads) as pool:
            args = [(locus_id, locus_info, outdir, meta_method) for locus_id, locus_info in grouped_loci]
            for result in pool.imap_unordered(meta_locus, args):
                # Update results
                for i, res in enumerate(result):
                    new_loci_info.loc[len(new_loci_info)] = res
                # Update progress
                progress.advance(task)

    new_loci_info.to_csv(f"{outdir}/loci_info.txt", sep="\t", index=False)
