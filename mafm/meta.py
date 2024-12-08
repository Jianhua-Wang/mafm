"""Meta analysis of multi-ancestry gwas data."""

import logging
import textwrap
from typing import Optional

import numpy as np
import pandas as pd

from mafm.constants import ColName
from mafm.ldmatrix import LDMatrix
from mafm.locus import Locus
from mafm.sumstats import munge, sort_alleles
from mafm.utils import io_in_tempdir, tool_manager

logger = logging.getLogger("MAFM")


@io_in_tempdir("./tmp/metal")
def run_metal(
    inputs: list[Locus],
    temp_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Run METAL for meta-analysis.

    TODO: replace metal with python implementation.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.
    temp_dir : Optional[str], optional
        Temporary directory, by default None.

    Returns
    -------
    pd.DataFrame
        Meta-analysis summary statistics.

    """
    logger.info(f"Running METAL for {len(inputs)} datasets")
    logger.info(f"Temp dir: {temp_dir}")

    metal_describe = """
    MARKER   SNPID
    DEFAULT  %i
    ALLELE   EA   NEA
    FREQ     EAF
    EFFECT   BETA
    STDERR   SE

    PROCESS %s
    """
    all_metal_describe = []
    for i, input in enumerate(inputs):
        sumstat = input.original_sumstats[
            [ColName.SNPID, ColName.EA, ColName.NEA, ColName.EAF, ColName.BETA, ColName.SE]
        ]
        sample_size = input.sample_size
        sumstat.to_csv(f"{temp_dir}/study{i}.sumstats", sep="\t", index=False)
        all_metal_describe.append(metal_describe % (sample_size, f"{temp_dir}/study{i}.sumstats"))
    all_metal_describe = "\n".join(all_metal_describe)
    out_file = f"{temp_dir}/metal_out"
    meta_script = f"""
    SCHEME   STDERR

    {all_metal_describe}

    OUTFILE {out_file} .txt
    ANALYZE"""
    with open(f"{temp_dir}/metal.metal", "w") as f:
        f.write(textwrap.dedent(meta_script))

    # Run METAL
    cmd = [f"{temp_dir}/metal.metal"]
    tool_manager.run_tool("metal", cmd, output_file_path=f"{out_file}1.txt")

    # Load METAL output
    sumstats = pd.read_csv(f"{out_file}1.txt", sep="\t")
    sumstats.columns = [ColName.SNPID, ColName.EA, ColName.NEA, ColName.BETA, ColName.SE, ColName.P, "Direction"]
    sumstats = sort_alleles(sumstats)
    sumstats.set_index(ColName.SNPID, inplace=True)

    return sumstats[[ColName.BETA, ColName.SE, ColName.P]]


def meta_sumstats(
    inputs: list[Locus],
    tool: str,
) -> pd.DataFrame:
    """
    Perform meta-analysis of summary statistics.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.
    tool : str
        Meta-analysis tool.
        You can choose from 'metal', 'metasoft'.

    Returns
    -------
    pd.DataFrame
        Meta-analysis summary statistics.

    """
    # concatenate input sumstats and meta EAF
    all_sumstat = []
    ith = 1
    for input in inputs:
        sumstat = input.original_sumstats
        sumstat = sumstat[[ColName.SNPID, ColName.CHR, ColName.BP, ColName.EA, ColName.NEA]]
        all_sumstat.append(sumstat)
        ith += 1
    all_sumstat = pd.concat(all_sumstat, axis=0)
    all_sumstat.drop_duplicates(subset=[ColName.SNPID], inplace=True)
    all_sumstat.set_index(ColName.SNPID, inplace=True)
    # meta EAF
    eaf_df = pd.DataFrame(index=all_sumstat.index)
    n_sum = sum([input.sample_size for input in inputs])
    weights = [input.sample_size / n_sum for input in inputs]
    for i, input in enumerate(inputs):
        df = input.original_sumstats.copy()
        df.set_index(ColName.SNPID, inplace=True)
        eaf_df[f"EAF_{i}"] = df[ColName.EAF] * weights[i]

    all_sumstat["EAF_META"] = eaf_df.sum(axis=1)

    if tool == "metal":
        meta_res = run_metal(inputs)
    elif tool == "metasoft":
        raise NotImplementedError("METASOFT is not implemented yet.")
    else:
        raise ValueError(f"Unsupported meta-analysis tool: {tool}")

    all_sumstat["BETA_META"] = meta_res[ColName.BETA]
    all_sumstat["SE_META"] = meta_res[ColName.SE]
    all_sumstat["P_META"] = meta_res[ColName.P]
    all_sumstat = all_sumstat[
        [ColName.CHR, ColName.BP, ColName.EA, ColName.NEA, "EAF_META", "BETA_META", "SE_META", "P_META"]
    ]
    all_sumstat.columns = [
        ColName.CHR,
        ColName.BP,
        ColName.EA,
        ColName.NEA,
        ColName.EAF,
        ColName.BETA,
        ColName.SE,
        ColName.P,
    ]
    all_sumstat.reset_index(inplace=True)
    return munge(all_sumstat)


def meta_lds(
    inputs: list[Locus],
) -> LDMatrix:
    """
    Perform meta-analysis of LD matrices.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.

    Returns
    -------
    LDMatrix
        Meta-analysis LD matrix.
    """
    # 1. Get unique variants across all studies
    # TODO: meta allele frequency of LD reference, if exists
    variant_dfs = [input.ld.map for input in inputs]
    ld_matrices = [input.ld.r for input in inputs]
    sample_sizes = [input.sample_size for input in inputs]

    merged_variants = pd.concat(variant_dfs, ignore_index=True)
    merged_variants.drop_duplicates(subset=["SNPID"], inplace=True)
    merged_variants.sort_values(["CHR", "BP"], inplace=True)
    merged_variants.reset_index(drop=True, inplace=True)
    all_variants = merged_variants["SNPID"].values
    variant_to_index = {snp: idx for idx, snp in enumerate(all_variants)}
    n_variants = len(all_variants)

    # 2. Initialize arrays using numpy operations
    merged_ld = np.zeros((n_variants, n_variants))
    weight_matrix = np.zeros((n_variants, n_variants))

    # 3. Process each study
    for ld_mat, variants_df, sample_size in zip(ld_matrices, variant_dfs, sample_sizes):
        # coverte float16 to float32, to avoid overflow
        ld_mat = ld_mat.astype(np.float32)

        # Get indices in the master matrix
        study_snps = variants_df["SNPID"].values
        study_indices = np.array([variant_to_index[snp] for snp in study_snps])

        # Create index meshgrid for faster indexing
        idx_i, idx_j = np.meshgrid(study_indices, study_indices)

        # Update matrices using vectorized operations
        merged_ld[idx_i, idx_j] += ld_mat * sample_size
        weight_matrix[idx_i, idx_j] += sample_size

    # 4. Compute weighted average
    mask = weight_matrix != 0
    merged_ld[mask] /= weight_matrix[mask]

    # 5. Prepare output variant information
    # Get complete variant information from the first occurrence of each variant
    # merged_variants = pd.concat(variant_dfs).drop_duplicates(subset="SNPID", keep="first")
    # merged_variants = merged_variants.set_index("SNPID").loc[all_variants].reset_index()
    return LDMatrix(merged_variants, merged_ld.astype(np.float16))


def meta_all(
    inputs: list[Locus],
    tool: str,
) -> Locus:
    """
    Perform meta-analysis of summary statistics and LD matrices.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.
    tool : str
        Meta-analysis tool.
        You can choose from 'metal', 'metasoft'.

    Returns
    -------
    Locus
        Meta-analysis result.

    """
    meta_sumstat = meta_sumstats(inputs, tool)
    meta_ld = meta_lds(inputs)
    sample_size = sum([input.sample_size for input in inputs])
    popu = set()
    for input in inputs:
        for pop in input.popu.split(","):
            popu.add(pop)
    popu = ",".join(sorted(popu))
    cohort = set()
    for input in inputs:
        for cohort_name in input.cohort.split(","):
            cohort.add(cohort_name)
    cohort = ",".join(sorted(cohort))

    return Locus(popu, cohort, sample_size, sumstats=meta_sumstat, ld=meta_ld, if_intersect=True)


def meta_by_population(
    inputs: list[Locus],
    tool: str,
) -> dict[str, Locus]:
    """
    Perform meta-analysis of summary statistics and LD matrices within each population.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.
    tool : str
        Meta-analysis tool.
        You can choose from 'metal', 'metasoft'.

    Returns
    -------
    Locus
        Meta-analysis result.

    """
    meta_popu = {}
    for input in inputs:
        popu = input.popu
        if popu not in meta_popu:
            meta_popu[popu] = [input]
        else:
            meta_popu[popu].append(input)

    for popu in meta_popu:
        if len(meta_popu[popu]) > 1:
            meta_popu[popu] = meta_all(meta_popu[popu], tool)
        else:
            meta_popu[popu] = meta_popu[popu][0]
    return meta_popu


def meta(
    inputs: list[Locus],
    tool: str,
    meta_method: str = "meta_all",
) -> list[Locus]:
    """
    Perform meta-analysis of summary statistics and LD matrices.

    Parameters
    ----------
    inputs : list[Locus]
        List of input data.
    tool : str
        Meta-analysis tool.
    meta_method : str, optional
        Meta-analysis method, by default "meta_all"
        Options: "meta_all", "meta_by_population", "no_meta".

    Returns
    -------
    Locus
        Meta-analysis result.

    """
    if meta_method == "meta_all":
        return [meta_all(inputs, tool)]
    elif meta_method == "meta_by_population":
        res = meta_by_population(inputs, tool)
        return [res[popu] for popu in res]
    elif meta_method == "no_meta":
        return inputs
    else:
        raise ValueError(f"Unsupported meta-analysis method: {meta_method}")
