"""Functions for processing summary statistics data."""

import logging
from typing import Optional

import pandas as pd

from .constants import ColName

logger = logging.getLogger("Sumstats")


def get_significant_snps(
    df: pd.DataFrame,
    pvalue_threshold: float = 5e-8,
    use_most_sig_if_no_sig: bool = True,
) -> pd.DataFrame:
    """
    Retrieve significant SNPs from the input DataFrame based on a p-value threshold.

    If no SNPs meet the significance threshold and `use_most_sig_if_no_sig` is True,
    the function returns the SNP with the smallest p-value.

    Parameters
    ----------
    df : pd.DataFrame
        The input summary statistics containing SNP information.
    pvalue_threshold : float, optional
        The p-value threshold for significance, by default 5e-8.
    use_most_sig_if_no_sig : bool, optional
        Whether to return the most significant SNP if no SNP meets the threshold, by default True.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing significant SNPs, sorted by p-value in ascending order.

    Raises
    ------
    ValueError
        If no significant SNPs are found and `use_most_sig_if_no_sig` is False.
    KeyError
        If required columns are not present in the input DataFrame.

    Examples
    --------
    >>> data = {
    ...     'SNPID': ['rs1', 'rs2', 'rs3'],
    ...     'P': [1e-9, 0.05, 1e-7]
    ... }
    >>> df = pd.DataFrame(data)
    >>> significant_snps = get_significant_snps(df, pvalue_threshold=5e-8)
    >>> print(significant_snps)
        snp_id  p_value
    0    rs1   1.0e-09
    1    rs3   1.0e-07
    """
    required_columns = {ColName.P, ColName.SNPID}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise KeyError(
            f"The following required columns are missing from the DataFrame: {missing_columns}"
        )

    sig_df = df.loc[df[ColName.P] < pvalue_threshold].copy()

    if sig_df.empty:
        if use_most_sig_if_no_sig:
            min_pvalue = df[ColName.P].min()
            sig_df = df.loc[df[ColName.P] == min_pvalue].copy()
            if sig_df.empty:
                raise ValueError("The DataFrame is empty. No SNPs available to select.")
            logging.debug(
                f"Using the most significant SNP: {sig_df.iloc[0][ColName.SNPID]}"
            )
            logging.debug(f"p-value: {sig_df.iloc[0][ColName.P]}")
        else:
            raise ValueError("No significant SNPs found.")
    else:
        sig_df.sort_values(by=ColName.P, inplace=True)
        sig_df.reset_index(drop=True, inplace=True)

    return sig_df


def make_SNPID_unique(
    sumstat: pd.DataFrame,
    remove_duplicates: bool = True,
    col_chr: str = ColName.CHR,
    col_bp: str = ColName.BP,
    col_ea: str = ColName.EA,
    col_nea: str = ColName.NEA,
    col_p: str = ColName.P,
) -> pd.DataFrame:
    """
    Generate unique SNP identifiers to facilitate the combination of multiple summary statistics datasets.

    This function constructs a unique SNPID by concatenating chromosome, base-pair position,
    and sorted alleles (EA and NEA). This unique identifier allows for efficient merging of
    multiple summary statistics without the need for extensive duplicate comparisons.

    The unique SNPID format: "chr-bp-sortedEA-sortedNEA"

    Parameters
    ----------
    sumstat : pd.DataFrame
        The input summary statistics containing SNP information.
    remove_duplicates : bool, optional
        Whether to remove duplicated SNPs, keeping the one with the smallest p-value, by default True.
    col_chr : str, optional
        The column name for chromosome information, by default ColName.CHR.
    col_bp : str, optional
        The column name for base-pair position information, by default ColName.BP.
    col_ea : str, optional
        The column name for effect allele information, by default ColName.EA.
    col_nea : str, optional
        The column name for non-effect allele information, by default ColName.NEA.
    col_p : str, optional
        The column name for p-value information, by default

    Returns
    -------
    pd.DataFrame
        The summary statistics DataFrame with unique SNPIDs, suitable for merging with other datasets.

    Raises
    ------
    KeyError
        If required columns are missing from the input DataFrame.
    ValueError
        If the input DataFrame is empty or becomes empty after processing.

    Examples
    --------
    >>> data = {
    ...     'CHR': ['1', '1', '2'],
    ...     'BP': [12345, 12345, 67890],
    ...     'EA': ['A', 'A', 'G'],
    ...     'NEA': ['G', 'G', 'A'],
    ...     'rsID': ['rs1', 'rs2', 'rs3'],
    ...     'P': [1e-5, 1e-6, 1e-7]
    ... }
    >>> df = pd.DataFrame(data)
    >>> unique_df = make_SNPID_unique(df, replace_rsIDcol=True, remove_duplicates=True)
    >>> print(unique_df)
        SNPID   CHR BP  EA  NEA rsID    P
    0  1-12345-A-G    1  12345  A   G  rs1  1.0e-05
    1  2-67890-A-G    2  67890  G   A  rs3  1.0e-07
    """
    # col_chr = col_chr or ColName.CHR
    # col_bp = col_bp or ColName.BP
    # col_ea = col_ea or ColName.EA
    # col_nea = col_nea or ColName.NEA
    required_columns = {
        col_chr,
        col_bp,
        col_ea,
        col_nea,
    }
    missing_columns = required_columns - set(sumstat.columns)
    if missing_columns:
        raise KeyError(
            f"The following required columns are missing from the DataFrame: {missing_columns}"
        )

    if sumstat.empty:
        raise ValueError("The input DataFrame is empty.")

    df = sumstat.copy()

    # Sort alleles to ensure unique representation (EA <= NEA)
    allele_df = df[[col_ea, col_nea]].apply(
        lambda row: sorted([row[col_ea], row[col_nea]]), axis=1, result_type="expand"
    )
    allele_df.columns = [col_ea, col_nea]

    # Create unique SNPID
    df[ColName.SNPID] = (
        df[col_chr].astype(str)
        + "-"
        + df[col_bp].astype(str)
        + "-"
        + allele_df[col_ea]
        + "-"
        + allele_df[col_nea]
    )

    # move SNPID to the first column
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index(ColName.SNPID)))
    df = df[cols]

    n_duplicated = df.duplicated(subset=[ColName.SNPID]).sum()

    if remove_duplicates and n_duplicated > 0:
        logger.debug(f"Number of duplicated SNPs: {n_duplicated}")
        # Sort by p-value to keep the SNP with the smallest p-value
        df.sort_values(by=col_p, inplace=True)
        df.drop_duplicates(subset=[ColName.SNPID], keep="first", inplace=True)
        # Sort DataFrame by chromosome and base-pair position
        df.sort_values(by=[col_chr, col_bp], inplace=True)
        df.reset_index(drop=True, inplace=True)
    elif n_duplicated > 0 and not remove_duplicates:
        logger.warning(
            """Duplicated SNPs detected. To remove duplicates, set `remove_duplicates=True`.
            Change the Unique SNP identifier to make it unique."""
        )
        # Change the Unique SNP identifier to make it unique. add a number to the end of the SNP identifier
        #  for example, 1-12345-A-G to 1-12345-A-G-1, 1-12345-A-G-2, etc. no alteration to the original SNP identifier
        dup_tail = "-" + df.groupby(ColName.SNPID).cumcount().astype(str)
        dup_tail = dup_tail.str.replace("-0", "")
        df[ColName.SNPID] = df[ColName.SNPID] + dup_tail

    if df.empty:
        raise ValueError("The resulting DataFrame is empty after processing.")

    logging.debug("Unique SNPIDs have been successfully created.")
    logging.debug(f"Total unique SNPs: {len(df)}")

    return df
