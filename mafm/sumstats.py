"""Functions for processing summary statistics data."""

import gzip
import logging
from typing import Optional

import pandas as pd

from .constants import ColName, ColRange, ColType

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
    ...     'P': [1e-9, 0.05, 1e-8]
    ... }
    >>> df = pd.DataFrame(data)
    >>> significant_snps = get_significant_snps(df, pvalue_threshold=5e-8)
    >>> print(significant_snps)
        SNPID  P
    0    rs1   1.0e-09
    1    rs3   1.0e-08
    """
    required_columns = {ColName.P, ColName.SNPID}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise KeyError(
            f"The following required columns are missing from the DataFrame: {missing_columns}"
        )

    sig_df = df.loc[df[ColName.P] <= pvalue_threshold].copy()

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
        if col_p in df.columns:
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

    logging.debug("Unique SNPIDs have been successfully created.")
    logging.debug(f"Total unique SNPs: {len(df)}")

    return df


def check_colnames(df: pd.DataFrame) -> pd.DataFrame:
    """
    Check column names in the DataFrame and fill missing columns with None.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to check for column names.

    Returns
    -------
    pd.DataFrame
        DataFrame with all required columns, filling missing ones with None.
    """
    outdf: pd.DataFrame = df.copy()
    for col in ColName.sumstat_cols:
        if col not in outdf.columns:
            outdf[col] = None
    return outdf[ColName.sumstat_cols]


def check_mandatory_cols(df: pd.DataFrame) -> None:
    """
    Check if the DataFrame contains all mandatory columns.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to check for mandatory columns.

    Raises
    ------
    ValueError
        If any mandatory columns are missing.
    """
    outdf = df.copy()
    missing_cols = set(ColName.mandatory_cols) - set(outdf.columns)
    if missing_cols:
        raise ValueError(f"Missing mandatory columns: {missing_cols}")
    return None


def rm_col_allna(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove columns from the DataFrame that are entirely NA.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame from which to remove columns.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns that are entirely NA removed.
    """
    outdf = df.copy()
    outdf = outdf.replace("", None)
    for col in outdf.columns:
        if outdf[col].isnull().all():
            logger.debug(f"Remove column {col} because it is all NA.")
            outdf.drop(col, axis=1, inplace=True)
    return outdf


def munge(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge the summary statistics DataFrame by performing a series of transformations.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing summary statistics.

    Returns
    -------
    pd.DataFrame
        The munged DataFrame with necessary transformations applied.

    Raises
    ------
    ValueError
        If any mandatory columns are missing.
    """
    check_mandatory_cols(df)
    outdf = df.copy()
    outdf = rm_col_allna(outdf)
    outdf = munge_chr(outdf)
    outdf = munge_bp(outdf)
    outdf = munge_allele(outdf)
    outdf = make_SNPID_unique(outdf)
    outdf = munge_pvalue(outdf)
    outdf = outdf.sort_values(by=[ColName.CHR, ColName.BP])
    outdf = munge_beta(outdf)
    outdf = munge_se(outdf)
    outdf = munge_eaf(outdf)
    outdf[ColName.MAF] = outdf[ColName.EAF]
    outdf = munge_maf(outdf)
    if ColName.RSID in outdf.columns:
        outdf = munge_rsid(outdf)
    outdf = check_colnames(outdf)
    return outdf


def munge_rsid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge rsID column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with rsID column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged rsID column.
    """
    outdf = df.copy()
    outdf[ColName.RSID] = outdf[ColName.RSID].astype(ColType.RSID)
    return outdf


def munge_chr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge chromosome column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with chromosome column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged chromosome column.
    """
    pre_n = df.shape[0]
    outdf = df[df[ColName.CHR].notnull()].copy()
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(str)
    outdf[ColName.CHR] = outdf[ColName.CHR].str.replace("chr", "")
    outdf[ColName.CHR] = outdf[ColName.CHR].replace(["X", "x"], 23)
    outdf[ColName.CHR] = pd.to_numeric(outdf[ColName.CHR], errors="coerce")
    outdf = outdf[outdf[ColName.CHR].notnull()]
    outdf = outdf[
        (outdf[ColName.CHR] >= ColRange.CHR_MIN)
        & (outdf[ColName.CHR] <= ColRange.CHR_MAX)
    ]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid chromosome.")
    outdf[ColName.CHR] = outdf[ColName.CHR].astype(ColType.CHR)
    return outdf


def munge_bp(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge position column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with position column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged position column.
    """
    pre_n = df.shape[0]
    outdf = df[df[ColName.BP].notnull()].copy()
    outdf[ColName.BP] = pd.to_numeric(outdf[ColName.BP], errors="coerce")
    outdf = outdf[outdf[ColName.BP].notnull()]
    outdf = outdf[
        (outdf[ColName.BP] > ColRange.BP_MIN) & (outdf[ColName.BP] < ColRange.BP_MAX)
    ]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid position.")
    outdf[ColName.BP] = outdf[ColName.BP].astype(ColType.BP)
    return outdf


def munge_allele(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge allele columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with allele columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged allele columns.
    """
    outdf = df.copy()
    for col in [ColName.EA, ColName.NEA]:
        pre_n = outdf.shape[0]
        outdf = outdf[outdf[col].notnull()]
        outdf[col] = outdf[col].astype(str).str.upper()
        outdf = outdf[outdf[col].str.match(r"^[ACGT]+$")]
        after_n = outdf.shape[0]
        logger.debug(f"Remove {pre_n - after_n} rows because of invalid {col}.")
    outdf = outdf[outdf[ColName.EA] != outdf[ColName.NEA]]
    return outdf


def munge_pvalue(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge p-value column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with p-value column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged p-value column.
    """
    outdf = df.copy()
    pre_n = outdf.shape[0]
    outdf[ColName.P] = pd.to_numeric(outdf[ColName.P], errors="coerce")
    outdf = outdf[outdf[ColName.P].notnull()]
    outdf = outdf[
        (outdf[ColName.P] > ColRange.P_MIN) & (outdf[ColName.P] < ColRange.P_MAX)
    ]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid p-value.")
    outdf[ColName.P] = outdf[ColName.P].astype(ColType.P)
    return outdf


def munge_beta(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge beta column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with beta column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged beta column.
    """
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.BETA] = pd.to_numeric(outdf[ColName.BETA], errors="coerce")
    outdf = outdf[outdf[ColName.BETA].notnull()]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid beta.")
    outdf[ColName.BETA] = outdf[ColName.BETA].astype(ColType.BETA)
    return outdf


def munge_se(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge standard error column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with standard error column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged standard error column.
    """
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.SE] = pd.to_numeric(outdf[ColName.SE], errors="coerce")
    outdf = outdf[outdf[ColName.SE].notnull()]
    outdf = outdf[outdf[ColName.SE] > ColRange.SE_MIN]
    after_n = outdf.shape[0]
    logger.debug(f"Remove {pre_n - after_n} rows because of invalid standard error.")
    outdf[ColName.SE] = outdf[ColName.SE].astype(ColType.SE)
    return outdf


def munge_eaf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge effect allele frequency column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with effect allele frequency column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged effect allele frequency column.
    """
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.EAF] = pd.to_numeric(outdf[ColName.EAF], errors="coerce")
    outdf = outdf[outdf[ColName.EAF].notnull()]
    outdf = outdf[
        (outdf[ColName.EAF] >= ColRange.EAF_MIN)
        & (outdf[ColName.EAF] <= ColRange.EAF_MAX)
    ]
    after_n = outdf.shape[0]
    logger.debug(
        f"Remove {pre_n - after_n} rows because of invalid effect allele frequency."
    )
    outdf[ColName.EAF] = outdf[ColName.EAF].astype(ColType.EAF)
    return outdf


def munge_maf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Munge minor allele frequency column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with minor allele frequency column.

    Returns
    -------
    pd.DataFrame
        DataFrame with munged minor allele frequency column.
    """
    pre_n = df.shape[0]
    outdf = df.copy()
    outdf[ColName.MAF] = pd.to_numeric(outdf[ColName.MAF], errors="coerce")
    outdf = outdf[outdf[ColName.MAF].notnull()]
    outdf[ColName.MAF] = outdf[ColName.MAF].apply(lambda x: 1 - x if x > 0.5 else x)
    outdf = outdf[
        (outdf[ColName.MAF] >= ColRange.MAF_MIN)
        & (outdf[ColName.MAF] <= ColRange.MAF_MAX)
    ]
    after_n = outdf.shape[0]
    logger.debug(
        f"Remove {pre_n - after_n} rows because of invalid minor allele frequency."
    )
    outdf[ColName.MAF] = outdf[ColName.MAF].astype(ColType.MAF)
    return outdf


def load_sumstats(
    filename: str,
    sep: Optional[str] = None,
    nrows: Optional[int] = None,
    skiprows: int = 0,
    comment: Optional[str] = None,
    gzipped: Optional[bool] = None,
) -> pd.DataFrame:
    """
    Load summary statistics from a file.

    Parameters
    ----------
    filename : str
        The path to the file containing the summary statistics.
        The header must contain the column names: CHR, BP, EA, NEA, EAF, BETA, SE, P.
    sep : str, optional
        The delimiter to use. If None, the delimiter is inferred from the file.
    nrows : int, optional
        Number of rows to read. If None, all rows are read.
    skiprows : int, default 0
        Number of lines to skip at the start of the file.
    comment : str, optional
        Character to split comments in the file.
    gzipped : bool, optional
        Whether the file is gzipped. If None, it is inferred from the file extension.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the loaded summary statistics.

    Notes
    -----
    The function infers the delimiter if not provided and handles gzipped files.
    """
    # determine whether the file is gzipped
    if gzipped is None:
        gzipped = filename.endswith("gz")

    # read the first line of the file to determine the separator
    if sep is None:
        if gzipped:
            f = gzip.open(filename, "rt")

        else:
            f = open(filename, "rt")
        if skiprows > 0:
            for _ in range(skiprows):
                f.readline()
        line = f.readline()
        f.close()
        if "\t" in line:
            sep = "\t"
        elif "," in line:
            sep = ","
        else:
            sep = " "
    logger.debug(f"File {filename} is gzipped: {gzipped}")
    logger.debug(f"Separator is {sep}")
    logger.debug(f"loading data from {filename}")
    # determine the separator, automatically if not specified
    sumstats = pd.read_csv(
        filename,
        sep=sep,
        nrows=nrows,
        skiprows=skiprows,
        comment=comment,
        compression="gzip" if gzipped else None,
    )
    sumstats = munge(sumstats)
    return sumstats