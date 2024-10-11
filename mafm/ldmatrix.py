"""Functions for reading and converting lower triangle matrices."""

import logging
from typing import Dict, Union

import numpy as np
import pandas as pd

from mafm.constants import ColName
from mafm.sumstats import make_SNPID_unique, munge_bp, munge_chr

logger = logging.getLogger("Sumstats")


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
                list(map(float, line.strip().split(delimiter))) for line in file if line.strip()
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
    outdf = make_SNPID_unique(outdf, col_ea=ColName.A1, col_nea=ColName.A2, remove_duplicates=False)
    return outdf


def load_ld(
    ld_path: str, map_path: str, delimiter: str = "\t"
) -> Dict[str, Union[pd.DataFrame, np.ndarray]]:
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

    Returns
    -------
    Dict[str, Union[pd.DataFrame, np.ndarray]]
        Dictionary containing the Variant IDs and the LD matrix.

    Raises
    ------
    ValueError
        If the number of variants in the map file does not match the number of rows in the LD matrix.
    """
    ld_df = load_ld_matrix(ld_path, delimiter)
    map_df = load_ld_map(map_path, delimiter)
    if ld_df.shape[0] != map_df.shape[0]:
        raise ValueError(
            "The number of variants in the map file does not match the number of rows in the LD matrix."
        )

    return {"map": map_df, "r": ld_df}
