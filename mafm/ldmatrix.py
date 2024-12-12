"""Functions for reading and converting lower triangle matrices."""

import logging

import numpy as np
import pandas as pd

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
        self.__check_length()

    def __repr__(self):
        """Return a string representation of the LDMatrix object."""
        return f"LDMatrix(map={self.map.shape}, r={self.r.shape})"

    def __check_length(self):
        """Check if the number of rows in the map file matches the number of rows in the LD matrix."""
        if len(self.map) != len(self.r):
            raise ValueError("The number of rows in the map file does not match the number of rows in the LD matrix.")

    def copy(self):
        """Return a copy of the LDMatrix object."""
        return LDMatrix(self.map.copy(), self.r.copy())


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
            rows = [list(map(float, line.strip().split(delimiter))) for line in file if line.strip()]
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    if not rows:
        raise ValueError("The input file is empty.")

    n = len(rows)
    lower_triangle = np.zeros((n, n))

    for i, row in enumerate(rows):
        if len(row) != i + 1:
            raise ValueError(f"Invalid number of elements in row {i + 1}. Expected {i + 1}, got {len(row)}.")
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
    if file_path.endswith(".npz"):
        ld_file_key = np.load(file_path).files[0]
        return np.load(file_path)[ld_file_key].astype(np.float32)
    lower_triangle = read_lower_triangle(file_path, delimiter)

    # Create the symmetric matrix
    symmetric_matrix = lower_triangle + lower_triangle.T

    # Fill the diagonal with 1
    np.fill_diagonal(symmetric_matrix, 1)

    # convert to float32
    symmetric_matrix = symmetric_matrix.astype(np.float32)
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
    # TODO: use REF/ALT instead of A1/A2
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
    outdf.reset_index(drop=True, inplace=True)
    # TODO: check if allele frequency is available
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


def load_ld(ld_path: str, map_path: str, delimiter: str = "\t", if_sort_alleles: bool = True) -> LDMatrix:
    """
    Read LD matrices and Variant IDs from files. Pair each matrix with its corresponding Variant IDs.

    Parameters
    ----------
    ld_path : str
        Path to the input text file containing the lower triangle matrix.
    map_path : str
        Path to the input text file containing the Variant IDs.
    delimiter : str, optional
        TODO: Support for npz files.
        TODO: Support for plink bin4 format.
        TODO: Support for ldstore bcor.
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
