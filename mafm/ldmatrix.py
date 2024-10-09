"""
This module contains functions for reading and converting lower triangle matrices.
"""

import warnings
from typing import List

import numpy as np


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


def lower_triangle_to_symmetric(file_path: str, delimiter: str = "\t") -> np.ndarray:
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
    # Read the lower triangle matrix from the file
    lower_triangle = read_lower_triangle(file_path, delimiter)

    # Determine the size of the square matrix
    n = lower_triangle.shape[0]

    # Create the symmetric matrix
    symmetric_matrix = lower_triangle + lower_triangle.T

    # Fill the diagonal with 1
    np.fill_diagonal(symmetric_matrix, 1)

    # convert to float16
    symmetric_matrix = symmetric_matrix.astype(np.float16)

    return symmetric_matrix
