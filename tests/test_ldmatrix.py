"""Unit tests for the ldmatrix module."""

import numpy as np
import pytest

from mafm.ldmatrix import lower_triangle_to_symmetric, read_lower_triangle


def test_read_lower_triangle_valid_file(tmp_path):
    """
    Test reading a valid lower triangle matrix from a file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected lower triangle matrix.
    """
    file_content = "1\n2\t3\n4\t5\t6\n"
    file_path = tmp_path / "lower_triangle.txt"
    file_path.write_text(file_content)

    expected_matrix = np.array([[1, 0, 0], [2, 3, 0], [4, 5, 6]], dtype=float)

    result = read_lower_triangle(file_path)
    assert np.array_equal(
        result, expected_matrix
    ), "The lower triangle matrix is not as expected."


def test_read_lower_triangle_invalid_file(tmp_path):
    """
    Test reading an invalid lower triangle matrix from a file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    file_content = "1\n2\t3\t4\n"
    file_path = tmp_path / "invalid_lower_triangle.txt"
    file_path.write_text(file_content)

    with pytest.raises(
        ValueError, match="Invalid number of elements in row 2. Expected 2, got 3."
    ):
        read_lower_triangle(file_path)


def test_read_lower_triangle_file_not_found():
    """
    Test reading a lower triangle matrix from a non-existent file.

    Raises
    ------
    AssertionError
        If the function does not raise the expected FileNotFoundError.
    """
    with pytest.raises(
        FileNotFoundError, match="The file 'non_existent_file.txt' does not exist."
    ):
        read_lower_triangle("non_existent_file.txt")


def test_lower_triangle_to_symmetric_valid_file(tmp_path):
    """
    Test converting a valid lower triangle matrix to a symmetric matrix.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected symmetric matrix.
    """
    file_content = "1\n0.1\t1\n0.2\t0.4\t1\n0.3\t0.5\t0.6\t1\n"
    file_path = tmp_path / "lower_triangle.txt"
    file_path.write_text(file_content)

    expected_matrix = np.array(
        [
            [1, 0.1, 0.2, 0.3],
            [0.1, 1, 0.4, 0.5],
            [0.2, 0.4, 1, 0.6],
            [0.3, 0.5, 0.6, 1],
        ],
        dtype=np.float16,
    )

    result = lower_triangle_to_symmetric(file_path)
    assert np.array_equal(
        result, expected_matrix
    ), "The symmetric matrix is not as expected."


def test_lower_triangle_to_symmetric_invalid_file(tmp_path):
    """
    Test converting an invalid lower triangle matrix to a symmetric matrix.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    file_content = "1\n2\t3\t4\n"
    file_path = tmp_path / "invalid_lower_triangle.txt"
    file_path.write_text(file_content)

    with pytest.raises(
        ValueError, match="Invalid number of elements in row 2. Expected 2, got 3."
    ):
        lower_triangle_to_symmetric(file_path)


def test_read_lower_triangle_empty_file(tmp_path):
    """
    Test reading an empty lower triangle matrix from a file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    file_path = tmp_path / "empty_lower_triangle.txt"
    file_path.write_text("")

    with pytest.raises(ValueError, match="The input file is empty."):
        read_lower_triangle(file_path)


def test_read_lower_triangle_custom_delimiter(tmp_path):
    """
    Test reading a lower triangle matrix from a file with a custom delimiter.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected lower triangle matrix.
    """
    file_content = "1\n2,3\n4,5,6\n"
    file_path = tmp_path / "lower_triangle_custom_delimiter.txt"
    file_path.write_text(file_content)

    expected_matrix = np.array([[1, 0, 0], [2, 3, 0], [4, 5, 6]], dtype=float)

    result = read_lower_triangle(file_path, delimiter=",")
    assert np.array_equal(
        result, expected_matrix
    ), "The lower triangle matrix is not as expected."


def test_read_lower_triangle_invalid_number_of_elements(tmp_path):
    """
    Test reading a lower triangle matrix with an invalid number of elements in a row.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    file_content = "1\n2\t3\t4\n"
    file_path = tmp_path / "invalid_number_of_elements.txt"
    file_path.write_text(file_content)

    with pytest.raises(
        ValueError, match="Invalid number of elements in row 2. Expected 2, got 3."
    ):
        read_lower_triangle(file_path)
