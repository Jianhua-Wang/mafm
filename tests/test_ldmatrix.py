"""Unit tests for the ldmatrix module."""

import numpy as np
import pandas as pd
import pytest

from mafm.ldmatrix import load_ld, load_ld_map, load_ld_matrix, read_lower_triangle


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
        dtype=np.float32,
    )

    result = load_ld_matrix(file_path)
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
        load_ld_matrix(file_path)


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


def test_load_ld_map_valid_file(tmp_path):
    """
    Test loading a valid LD map file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected DataFrame.
    """
    file_content = "CHR\tBP\tA1\tA2\n1\t1000\tA\tG\n1\t2000\tC\tT\n2\t3000\tT\tC\n"
    file_path = tmp_path / "valid_map.txt"
    file_path.write_text(file_content)

    expected_df = pd.DataFrame(
        {
            "SNPID": ["1-1000-A-G", "1-2000-C-T", "2-3000-C-T"],
            "CHR": [1, 1, 2],
            "BP": [1000, 2000, 3000],
            "A1": ["A", "C", "T"],
            "A2": ["G", "T", "C"],
        }
    )

    result = load_ld_map(file_path)
    pd.testing.assert_frame_equal(result, expected_df)


def test_load_ld_map_missing_columns(tmp_path):
    """
    Test loading an LD map file with missing columns.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    file_content = "CHR\tBP\tA1\n1\t1000\tA\n1\t2000\tC\n2\t3000\tT\n"
    file_path = tmp_path / "missing_columns_map.txt"
    file_path.write_text(file_content)

    with pytest.raises(
        ValueError, match=r"Missing columns in the input file: \['A2'\]"
    ):
        load_ld_map(file_path)


def test_load_ld_map_invalid_alleles(tmp_path):
    """
    Test loading an LD map file with invalid alleles.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected DataFrame.
    """
    file_content = "CHR\tBP\tA1\tA2\n1\t1000\tA\tG\n1\t2000\tC\tT\n2\t3000\tT\tN\n"
    file_path = tmp_path / "invalid_alleles_map.txt"
    file_path.write_text(file_content)

    expected_df = pd.DataFrame(
        {
            "SNPID": ["1-1000-A-G", "1-2000-C-T"],
            "CHR": [1, 1],
            "BP": [1000, 2000],
            "A1": ["A", "C"],
            "A2": ["G", "T"],
        }
    )

    result = load_ld_map(file_path)
    pd.testing.assert_frame_equal(result, expected_df)


def test_load_ld_valid_files(tmp_path):
    """
    Test loading LD matrices and Variant IDs from valid files.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not return the expected dictionary.
    """
    ld_content = "1\n0.1\t1\n0.2\t0.4\t1\n"
    ld_path = tmp_path / "lower_triangle.txt"
    ld_path.write_text(ld_content)

    map_content = "CHR\tBP\tA1\tA2\n1\t1000\tA\tG\n1\t2000\tC\tT\n2\t3000\tT\tC\n"
    map_path = tmp_path / "valid_map.txt"
    map_path.write_text(map_content)

    expected_ld_matrix = np.array(
        [
            [1, 0.1, -0.2],
            [0.1, 1, -0.4],
            [-0.2, -0.4, 1],
        ],
        dtype=np.float32,
    )

    expected_map_df = pd.DataFrame(
        {
            "SNPID": ["1-1000-A-G", "1-2000-C-T", "2-3000-C-T"],
            "CHR": [1, 1, 2],
            "BP": [1000, 2000, 3000],
            "A1": ["A", "C", "C"],
            "A2": ["G", "T", "T"],
        }
    )

    result = load_ld(ld_path, map_path)
    pd.testing.assert_frame_equal(result.map, expected_map_df)  # type: ignore
    assert np.array_equal(
        result.r, expected_ld_matrix
    ), "The LD matrix is not as expected."


def test_load_ld_mismatched_files(tmp_path):
    """
    Test loading LD matrices and Variant IDs from files with mismatched number of variants.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the function does not raise the expected ValueError.
    """
    ld_content = "1\n0.1\t1\n0.2\t0.4\t1\n0.3\t0.5\t0.6\t1\n"
    ld_path = tmp_path / "lower_triangle.txt"
    ld_path.write_text(ld_content)

    map_content = "CHR\tBP\tA1\tA2\n1\t1000\tA\tG\n1\t2000\tC\tT\n2\t3000\tT\tC\n"
    map_path = tmp_path / "valid_map.txt"
    map_path.write_text(map_content)

    with pytest.raises(
        ValueError,
        match="The number of variants in the map file does not match the number of rows in the LD matrix.",
    ):
        load_ld(ld_path, map_path)
