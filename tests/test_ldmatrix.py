import os
import tempfile

import numpy as np
import pytest

from mafm.ldmatrix import lower_triangle_to_symmetric, read_lower_triangle


@pytest.fixture
def sample_lower_triangle_file():
    content = "0.1\n0.2\t0.3\n0.4\t0.5\t0.6\n"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(content)
    yield tmp.name
    os.unlink(tmp.name)


def test_read_lower_triangle(sample_lower_triangle_file):
    result = read_lower_triangle(sample_lower_triangle_file)
    expected = np.array([[0.1, 0.0, 0.0], [0.2, 0.3, 0.0], [0.4, 0.5, 0.6]])
    np.testing.assert_array_almost_equal(result, expected)


def test_read_lower_triangle_custom_delimiter(sample_lower_triangle_file):
    content = "0.1\n0.2;0.3\n0.4;0.5;0.6"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(content)
    result = read_lower_triangle(tmp.name, delimiter=";")
    expected = np.array([[0.1, 0.0, 0.0], [0.2, 0.3, 0.0], [0.4, 0.5, 0.6]])
    np.testing.assert_array_almost_equal(result, expected)
    os.unlink(tmp.name)


def test_read_lower_triangle_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_lower_triangle("non_existent_file.txt")


def test_read_lower_triangle_empty_file():
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        pass
    with pytest.raises(ValueError, match="The input file is empty."):
        read_lower_triangle(tmp.name)
    os.unlink(tmp.name)


def test_read_lower_triangle_invalid_data():
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write("0.1\n0.2\t0.3\n0.4\t0.5\t0.6\t0.7\n")
    with pytest.raises(ValueError, match="Invalid number of elements in row"):
        read_lower_triangle(tmp.name)
    os.unlink(tmp.name)


def test_lower_triangle_to_symmetric(sample_lower_triangle_file):
    result = lower_triangle_to_symmetric(sample_lower_triangle_file)
    expected = np.array(
        [[1.0, 0.2, 0.4], [0.2, 1.0, 0.5], [0.4, 0.5, 1.0]],
        dtype=np.float16,
    )
    np.testing.assert_array_almost_equal(result, expected)


def test_lower_triangle_to_symmetric_custom_delimiter():
    content = "0.1\n0.2;0.3\n0.4;0.5;0.6"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(content)
    result = lower_triangle_to_symmetric(tmp.name, delimiter=";")
    expected = np.array(
        [[1.0, 0.2, 0.4], [0.2, 1.0, 0.5], [0.4, 0.5, 1.0]],
        dtype=np.float16,
    )
    np.testing.assert_array_almost_equal(result, expected)
    os.unlink(tmp.name)


def test_lower_triangle_to_symmetric_file_not_found():
    with pytest.raises(FileNotFoundError):
        lower_triangle_to_symmetric("non_existent_file.txt")


def test_lower_triangle_to_symmetric_empty_file():
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        pass
    with pytest.raises(ValueError, match="The input file is empty."):
        lower_triangle_to_symmetric(tmp.name)
    os.unlink(tmp.name)


def test_lower_triangle_to_symmetric_invalid_data():
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write("0.1\n0.2\t0.3\n0.4\t0.5\t0.6\t0.7\n")
    with pytest.raises(ValueError, match="Invalid number of elements in row"):
        lower_triangle_to_symmetric(tmp.name)
    os.unlink(tmp.name)
