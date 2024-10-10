import logging
import os
import shutil

import pytest

from mafm.utils import io_in_tempdir

# Setup logging for tests
logging.basicConfig(level=logging.DEBUG)


def test_io_in_tempdir_creates_temp_dir():
    @io_in_tempdir()
    def test_func(temp_dir):
        assert os.path.exists(temp_dir)
        return temp_dir

    result = test_func()  # type: ignore
    assert not os.path.exists(
        result
    )  # Temp dir should be deleted after function execution


def test_io_in_tempdir_custom_parent_dir():
    custom_dir = "./custom_tmp"

    @io_in_tempdir(dir=custom_dir)
    def test_func(temp_dir):
        # assert os.path.dirname(temp_dir) == custom_dir
        return temp_dir

    result = test_func()  # type: ignore
    assert not os.path.exists(result)
    assert os.path.exists(custom_dir)  # Parent directory should still exist
    shutil.rmtree(custom_dir)  # Clean up manually


def test_io_in_tempdir_retains_dir_on_debug():
    logger = logging.getLogger("IO")
    original_level = logger.level
    logger.setLevel(logging.DEBUG)

    @io_in_tempdir()
    def test_func(temp_dir):
        return temp_dir

    result = test_func()  # type: ignore
    assert os.path.exists(result)  # Temp dir should be retained
    os.rmdir(result)  # Clean up manually

    logger.setLevel(original_level)  # Restore original logging level


def test_io_in_tempdir_handles_exceptions():
    @io_in_tempdir()
    def test_func(temp_dir):
        raise ValueError("Test exception")

    with pytest.raises(ValueError):
        test_func()  # type: ignore


def test_io_in_tempdir_passes_args_and_kwargs():
    @io_in_tempdir()
    def test_func(arg1, temp_dir, kwarg1=None):
        assert arg1 == "test_arg"
        assert os.path.exists(temp_dir)
        assert kwarg1 == "test_kwarg"
        return arg1, temp_dir, kwarg1

    result = test_func("test_arg", kwarg1="test_kwarg")  # type: ignore
    assert result[0] == "test_arg"
    assert result[2] == "test_kwarg"
    assert not os.path.exists(result[1])  # Temp dir should be deleted


@pytest.mark.parametrize("log_level", [logging.INFO, logging.WARNING, logging.ERROR])
def test_io_in_tempdir_removes_dir_above_info(log_level):
    logger = logging.getLogger("IO")
    original_level = logger.level
    logger.setLevel(log_level)

    @io_in_tempdir()
    def test_func(temp_dir):
        return temp_dir

    result = test_func()  # type: ignore
    assert not os.path.exists(result)  # Temp dir should be deleted for INFO and above

    logger.setLevel(original_level)  # Restore original logging level


def test_io_in_tempdir_cleanup_error_handling(mocker):
    # Mock shutil.rmtree to simulate a cleanup error
    mock_rmtree = mocker.patch(
        "shutil.rmtree", side_effect=OSError("Mocked cleanup error")
    )

    # Create a mock logger
    mock_logger = mocker.Mock()
    mock_logger.getEffectiveLevel.return_value = logging.INFO  # Set the log level
    mock_logger.warning = mocker.Mock()  # Add a mock warning method

    # Patch the getLogger function to return our mock logger
    mocker.patch("logging.getLogger", return_value=mock_logger)

    @io_in_tempdir()
    def test_func(temp_dir):
        return temp_dir

    result = test_func()  # type: ignore

    # Verify that rmtree was called
    mock_rmtree.assert_called_once()

    # Verify that a warning was logged
    mock_logger.warning.assert_called_once_with(mocker.ANY)

    # Optional: You can check the content of the warning message if needed
    warning_message = mock_logger.warning.call_args[0][0]
    assert "Failed to remove temporary directory" in warning_message
