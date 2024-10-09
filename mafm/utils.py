"""Functions and decorators for common tasks in Python programming."""

import logging
import shutil
import tempfile
from functools import wraps
from typing import Callable, TypeVar

logger = logging.getLogger("Utils")


# Type variable for decorator
F = TypeVar("F", bound=Callable[..., any])  # type: ignore


def io_in_tempdir(dir: str = "./tmp") -> Callable[[F], F]:
    """
    Create a temporary directory for I/O operations during function execution.

    This decorator creates a temporary directory before executing the decorated function and
    provides the path to this directory via the `temp_dir` keyword argument. After the function
    execution, the temporary directory is removed based on the logging level:
    - If the logging level is set to `INFO` or higher, the temporary directory is deleted.
    - If the logging level is lower than `INFO` (e.g., `DEBUG`), the directory is retained for inspection.

    Parameters
    ----------
    dir : str, optional
        The parent directory where the temporary directory will be created, by default "./tmp".

    Returns
    -------
    Callable[[F], F]
        A decorator that manages a temporary directory for the decorated function.

    Raises
    ------
    OSError
        If the temporary directory cannot be created.

    Examples
    --------
    ```python
    @io_in_tempdir(dir="./temporary")
    def process_data(temp_dir: str, data: str) -> None:
        # Perform I/O operations using temp_dir
        with open(f"{temp_dir}/data.txt", "w") as file:
            file.write(data)

    process_data(data="Sample data")
    ```
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs) -> any:  # type: ignore
            temp_dir = tempfile.mkdtemp(dir=dir)
            logger = logging.getLogger("IO")
            logger.debug(f"Created temporary directory: {temp_dir}")

            try:
                # Inject temp_dir into the function's keyword arguments
                result = func(*args, temp_dir=temp_dir, **kwargs)
            except Exception as e:
                logger.error(f"An error occurred in function '{func.__name__}': {e}")
                raise
            else:
                # Determine whether to remove the temporary directory based on the logging level
                if logger.getEffectiveLevel() >= logging.INFO:
                    try:
                        shutil.rmtree(temp_dir)
                        logger.debug(f"Removed temporary directory: {temp_dir}")
                    except Exception as cleanup_error:
                        logger.warning(
                            f"Failed to remove temporary directory '{temp_dir}': {cleanup_error}"
                        )
                else:
                    logger.debug(
                        f"Retaining temporary directory '{temp_dir}' for inspection due to logging level."
                    )
                return result

        return wrapper  # type: ignore

    return decorator
