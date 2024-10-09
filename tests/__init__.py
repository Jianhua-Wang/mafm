"""Unit test package for mafm."""

"""Top-level package for easy_finemap."""

import logging

from rich.logging import RichHandler

__author__ = """Jianhua Wang"""
__email__ = "jianhua.mert@gmail.com"
__version__ = "0.0.1"


logging.basicConfig(
    level=logging.WARNING,
    format="%(name)s - %(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
