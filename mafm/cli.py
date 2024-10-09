"""Command line interface."""

import typer

app = typer.Typer()


def main():
    """MAFM: A user-friendly tool for fine-mapping."""
    typer.echo("mafm")
    typer.echo("=" * len("mafm"))
    typer.echo("multi-ancestry fine-mapping pipeline")


if __name__ == "__main__":
    app(main)
