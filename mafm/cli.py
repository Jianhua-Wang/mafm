"""Console script for mafm."""

import typer

app = typer.Typer()


def main():
    """Main entrypoint."""
    typer.echo("mafm")
    typer.echo("=" * len("mafm"))
    typer.echo("multi-ancestry fine-mapping pipeline")


if __name__ == "__main__":
    app(main)
