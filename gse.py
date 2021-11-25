import click
import os

from src.GEO2Hegemon import GEO2Hegemon
from src.Counts2Expr import Counts2Expr

__author__ = "Oliver Tucher"


@click.group()
def main() -> None:
    """
    CLI for processing GSE datasets
    """
    pass


@main.command()
@click.argument("accession_id")
@click.option("--output", help="destination for files")
@click.option("--counts", help="use cleaned counts file")
def geo2hegemon(accession_id: str, output: str = None, counts: str = None) -> None:
    """Create hegemon files from NCBI GEO Accession ID

    Args:
        accessionID (str): NCBI GEO AccessionID to process
        output (str, optional): directory to save created files. Defaults to None
        counts (str, optional): a raw counts file if required. Defaults to None.
    """
    if output == None:
        output = "./" + accession_id

    if not os.path.exists(output):
        os.mkdir(output)

    if counts != None:
        Counts2Expr(counts).export()

    os.chdir(output)

    GEO2Hegemon(accession_id).export_all()
    click.echo(f"Output filed in: {output}")


if __name__ == "__main__":
    main()
