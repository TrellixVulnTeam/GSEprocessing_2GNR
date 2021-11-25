import click
import os

from src.NanoStringRCC import NanoStringRCC
from src.GEO2Hegemon import GEO2Hegemon
from src.Counts2Expr import Counts2Expr
from src.tar2rcc import make_rcc_dir, gunzip_files

__author__ = "Oliver Tucher"


@click.group()
def main() -> None:
    """
    CLI for processing GSE datasets and NanoString RCC data
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


@main.command()
@click.argument("tar_file")
def tar2rcc(tar_file: str) -> None:
    """Parse a tar file into directory containing RCC files

    Args:
        tar_file (str): tar zipped file
    """
    rcc_dir = make_rcc_dir(tar_file)
    gunzip_files(rcc_dir)
    click.echo(f"Output filed: {rcc_dir}")


if __name__ == "__main__":
    main()


@main.command()
@click.argument("input_dir")
@click.argument("output_dir")
def rcc2hegemon(input_dir: str, output_dir: str) -> None:
    """Creates NanoString files from directory of .RCC

    Args:
        input_dir (str): output directory for nano files
    """

    NanoStringRCC(input_dir, output_dir).export_all()
    click.echo(f"Output file in: {output_dir}")
