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
    Simple CLI for processing GSE datasets
    """
    pass


@main.command()
@click.argument("input_dir")
@click.argument("output_dir")
def rcc2hegemon(input_dir: str, output_dir: str) -> None:
    """Creates NanoString files from directory of .RCC

    Args:
        input_dir (str): output directory for nano files
    """
    output_dir = NanoStringRCC(input_dir, output_dir).export_all()
    click.echo(f"Output file created: {output_dir}")


@main.command()
@click.argument("accession_id")
@click.argument("output_dir")
@click.option("--counts_file", help="counts_file")
def geo2hegemon(accession_id: str, output_dir: str, counts_file: str = None) -> None:
    """Create hegemon files from NCBI GEO Accession ID

    Args:
        accessionID (str): NCBI GEO AccessionID to process
        output_dir (str): directory to save created files
        counts_file (str, optional): a raw counts file if required. Defaults to None.
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)

    if counts_file != None:
        Counts2Expr(accession_id, counts_file, output_dir)

    GEO2Hegemon(accession_id).export_all()


@main.command()
@click.argument("tar_file")
def tar2rcc(tar_file: str) -> None:
    """Parse a tar file into directory containing RCC files

    Args:
        tar_file (str): tar zipped file
    """
    rcc_dir = make_rcc_dir(tar_file)
    gunzip_files(rcc_dir)


if __name__ == "__main__":
    main()
