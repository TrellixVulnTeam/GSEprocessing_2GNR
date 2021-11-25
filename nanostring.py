import click

from src.NanoStringSample import NanoStringSample
from src.NanoStringRCC import NanoStringRCC
from src.tar2rcc import make_rcc_dir, gunzip_files

__author__ = "Oliver Tucher"


@click.group()
def main() -> None:
    """CLI for processing NanoString RCC data"""
    pass


@main.command()
@click.argument("input_file")
def tar2rcc(input_file: str) -> None:
    """Extract all from .tar file and gunzip any .gz files

    Args:
        input_file (str): tarfile
    """
    rcc_dir = make_rcc_dir(input_file)
    gunzip_files(rcc_dir)
    click.echo(f"Output filed: {rcc_dir}")


@main.command()
@click.argument("input", help="directory of .RCC files")
@click.argument("output", help="directory to file hegemon files")
def rcc2hegemon(input: str, output: str) -> None:
    """Create hegemon files from directory of RCC files

    Args:
        input (str): Directory containing RCC files
        output (str): Directory where hegemon files will be filed
    """
    NanoStringRCC(input, output).export_all()
    click.echo(f"Output file in: {output}")


if __name__ == "__main__":
    main()
