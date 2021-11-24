import click

from src.NanoString import NanoString
from src.NCBIGeo import NCBIGeo
from src.tar2rcc import make_rcc_dir, gunzip_files

# from src.RNASeqExpr import RNASeqExpr

__author__ = "Oliver Tucher"


@click.group()
def main() -> None:
    """
    Simple CLI for processing GSE datasets
    """
    pass


@main.command()
@click.argument("output_dir")
def nano(output_dir: str) -> None:
    """Creates nano files

    Args:
        output_dir (str): output directory for nano files
    """
    output_file = NanoString(output_dir).export("code_summary")
    click.echo(f"Output file created: {output_file}")


@main.command()
@click.argument("accession_id")
def ncbi(accession_id: str) -> None:
    """Creates NCBI files

    Args:
        accession_id (str): id of accession
    """
    NCBIGeo(accession_id).export_all()


@main.command()
@click.argument("tar_file")
def tar2rcc(tar_file: str) -> None:
    """Parse a tar file into directory containing RCC files

    Args:
        tar_file (str): tar zipped file
    """
    rcc_dir = make_rcc_dir(tar_file)
    gunzip_files(rcc_dir)


@main.command()
@click.argument("input_file")
def rna(input_file: str) -> None:
    """Generate RNA sequence counts

    Args:
        input_file (str): input RNA file
    """
    pass
    # RNASeqExpr(input_file).export()


if __name__ == "__main__":
    main()
