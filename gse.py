#!/bin/env python3

import argparse
import os

from src.GEO2Hegemon import GEO2Hegemon
from src.Counts2Expr import Counts2Expr

__author__ = "Oliver Tucher"

parser = argparse.ArgumentParser(description="CLI For NCBI GEO Accenssion Processing")
parser.add_argument(
    "accessionID",
    metavar="GSE ID",
    type=str,
    help="NCBI GEO Accession ID (GSE ID) to parse",
)
parser.add_argument(
    "-o",
    "--output_dir",
    metavar="Output Directory",
    help="Directory location to file parsed data",
)
parser.add_argument(
    "-c",
    "--counts",
    metavar="Raw 'GSE*-GPL*-counts.txt file",
    help="optional raw counts.txt file if required",
)
args = parser.parse_args()

print(f"Parsing {args.accessionID}")


def geo2hegemon(accessionID: str, output_dir: str = None, counts: str = None) -> None:
    """Create hegemon files from NCBI GEO Accession ID

    Args:
        accessionID (str): NCBI GEO AccessionID to process
        output (str, optional): directory to save created files. Defaults to None
        counts (str, optional): a raw counts file if required. Defaults to None.
    """
    if output_dir == None:
        output_dir = "./" + accessionID

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if counts != None:
        Counts2Expr(counts).export()

    os.chdir(output_dir)

    GEO2Hegemon(accessionID).export_all()


geo2hegemon(args.accessionID, output_dir=args.output_dir, counts=args.counts)
