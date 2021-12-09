#!/usr/bin/env python3

import argparse
import os

from src.MakeHegemon import MakeHegemon

__author__ = "Oliver Tucher"

parser = argparse.ArgumentParser(description="CLI to make '-idx.txt' file file")
parser.add_argument(
    "expr",
    type=str,
    help="Expression file to parse",
)
parser.add_argument(
    "-o",
    "--output_dir",
    metavar="Output Directory",
    help="Directory location to file parsed data",
)
args = parser.parse_args()
print(f"Creating IDX file from {args.expr}")


def make_idx(expr_file: str, output_dir: str = None) -> None:
    """Create a '-idx.txt' file from an '-expr.txt' file

    Args:
        expr_file (str): expression file to parse
        output_dir (str, optional): directory to file output. Defaults to None.
    """
    if output_dir == None:
        output_dir = os.path.join(os.getcwd(), expr_file)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    os.chdir(output_dir)

    idx_file = MakeHegemon().idx(expr_file=expr_file)
    idx_export = expr_file[:-8] + "idx.txt"
    idx_file.to_csv(idx_export, sep="\t")

    print(f"{idx_export} file created.")


make_idx(args.expr, output_dir=args.output_dir)
