from dataclasses import dataclass
import GEOparse
from GEOparse.GEOTypes import GSE
import pandas as pd

# import sqlite3
# import csv
import os
import sys

# import math
import re

# import tempfile
# import statistics
# import argparse
# import pymongo
import gzip
import shutil

# import numpy as np
# import io
# import scanpy as sc
# import glob
# from anndata import AnnData


@dataclass
class NCBIGeo:
    accessionID: str

    def __post_init__(self):
        path = os.path.join(".", self.accessionID)
        gse = GEOparse.get_GEO(geo=str(self.accessionID), destdir=path, silent=True)
        setattr(self, "gse", gse)

        soft_gz = str(self.accessionID) + "_family.soft.gz"
        soft_gz = os.path.join(path, soft_gz)
        soft_txt = soft_gz[:-3] + ".txt"
        with gzip.open(soft_gz, "rb") as f_in, open(soft_txt, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    def m_analyze(self):
        for gsm_name, gsm in self.gse.gsms.items():
            print("Name: ", gsm_name)
            print(
                "Metadata:",
            )


if __name__ == "__main__":
    gse = sys.argv[1]
    ncbi = NCBIGeo(gse)
    ncbi.m_analyze()
