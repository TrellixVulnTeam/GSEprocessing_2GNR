"""The following script takes a GSE identification and using GEOParse downloads the .soft.gz file. From this file, 
four files are created to initiate Boolean Implication Analysis network. These files are used in the supplementary 
gse_processing bash script to generate the network.
"""
from dataclasses import dataclass
import GEOparse
from GEOparse.GEOTypes import GSE
import pandas as pd
import sys
import re
import gzip
import shutil
import numpy as np


@dataclass
class NCBIGeo:
    """Takes a GSE Accession ID string and uses GEOParse module to pull .soft.gz file from NCBI Geo. The gse object is
    assigned as an attribute to the class. Five methods are available to create five files for export, four of which are
    required for the boolean implication analysis network.
    """

    accessionID: str

    def __post_init__(self):
        """Adds GEOparse gse class as attribute"""
        gse = GEOparse.get_GEO(geo=str(self.accessionID), silent=True)
        setattr(self, "gse", gse)

    def soft_file_txt(self):
        """Converts .soft.gz file into easily readable .txt file"""
        soft_gz = str(self.accessionID) + "_family.soft.gz"
        soft_txt = soft_gz[:-3] + ".txt"
        with gzip.open(soft_gz, "rb") as f_in, open(soft_txt, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    def make_expr(self, gpl, takeLog=False, export=False):
        """Pulls expression data from .soft file

        Args:
            gpl (str): a GEOparse gpl machine name
            takeLog (bool, optional): If True, function takes Log2 of all values. Defaults to False.
            export (bool, optional): If Ture, function exports to .txt file. Defaults to False.

        Returns:
            pandas.DataFrame: DataFrame of .soft file expression data
        """
        for name, gsm in self.gse.gsms.items():
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl.name:
                continue

            gsm_df = gsm.table.set_index("ID_REF")
            gsm_df.columns = [name]

            if takeLog:
                gsm_df[name] = np.where(gsm_df[name] > 0, np.log2(gsm_df[name]), -1)

            if "expr_df" not in locals():
                expr_df = gsm_df
            else:
                expr_df = expr_df.merge(gsm_df, left_index=True, right_index=True)

        if export:
            expr_filename = f"{self.accessionID}-{gpl.name}-expr.txt"
            expr_df.to_csv(expr_filename, sep="\t")
        return expr_df

    def make_idx(self, gpl, export=False):
        """Makes idx dataframe including binary expression information for Boolean Network

        Args:
            gpl (str): GEO GPL name
            export (bool, optional): If True, method exports dataframe to .txt. Defaults to False.

        Returns:
            pandas.DataFrame: DataFrame including idx information
        """
        pos = 0
        idx = {
            "Ptr": [],
            "ProbeID": [],
            "Name": [],
            "Description": [],
        }

        expr = f"{self.accessionID}-{gpl.name}-expr.txt"
        with open(expr, "rb") as f:
            for line in f:
                if pos == 0:
                    pos += len(line)
                else:
                    idx["Ptr"].append(pos)
                    pos += len(line)
                    split = line.decode("utf-8").split("\t")
                    idx["ProbeID"].append(split[0])
                    idx["Name"].append(split[1].split(":")[0])
                    idx["Description"].append(":".join(split[1].split(":")[1:]))

        idx_df = pd.DataFrame(idx).set_index("ProbeID")

        if export:
            idx_filename = f"{self.accessionID}-{gpl.name}-idx.txt"
            idx_df.to_csv(idx_filename, sep="\t")
        return idx_df

    def make_survival(self, gpl, export=False):
        to_drop = [
            "geo_accession",
            "status",
            "date",
            "protocol",
            "proccesing",
            "data_processing",
            "contact",
            "supplementary",
            "platform_id",
            "series_id",
        ]

        all_metadata = {}
        for name, gsm in self.gse.gsms.items():
            # confirm gsm is associated with input gpl
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl.name:
                continue

            # remove keys from metadata that aren't desired in survival
            metadata = gsm.metadata.copy()
            for key in gsm.metadata:
                for drop in to_drop:
                    if re.search(drop, key):
                        metadata.pop(key, None)
            all_metadata[name] = metadata
        all_metadata = pd.DataFrame(all_metadata).T

        for column in metadata.columns:
            # split columns with list values into seperate columns
            all_metadata[column] = all_metadata[column].apply(lambda x: "\t".join(x))
            to_merge = all_metadata[column].str.split("\t", expand=True)
            if len(to_merge.columns) > 1:
                col_names = [str(i + 1) for i in range(len(to_merge.columns))]
                col_names = [column + "_" + name for name in col_names]
                to_merge.columns = col_names
            else:
                to_merge.columns = [column]

            if "df" not in locals():
                df = to_merge
            else:
                df = df.merge(to_merge, left_index=True, right_index=True)
        df.index.name = "ArrayID"

        for column in df.columns:
            # rename columns with cell value label
            if df[column].str.contains(":").all():
                value = df[column].str.extract(r"(.*:)").iloc[0, 0][:-1]
                df = df.rename(columns={column: value})
                df[value] = df[value].str.extract(r".*: (.*)")

        if export:
            survival_filename = f"{self.accessionID}-{gpl.name}-survival.txt"
            df.to_csv(survival_filename, sep="\t")
        return df

    def make_ih(self, gpl, export=False):
        survival_df = self.make_survival(gpl).reset_index()
        ih_df = survival_df[["ArrayID", "title"]]
        ih_df.insert(1, "ArrayHeader", ih_df["ArrayID"])
        ih_df = ih_df.set_index("ArrayID")

        if export:
            ih_filename = f"{self.accessionID}-{gpl.name}-ih.txt"
            ih_df.to_csv(ih_filename, sep="\t")
        return ih_df

    def export_all(self):
        for _, gpl in self.gse.gpls.items():
            # self.soft_file_txt()
            self.make_expr(gpl, takeLog=True, export=True)
            self.make_idx(gpl, export=True)
            self.make_survival(gpl, export=True)
            self.make_ih(gpl, export=True)


if __name__ == "__main__":
    gse = sys.argv[1]
    ncbi = NCBIGeo(gse)
    ncbi.export_all()
