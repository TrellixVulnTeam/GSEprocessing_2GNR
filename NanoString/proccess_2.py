from dataclasses import dataclass
import GEOparse
from GEOparse.GEOTypes import GSE
import pandas as pd
import os
import sys
import re
import gzip
import shutil
import numpy as np


@dataclass
class NCBIGeo:
    accessionID: str

    def __post_init__(self):
        gse = GEOparse.get_GEO(geo=str(self.accessionID), silent=True)
        setattr(self, "gse", gse)

    def soft_file_txt(self):
        soft_gz = str(self.accessionID) + "_family.soft.gz"
        soft_txt = soft_gz[:-3] + ".txt"
        with gzip.open(soft_gz, "rb") as f_in, open(soft_txt, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    def make_expr(self, gpl, takeLog=False, export=False):
        expr_df = pd.DataFrame()
        for name, gsm in self.gse.gsms.items():
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl.name:
                continue

            gsm_df = gsm.table.set_index("ID_REF")
            gsm_df.columns = [name]

            if takeLog:
                gsm_df[name] = np.where(gsm_df[name] > 0, np.log2(gsm_df[name]), -1)
            if expr_df.empty:
                expr_df = gsm_df
            else:
                expr_df = expr_df.merge(gsm_df, left_index=True, right_index=True)

        if export:
            expr_filename = f"{self.accessionID}-{gpl.name}-expr.txt"
            expr_df.to_csv(expr_filename, sep="\t")
        return expr_df

    def make_idx(self, gpl, export=False):
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
            metadata = gsm.metadata.copy()
            for key in gsm.metadata:
                for drop in to_drop:
                    if re.search(drop, key):
                        metadata.pop(key, None)
            all_metadata[name] = metadata

        df = pd.DataFrame(all_metadata).T
        survival_df = pd.DataFrame()
        for column in df.columns:
            df[column] = df[column].apply(lambda x: "\t".join(x))
            to_merge = df[column].str.split("\t", expand=True)
            if len(to_merge.columns) > 1:
                col_names = []
                for col in to_merge.columns:
                    value = to_merge[col].str.extract(r"(.*:)").iloc[0, 0][:-1]
                    col_names.append(value)
                to_merge.columns = col_names
            else:
                to_merge.columns = [column]
            if survival_df.empty:
                survival_df = to_merge
            else:
                survival_df = survival_df.merge(
                    to_merge, left_index=True, right_index=True
                )
        survival_df.index.name = "ArrayID"

        if export:
            survival_filename = f"{self.accessionID}-{gpl.name}-survival.txt"
            survival_df.to_csv(survival_filename, sep="\t")
        return survival_df

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
            self.soft_file_txt()
            self.make_expr(gpl, takeLog=True, export=True)
            self.make_idx(gpl, export=True)
            self.make_survival(gpl, export=True)
            self.make_ih(gpl, export=True)


if __name__ == "__main__":
    gse = sys.argv[1]
    ncbi = NCBIGeo(gse)
    ncbi.export_all()
