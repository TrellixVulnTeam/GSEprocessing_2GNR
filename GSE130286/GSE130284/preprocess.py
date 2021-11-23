from dataclasses import dataclass
import GEOparse
import pandas as pd
from anndata import AnnData
import scanpy as sc
import sys


@dataclass
class Preprocess:
    accessionID: str
    file_in: str

    def __post_init__(self):
        """Adds GEOparse gse class as attribute"""
        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        setattr(self, "gse", gse)

    def rename_columns(self, df):
        col_names = []
        for col in df.columns:
            col = col.split("_")
            col = " donor ".join(col)
            col = col[:2] + " " + col[2:] + "_RNA-seq"
            if "IL" in col:
                col = col.replace("IL", "IL-")
            if "INB" in col:
                col = col.replace("INB16", "CTV-1")
            col_names.append(col)
        return col_names

    def preprocess(self, gpl):
        df = pd.read_csv(self.file_in)
        df = df.set_index("Name")

        df.columns = self.rename_columns(df)

        if "ProbeID" not in df.columns:
            ensg_file = "HomoSapiens_ENST,ProbeID,Name.txt"
            ensg_df = pd.read_csv(ensg_file, sep="\t", header=0)
            ensg_df = ensg_df.drop("ENST", axis=1)

            df = df.merge(ensg_df, on="Name")
            df = df.set_index(["ProbeID", "Name"])

        col_rename = {}
        for name, gsm in self.gse.gsms.items():
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl:
                continue

            col_rename[gsm.metadata["title"][0]] = name
        df = df.rename(columns=col_rename)

        file_out = f"{self.accessionID}-{gpl}-counts.txt"
        df.to_csv(file_out, sep="\t")

        return df

    def make_expr(self, gpl):
        expr = self.preprocess(gpl)
        expr = expr.drop(["ProbeID", "Name"], axis=1)

        print("Normalizing")
        adata = AnnData(expr.T)
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata, base=2)

        norm_df = pd.DataFrame(adata.X)
        norm_df = norm_df.T

        # norm_df.insert(0, 'ProbeID', list(df['ProbeID']))
        # norm_df.insert(1, 'Name', list(df['Name']))
        # norm_df.columns = list(df.columns)

        print(norm_df)


if __name__ == "__main__":
    accessionID = sys.argv[1]
    file_in = sys.argv[2]
    my_gse = Preprocess(accessionID, file_in)
    gpl = list(my_gse.gse.gpls.keys())[0]
    my_gse.make_expr(gpl)
