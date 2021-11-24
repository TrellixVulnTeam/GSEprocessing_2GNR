from dataclasses import dataclass
import GEOparse
import pandas as pd
import scanpy as sc
import os
import sys


@dataclass
class RNASeqExpr:
    file_in: str

    def __post_init__(self):
        accessionID, gpl, _ = os.path.basename(file_in).split("-")
        self.accessionID = accessionID
        self.gpl = gpl

        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gse = gse

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

    def counts(self):
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
            col_rename[gsm.metadata["title"][0]] = name
        df = df.rename(columns=col_rename)
        return df

    def expr(self):
        adata = sc.AnnData(self.counts().T)
        # adata.var_names_make_unique()
        # 1e6 = counts per million (cpm) normalization
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata, base=2)

        df = adata.to_df().T
        return df

    def export(self):
        expr = f"{self.accessionID}-{self.gpl}-expr.txt"
        self.expr().to_csv(expr, sep="\t")


if __name__ == "__main__":
    file_in = sys.argv[1]
    my_gse = RNASeqExpr(file_in)
    my_gse.export()
