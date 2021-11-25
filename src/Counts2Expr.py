from dataclasses import dataclass
import os
import pandas as pd
import scanpy as sc


@dataclass
class Counts2Expr:
    file_in: str

    def __post_init__(self):
        if ".txt" not in self.file_in:
            raise ValueError(f"{self.file_in} must be a .txt file")

        if "GSE" not in self.file_in or "GPL" not in self.file_in:
            raise ValueError(
                "Counts file must be formated as 'GSEXXX-GPLXXX-counts.txt"
            )

        gse, gpl, _ = os.path.split(self.file_in)[1].split("-")
        self.accessionID = gse
        self.gpl = gpl

    def expr(self):
        df = pd.read_csv(self.file_in, sep="\t")
        df = df.set_index(["ProbeID", "Name"])

        adata = sc.AnnData(df.T)
        # adata.var_names_make_unique()
        # 1e6 = counts per million (cpm) normalization
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata, base=2)

        return adata.to_df().T

    def export(self):
        expr = f"{self.accessionID}-{self.gpl}-expr.txt"
        self.expr().to_csv(expr, sep="\t")
