from dataclasses import dataclass
import GEOparse
import pandas as pd
import scanpy as sc


@dataclass
class Counts2Expr:
    accessionID: str
    file_in: str
    output_dir: str

    def __post_init__(self):
        if ".txt" not in self.file_in:
            raise ValueError(f"{self.file_in} must be a .txt file")

        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        gpl = gse.gpls.items()[0][1]
        self.gpl = gpl

        df = pd.read_csv(self.file_in, sep="\t")
        df = df.set_index(["Probe ID", "Name"])

        adata = sc.AnnData(df.T)
        # adata.var_names_make_unique()
        # 1e6 = counts per million (cpm) normalization
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata, base=2)

        return adata.to_df().T

    def export(self):
        expr = f"{self.accessionID}-{self.gpl}-expr.txt"
        self.expr().to_csv(expr, sep="\t")
