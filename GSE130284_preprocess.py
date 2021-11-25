import pandas as pd
import GEOparse


class GSE130284:
    file_in: str

    def __post_init__(self):
        self.accessionID = "GSE130284"
        self.gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gpl = self.gse.gpls.values()[0]

    def preprocess(self):
        """Particular preprocess for GSE130284

        Returns:
            pd.DataFrame: Raw Counts dataframe
        """
        df = pd.read_exceL(self.file_in)
        df = df.rename(columns={"Id": "Name", "KNalone_2": "NKalone_2"})
        df = df.set_index("Name")
        df = df.iloc[:, :12]

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
        df.columns = col_names

    def counts(self):
        df = self.preprocess()
        if "ProbeID" not in df.columns:
            # merge with reference file containing gene names and ENSG
            ensg_file = "~/GSEprocessing/ensg.txt"
            ensg_df = pd.read_csv(ensg_file, sep="\t", header=0)
            ensg_df = ensg_df.drop("ENST", axis=1)

            df = df.merge(ensg_df, on="Name")
            df = df.set_index(["ProbeID", "Name"])

        col_rename = {}
        for name, gsm in self.gse.gsms.items():
            # map gsm name to sample title
            col_rename[gsm.metadata["title"][0]] = name
        df = df.rename(columns=col_rename)

        return df

    def export(self):
        df = self.counts()
        filename = f"{self.accessionID}-{self.gpl}-counts.txt"
        df.to_csv(filename)
