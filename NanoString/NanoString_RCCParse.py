"""Parse a directory of raw .RCC files (12 per NanoString read). The directory is converted toa NanoString class
where each of the 12 files are converted to a NanoStringSample class and parsed into one DataFrame per html tag.
The NanoString class then compiles the individual attributes of each NanoStringSample into one dataframe per overall
attribute

Returns:
    NanoString: A class that includes all 12 NanoString samples and relevant attributes that compile all sample info
    into compiled dataframes. Methods include raw_count and counts normalized.
"""
from dataclasses import dataclass
from bs4 import BeautifulSoup
from io import StringIO
import pandas as pd
import numpy as np
from scipy import stats
import glob
import pathlib
import os
import sys


@dataclass
class NanoStringSample:
    sample_path: str

    def __post_init__(self):
        with open(self.sample_path) as file_in:
            soup = BeautifulSoup(file_in.read(), "html.parser")
        tags = [tag.name for tag in soup.find_all()]
        for tag in tags:
            tag_string = getattr(soup, tag).string.strip()
            tag_string = StringIO(tag_string)
            if tag == "code_summary":
                df = pd.read_csv(tag_string)
                df = df[df["Name"].notnull()]
                df = df.set_index(["CodeClass", "Name", "Accession"])
                df = df.astype("float")
            else:
                df = pd.read_csv(tag_string, names=["Attribute", "Value"])
                df = df.set_index("Attribute")
            df.columns.name = tag
            setattr(self, tag, df)
        id = self.lane_attributes.loc["ID"][0]
        setattr(self, "ID", id)


@dataclass
class NanoString:
    rcc_dir: str

    def compile_samples(self, attribute):
        for sample in self.samples:
            to_merge = getattr(sample, attribute)
            to_merge.columns = ["Sample " + sample.ID]
            if "df" not in locals():
                df = to_merge
            else:
                df = df.merge(to_merge, left_index=True, right_index=True)

        cols = sorted(df.columns, key=lambda x: int(x.split(" ")[1]))
        df = df.reindex(cols, axis=1)
        return df

    def __post_init__(self):
        samples = []
        for rcc_file in glob.glob(os.path.join(self.rcc_dir, "*.RCC")):
            sample = NanoStringSample(rcc_file)
            samples.append(sample)
        setattr(self, "samples", samples)

        attrs = []
        for sample in self.samples:
            for attr in vars(sample).keys():
                if attr not in attrs:
                    attrs.append(attr)
        attrs.remove("sample_path")
        attrs.remove("ID")
        for attr in attrs:
            df = self.compile_samples(attr)
            setattr(self, attr, df)

    def counts_norm(self, type="positive", takeLog=True, export=False):
        valid = ["Positive", "Housekeeping"]
        if type.title() not in valid:
            raise ValueError(f"counts_norm: type must be one of {valid}.")

        df = self.code_summary
        df_pos = df.loc[df.index.get_level_values("CodeClass") == type.title()]
        geo_mean = stats.gmean(df_pos)
        a_mean = geo_mean.mean()
        norm_factor = [a_mean / geo for geo in geo_mean]

        for i, factor in enumerate(norm_factor):
            df.iloc[:, i] = df.iloc[:, i].apply(lambda x: x * factor)
            if takeLog:
                df.iloc[:, i] = np.where(df.iloc[:, i] > 0, np.log2(df.iloc[:, i]), -1)
        if type.title() == "Positive":
            df.columns.name = "Positive Control Normalization"
        elif type.title() == "Housekeeping":
            df.columns.name = "CodeSet Content Normalization"
        setattr(self, "counts_norm", df)
        return df

    def export(self, attr):
        df = getattr(self, attr)
        dir_name = pathlib.Path(self.rcc_dir).stem
        filename = f"{dir_name}-{attr}.csv"
        df.to_csv(filename)

    def export_all(self):
        dont_export = ["rcc_dir", "samples", "header", "messages", "counts_norm"]
        for attr in vars(self).keys():
            if attr not in dont_export:
                self.export(attr)


if __name__ == "__main__":
    rcc_dir = sys.argv[1]
    nanostring = NanoString(rcc_dir)
    nanostring.export("code_summary")
