from dataclasses import dataclass
from bs4 import BeautifulSoup
from io import StringIO
import pandas as pd
import numpy as np
from scipy import stats
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


@dataclass
class NanoString:
    rcc_dir: str

    def __post_init__(self):
        samples = []
        for rcc_filename in os.listdir(self.rcc_dir):
            rcc_path = os.path.join(self.rcc_dir, rcc_filename)
            sample = NanoStringSample(rcc_path)
            samples.append(sample)
        setattr(self, "samples", samples)

    def compile_samples(self, attribute):
        df = getattr(self.samples[0], attribute)
        df.columns = ["Sample 1"]
        for i, sample in enumerate(self.samples[1:]):
            to_merge = getattr(sample, attribute)
            to_merge.columns = ["Sample " + str(i + 2)]
            df = df.merge(to_merge, left_index=True, right_index=True)
        return df

    def sample_attributes(self, export=False):
        df = self.compile_samples("sample_attributes")
        if export:
            filename = f"{self.rcc_dir}-sample_attributes.csv"
            df.to_csv(filename)
        return df

    def lane_attributes(self, export=False):
        df = self.compile_samples("lane_attributes")
        if export:
            filename = f"{self.rcc_dir}-lane_attributes.csv"
            df.to_csv(filename)
        return df

    def raw_counts(self, export=False):
        df = self.compile_samples("code_summary")
        if export:
            filename = f"{self.rcc_dir}-raw_counts.csv"
            df.to_csv(filename)
        return df

    def counts_norm(self, type="positive", takeLog=True, export=False):
        valid = ["Positive", "Housekeeping"]
        if type.title() not in valid:
            raise ValueError(f"counts_norm: type must be one of {valid}.")

        df = self.raw_counts()
        df_pos = df.loc[df.index.get_level_values("CodeClass") == type.title()]
        geo_mean = stats.gmean(df_pos)
        amean = geo_mean.mean()
        norm_factor = [amean / geo for geo in geo_mean]

        for i, factor in enumerate(norm_factor):
            df.iloc[:, i] = df.iloc[:, i].apply(lambda x: x * factor)
            if takeLog:
                df.iloc[:, i] = np.where(df.iloc[:, i] > 0, np.log2(df.iloc[:, i]), -1)
        if type.title() == "Positive":
            df.columns.name = "Positive Control Normalization"
        elif type.title() == "Housekeeping":
            df.columns.name = "CodeSet Content Normalization"

        if export:
            filename = f"{self.rcc_dir}-counts_norm-{type.title()}.csv"
            df.to_csv(filename)
        return df


if __name__ == "__main__":
    rcc_dir = sys.argv[1]
    nanostring = NanoString(rcc_dir)
    print(nanostring.counts_norm())
    # nanostring.sample_attributes(export=True)
    # nanostring.lane_attributes(export=True)
    # nanostring.raw_counts(export=True)
    # nanostring.counts_norm(takeLog=True, export=True)
