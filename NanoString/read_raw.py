from dataclasses import dataclass
import os
import re
import pandas as pd
import numpy as np
from scipy import stats
import sys


@dataclass
class NanoStringSample:
    sample_path: str

    def __post_init__(self):
        with open(self.sample_path) as file_in:
            file_in = file_in.read()
        file_in = re.split("</.*>", file_in)
        file_in = [line.split("\n") for line in file_in]
        for group in file_in:
            try:
                r = re.compile("<.*>")
                tag = list(filter(r.match, group))[0]
                attr = tag[1:-1].lower()
                data = group[group.index(tag) + 1 :]
                data = [d.split(",") for d in data]
                if attr == "code_summary":
                    df = pd.DataFrame(data[1:], columns=data[0])
                    df = df.set_index(["CodeClass", "Name", "Accession"])
                    df = df.astype("float")
                else:
                    df = pd.DataFrame(data)
                    df.columns = ["Attribute", "Value"]
                    df = df.set_index("Attribute")
                df.columns.name = attr
                setattr(self, attr, df)
            except Exception as e:
                pass


@dataclass
class NanoString:
    rcc_dir: str

    def __post_init__(self):
        samples = []
        for sample_filename in os.listdir(self.rcc_dir):
            sample_path = os.path.join(self.rcc_dir, sample_filename)
            sample = NanoStringSample(sample_path)
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

    def sample_attributes(self):
        df = self.compile_samples("sample_attributes")
        return df

    def lane_attributes(self):
        df = self.compile_samples("lane_attributes")
        return df

    def raw_counts(self):
        df = self.compile_samples("code_summary")
        return df

    def counts_norm(self, type="positive", take_log=True):
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
            if take_log:
                df.iloc[:, i] = np.log2(df.iloc[:, i])
        if type.title() == "Positive":
            df.columns.name = "Positive Control Normalization"
        elif type.title() == "Housekeeping":
            df.columns.name = "CodeSet Content Normalization"
        return df


if __name__ == "__main__":
    rcc_dir = sys.argv[1]
    nanostring = NanoString(rcc_dir)
    df = nanostring.counts_norm()
    df.to_csv("test.csv")
