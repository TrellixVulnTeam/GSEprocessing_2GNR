"""youtube video: https://www.youtube.com/watch?v=AXAPR44zitQ
data analysis guidelines: https://www.nanostring.com/wp-content/uploads/2020/12/Gene_Expression_Data_Analysis_Guidelines.pdf

    Returns:
        [type]: [description]
    """

from dataclasses import dataclass
import GEOparse
import tarfile
import gzip
import os
import re
import pandas as pd
from scipy import stats
import sys

# import shutil


def get_gse(accessionID):
    path = os.path.join(".", accessionID)

    gse = GEOparse.get_GEO(geo=str(accessionID), destdir=path)

    # soft_filename = str(accessionID) + "_family.soft.gz"
    # with gzip.open(soft_filename, 'rb') as f_in and open(soft_file_name[:3] + ".txt", 'wb') as f_out:
    #     shutil.copyfileobj(f_in, f_out)

    return gse


@dataclass
class NanoStringSample:
    sample_path: str

    def __post_init__(self):
        with gzip.open(self.sample_path, "rt") as file_in:
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
    accessionID: str

    def __post_init__(self):
        tar_path = self.accessionID + "_RAW.tar"
        with tarfile.open(tar_path) as tar:
            try:
                rcc_dir = self.accessionID + "_RCC_Files"
                os.mkdir(rcc_dir)
                tar.extractall(rcc_dir)
            except:
                pass
            samples = []
            for sample_filename in os.listdir(rcc_dir):
                sample_path = os.path.join(rcc_dir, sample_filename)
                sample = NanoStringSample(sample_path)
                samples.append(sample)
            setattr(self, "samples", samples)

    def raw_counts(self):
        df = self.samples[0].code_summary
        df.columns = ["Sample 1"]
        for i, sample in enumerate(self.samples[1:]):
            to_merge = sample.code_summary
            to_merge.columns = ["Sample " + str(i + 2)]
            df = df.merge(to_merge, left_index=True, right_index=True)
        return df

    def counts_norm(self, type):
        valid = ["Positive", "Housekeeping"]
        if type.title() not in valid:
            raise ValueError(f"counts_norm: type must be one of {valid}.")

        df = self.raw_counts()
        df_pos = df.reset_index()
        df_pos = df_pos[df_pos.CodeClass == type.title()]
        df_pos = df_pos.set_index(["CodeClass", "Name", "Accession"])
        geo_mean = stats.gmean(df_pos)
        amean = geo_mean.mean()
        norm_factor = [amean / geo for geo in geo_mean]
        for i, factor in enumerate(norm_factor):
            df.iloc[:, i] = df.iloc[:, i].apply(lambda x: x * factor)
        if type.title() == "Positive":
            df.columns.name = "Positive Control Normalization"
        elif type.title() == "Housekeeping":
            df.columns.name = "CodeSet Content Normalization"
        return df


if __name__ == "__main__":
    gse = sys.argv[1]
    nanostring = NanoString(gse)
    print(nanostring.counts_norm("positive"))
