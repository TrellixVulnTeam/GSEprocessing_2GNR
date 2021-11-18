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
import sys

# import shutil


def get_gse(accessionID):
    """Wanted to make function to pull RAW file from GEO but couldn't figure it out

    Args:
        accessionID (str): GSE ID

    Returns:
        gse: GEOParse.BaseGEO
    """
    path = os.path.join(".", accessionID)

    gse = GEOparse.get_GEO(geo=str(accessionID), destdir=path)

    # soft_filename = str(accessionID) + "_family.soft.gz"
    # with gzip.open(soft_filename, 'rb') as f_in and open(soft_file_name[:3] + ".txt", 'wb') as f_out:
    #     shutil.copyfileobj(f_in, f_out)

    return gse


@dataclass
class NanoString:
    """Class to hold tar metadata NanoString files"""

    path: str
    filename: str

    def __post_init__(self):
        """Reads each file from .tar, splits into lists per <tag> , converts to dataframe
        and adds as attribute using <tag> name.
        """
        with gzip.open(os.path.join(self.path, self.filename), "rt") as file_in:
            file_in = file_in.read()
            file_in = re.split("</.*>", file_in)
            file_in = [line.split("\n") for line in file_in]
            for group in file_in:
                try:
                    r = re.compile("<.*>")
                    tag = list(filter(r.match, group))[0]
                    data = group[group.index(tag) + 1 :]
                    data = [d.split(",") for d in data]
                    df = pd.DataFrame(data)
                    df.index.name = tag
                    setattr(self, tag[1:-1], df)
                except Exception as e:
                    pass
        setattr(self, "metadata_keys", list(vars(self).keys()))


def read_nanostring_tar(filename):
    tar = tarfile.open(filename)
    try:
        os.makedir("TarFiles")
        tar.extractall("./TarFiles")
    except:
        pass
    samples = []
    for sample_filename in os.listdir("./TarFiles"):
        sample = NanoString("./TarFiles", sample_filename)
        samples.append(sample)
    return samples


def final():
    doc = sys.argv[1] + "_RAW.tar"
    samples = read_nanostring_tar(doc)
    df = samples[0].Code_Summary
    for sample in samples[1:]:
        to_merge = sample.Code_Summary
        df = df.merge(to_merge, how="left", on=[0, 1, 2])
    columns = ["Sample " + str(i) for i in range(1, 13)]
    df.columns = ["CodeClass", "Name", "Accession"] + columns
    df = df.drop(0)[:-1]
    df = df.astype({column: "float" for column in columns})
    df = df.set_index(["CodeClass", "Name", "Accession"])
    df["Total"] = df.sum(axis=1)
    print(df.reset_index().groupby("CodeClass").sum()["Total"])


if __name__ == "__main__":
    final()
