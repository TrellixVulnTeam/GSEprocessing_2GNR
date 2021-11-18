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
class NanoStringSample:
    """Class to hold tar metadata NanoString files"""

    sample_path: str

    def __post_init__(self):
        """Reads each file from .tar, splits into lists per <tag> , converts to dataframe
        and adds as attribute using <tag> name.
        """
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
                    df.columns.name = tag
                    setattr(self, attr, df)
                except Exception as e:
                    pass


def read_nanostring_tar(filename):
    """Reads .tar file and creates directory TarFile where all tar files are extracted
    Then iterates through all extracted files and creates NanoString objects

    Args:
        filename (str): .tar file for conversion to NanoString

    Returns:
        list: a list of NanoString objects. One NanoString is created per file
    """
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


def raw_counts(samples):
    """Iterates through all NanoString objects and pulls the Code_Summary attribute
    dataframe which includes raw read counts. These are then all merged together to create
    one dataframe with all sample reads.

    Args:
        samples (list): list of NanoString objects

    Returns:
        pd.DataFrame: Singular dataframe containing Code_Summary read counts for all files
    """
    df = samples[0].code_summary
    df.columns = ["Sample 1"]
    for i, sample in enumerate(samples[1:]):
        to_merge = sample.code_summary
        to_merge.columns = ["Sample " + str(i + 2)]
        df = df.merge(to_merge, left_index=True, right_index=True)
    return df


def raw_counts_norm(df):
    """Normalizes input df via (1) calculating sample geometric mean (gmean), (2) calculating the
    arithmetic mean of all gmeans and (3) Dividing amean by each gmean to generate a normalization
    factor for each sample. This factor is then multipled across the entire sample column. Process
    pulled from NanoString 'Gene Expression Data Analysis Guideline' pg. 13.

    Args:
        df (pd.DataFrame): A dataframe containing raw read counts for all samples

    Returns:
        pd.DataFrame: Normalized dataframe with the respective normalization factor
        multiplied to every cell.
    """
    df_pos = df.reset_index()
    df_pos = df_pos[df_pos.CodeClass == "Positive"]
    df_pos = df_pos.set_index(["CodeClass", "Name", "Accession"])
    geo_mean = stats.gmean(df_pos)
    amean = geo_mean.mean()
    norm_factor = [amean / geo for geo in geo_mean]
    for i, factor in enumerate(norm_factor):
        df.iloc[:, i] = df.iloc[:, i].apply(lambda x: x * factor)
    return df


if __name__ == "__main__":
    filename = sys.argv[1] + "_RAW.tar"
    samples = read_nanostring_tar(filename)
    code_summary_df = raw_counts(samples)
    df_norm = raw_counts_norm(code_summary_df)
