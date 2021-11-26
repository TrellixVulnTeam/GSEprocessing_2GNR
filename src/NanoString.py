from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy import stats
import glob
import pathlib
import os

from .NanoStringSample import NanoStringSample


@dataclass
class NanoString:
    input_dir: str

    def __post_init__(self):
        samples = []
        for rcc_file in glob.glob(os.path.join(self.input_dir, "*.RCC")):
            sample = NanoStringSample(rcc_file)
            samples.append(sample)
        self.samples = samples

        # create list of all attributes in any sample
        attrs = ["file_in", "ID"]
        for sample in self.samples:
            for attr in vars(sample).keys():
                if attr not in attrs:
                    attrs.append(attr)

        # remove attrs that wont' be compiled
        for attr in attrs[2:]:
            df = self.compile_samples(attr)
            setattr(self, attr, df)

    def compile_samples(self, attribute: str) -> pd.DataFrame:
        """Compile each samples attribute into one dataframe including all samples

        Args:
            attribute (str): NanoStringSample attribute to be compiled

        Returns:
            pd.DataFrame: DataFrame with all NanoStringSample info
        """
        for sample in self.samples:
            # try statement allows attributes in only some samples
            try:
                to_merge = getattr(sample, attribute)
                to_merge.columns = [sample.ID]
                if "df" not in locals():
                    df = to_merge
                else:
                    df = df.merge(to_merge, left_index=True, right_index=True)
            except:
                pass

        cols = sorted(df.columns, key=lambda x: int(x.split(" ")[1]))
        df = df.reindex(cols, axis=1)
        return df

    def counts_norm(
        self,
        takeLog: bool = True,
        type: str = "positive",
    ) -> pd.DataFrame:
        """Normalize the NanoString's raw counts

        Args:
            takeLog (bool, optional): take log of every count value prior to Norm. Defaults to True.
            type (str, optional): type of normalization to be used. Defaults to "positive".

        Raises:
            ValueError: if normalization type called doesn't exist

        Returns:
            pd.DataFrame: Normalized dataframe same as input shape
        """

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

    def export(self, attr: str) -> None:
        """Export to given output_dir

        Args:
            attr (str): attribute to export
        """
        df = getattr(self, attr)
        dir_name = pathlib.Path(self.input_dir).stem
        filename = f"{dir_name}-{attr}.csv"
        df.to_csv(filename)

    def export_all(self) -> None:
        """Export all attributes"""
        dont_export = [
            "input_dir",
            "output_dir",
            "samples",
            "header",
            "messages",
            "counts_norm",
        ]
        for attr in vars(self).keys():
            if attr not in dont_export:
                self.export(attr)
