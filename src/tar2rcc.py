import tarfile
import re
import glob
import gzip
import pathlib
import os


def make_rcc_dir(tar_file):
    """Extracts all files from a tarfile into a new directory with '_RCC' appended

    Args:
        tar_file (str): must be a tarfile

    Returns:
        str: a path to new "_RCC" directory
    """
    if not tarfile.is_tarfile(tar_file):
        raise ValueError(f"{tar_file} is not a tarfile")
    tar_path, filename = os.path.split(os.path.abspath(tar_file))
    rcc_dir_name = re.search(".*[^.tar]", filename).group(0)
    rcc_dir_name = rcc_dir_name + "_RCC"
    rcc_dir = os.path.join(tar_path, rcc_dir_name)

    try:
        os.mkdir(rcc_dir)
    except:
        pass

    with tarfile.open(tar_file) as tar:
        tar.extractall(rcc_dir)

    return rcc_dir


def gunzip_files(my_dir: str) -> None:
    """Unzips any gzipped files in the given directory

    Args:
        my_dir (str): path to directory for files to be gunzipped
    """
    for gzip_file in glob.glob(os.path.join(my_dir, "*.gz")):
        file_stem = pathlib.Path(gzip_file).stem
        gunzip_file = os.path.join(my_dir, file_stem)
        try:
            if not os.path.exists(gunzip_file):
                with gzip.open(gzip_file, "rb") as file_in, open(
                    gunzip_file, "wb"
                ) as file_out:
                    for line in file_in:
                        file_out.write(line)
            os.remove(gzip_file)
        except Exception as e:
            print(e)
