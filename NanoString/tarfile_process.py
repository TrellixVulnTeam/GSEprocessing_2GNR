import tarfile
import re
import glob
import gzip
import os
import sys


def make_tarfile(source_dir, output_filename):
    try:
        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))
        print(f"{source_dir} converted to .tar.gz file")
    except Exception as e:
        print(e)


def make_rcc_dir(tar_file):
    tar_path, filename = os.path.split(os.path.abspath(tar_file))
    rcc_dir_name = re.search(".*[^.tar]", filename).group(0)
    rcc_dir_name = rcc_dir_name + "_RCC_files"
    rcc_dir = os.path.join(tar_path, rcc_dir_name)
    try:
        os.mkdir(rcc_dir)
    except:
        pass

    with tarfile.open(tar_file) as tar:
        tar.extractall(rcc_dir)
    return rcc_dir


def gunzip_files(my_dir):
    for gzip_file in glob.glob(os.path.join(my_dir, "*.gz")):

        base = os.path.basename(gzip_file)
        gunzip_file = os.path.join(my_dir, base[:-3])
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


if __name__ == "__main__":
    file = sys.argv[1]
    rcc_dir = make_rcc_dir(file)
    gunzip_files(rcc_dir)
