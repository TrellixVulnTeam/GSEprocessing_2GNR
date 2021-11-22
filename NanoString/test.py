import rpy2.robjects.packages as rpackages

base = rpackages.importr("base")
utils = rpackages.importr("utils")
utils.chooseCRANmirror(ind=1)  # select the first mirror in the list

from rpy2.robjects.vectors import StrVector

if not rpackages.isinstalled("ruv"):
    utils.install_packages("ruv")

## Create some simulated data
