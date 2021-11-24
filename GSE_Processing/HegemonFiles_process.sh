#!/bin/bash
args=("$@")
GSE=${args[0]}

mkdir -p $GSE

cd $GSE

cp -Rup ../HegemonFiles.py .
cp -Rup ../HomoSapiens_ENST,ProbeID,Name.txt .

python3 HegemonFiles.py $GSE

rm HomoSapiens_ENST,ProbeID,Name.txt
rm HegemonFiles.py 

cd ..
