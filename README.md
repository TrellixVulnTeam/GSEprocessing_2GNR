# GSEprocessing

This repo contains proccessing scripts to parse either NCBI GEO RNA-seq data, NCBI GEO MicroArray data, NCBI GEO NanoString data, or raw NanoString data into four working files for the Boolean Lab Hegemon. The working files include (1) expression file containing all genes and normalized expressoin data (2) index file containing binary location for the expression file (3) survival file containing metadata for the expression file and (4) indexHeader file linking NCBI GSM name to sample names.

Below is a list of all files in the repo

(1) HegemonFiles.py takes a GEO AccessionID and parses the .soft.gz file into relevant expression and metadata. This is converted into the four hegemon files

(2) HegemonFiles_process bash script to run all files in Boolean Lab ssh

(3) HomoSapiens_ENST,ProbeID,Name.txt reference file to map columns for expression file
(4) NanoString_RCCparse.py parses a directory of raw NanoString.RCC files to create four hegemon files (pending)

(5) RNAseq-expr.py parses a supplementary RNA-seq counts file to create a normalized expression file using scipy. Must be customized to fit input count file

(6) Tarfile2RCC.py extracts all files from .tar file and gunzips any .gzip files
