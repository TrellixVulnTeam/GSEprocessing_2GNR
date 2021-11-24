# GSEprocessing

This repo contains proccessing scripts to parse either NCBI GEO RNA-seq data, NCBI GEO MicroArray data, NCBI GEO NanoString data, or raw NanoString data into four working .txt files for the Boolean Lab Hegemon. The working files include:
> 1. expression (expr) file containing all genes and normalized expressoin data 
> 2. index file containing binary location for the expression file 
> 3. survival file containing metadata for the expression file
> 4. indexHeader file linking NCBI GSM name to sample names

## Descriptio of all files in the repo

1. **HegemonFiles.py**: Takes a GEO AccessionID as input and creates the four hegemon files above
2. **HegemonFiles_process**: Bash script to run all files in Boolean Lab ssh
3. **HomoSapiens_ENST,ProbeID,Name.txt**: Reference file to map columns for expression file
4. **NanoString_RCCparse.py**: Parse a directory of raw NanoString.RCC files to create four hegemon files (usable, but expr file is pending normalization yet complete)
5. **RNAseq-expr.py**: Parse a supplementary RNA-seq counts file to create a normalized expression file using scipy. Must be customized to fit input count file
6. **Tarfile2RCC.py**: Extract all files from .tar file and gunzips any .gz files
