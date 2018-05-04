
# This is an update to the Arg!Seq 16S processing script for processing Argonne's EMP formatted reads.
### Author and version information:
version = "0.9"
author = "JGrif"
date = "10/18/17"
# git_repo: finish

### EMP Format Description:
# Argonne deliver reads in a 3 file format - 1 forward and reverse file and a separate file with barcodes.

### Environment:
# Macqiime 1.9
# conda python 3.5
# vsearch 1.11.1
# cutadapt 1.14
# Arg_seq_functions.py in python path
# Taxonomy database: Greengenes, Silva, or RDP currently (Requires downloading an external database)
# Silva db notes (Super Useful): https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_notes.txt
# For now, I am trying the

### Seq processing workflow using subprocess:
# 1. Pair Ends with qiime join_paired_ends
# 2. Demultiplex wih qiime split_libraries
# 3. Remove primers with cutadapt
# 4. QC using vsearch
# 5. Relabel sequences for vsearch using custom script
# 6. Vsearch cluster reads
# 7. Assign taxonomy to OTUs using a specififed database
# 8. Create Qiime1 style BIOM table, add taxonomy, and summarize
# 9. Align seqs and construct a tree using Qiime1's default options

### TODO:
# Enable logging of all terminal output and errors
# Fix assign taxonomy with GG and RDP. Add support for the vsearch
# consensus classifier.
# Migrate code to Qiime2 if they fix vectorization issues with DADA2 performance and make better tools for interacting with artifacts.
# Write tests. LOL probably not...

### 0.9 Changelog
# Added parameter file import option
# Merged most custom subprocess calls into a standardized function
# Reorganized upper level modules into main module and functions calls
# Added parameter logging and user input at run-time to avoid
# Removed diversity analyses pipelines since these are not one size fits all and I want more control over these.
