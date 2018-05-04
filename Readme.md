
## ArgSeq & UICSeq ReadMe
This repo contains workflow scripts for processing Argonne's EMP formatted reads.

### EMP Format Description:
- Argonne deliver reads in a 3 file format - 1 forward and reverse file and a separate file with barcodes.
- UIC demultiplexes files and sends an individual file for each sample barcode.
Use the appropriate workflow script (Arg/UIC) depending on what format your reads are in.

### Environment:
1. Macqiime 1.9+
2. conda python 3.5
3. vsearch 1.11.1
4. cutadapt 1.14
5. Taxonomy database: Greengenes, Silva, or RDP. Silva & RDP are preferred to greengenes but don't come with QIIME/requires downloading an external database

### Seq processing workflow using subprocess for Argonne reads:
1. Pair Ends with **qiime** *join_paired_ends*
2. Demultiplex wih **qiime** *split_libraries*
3. Remove primers with **cutadapt**
4. QC using **vsearch**
5. Relabel sequences for **vsearch** using custom script
6. **Vsearch** cluster reads
7. Assign taxonomy to OTUs using a specififed database (**qiime**)
8. Create **Qiime** style BIOM table, add taxonomy, and summarize
9. Align seqs and construct a tree using **Qiime**1's default options

### TODO:
1. Enable logging of all terminal output and errors
2. Fix assign taxonomy with GG and RDP. Add support for the vsearch
3. consensus classifier.
4. Migrate code to Qiime2 if they fix vectorization issues with DADA2 performance and make better tools for interacting with artifacts.

### 0.9 Changelog
Added parameter file import option
Merged most custom subprocess calls into a standardized function
Reorganized upper level modules into main module and functions calls
Added parameter logging and user input at run-time
Removed diversity analyses pipelines since these are not one size fits all and I want more control over these.
