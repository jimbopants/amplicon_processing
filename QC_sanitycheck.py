#!/Users/jimbo/anaconda/bin/python
# This script checks the number of sequences lost at various analysis points and prints the output to terminal.

### Author and version information:
version = "1.0"
author = "JGrif"
date = "10/20/17"

### Imports
import pandas as pd
import argparse
import sys
import os.path

def main():
    args = parse_arguments()

    # header
    print("\n\n... Calculating Reads Passing QC")
    print("Written by {0} v{1} Last Edited:{2}\n".format(author, version, date) )
    
    files_df = expected_files(args.base_dir)
    check_directory(args.base_dir, files_df)
    files_df, summary_lines = return_seq_counts(files_df)
    print(files_df)
    print("BIOM Summary")
    for i in summary_lines:
        print(i)

# Functions:

def return_seq_counts(files_df):
    # call seq_counter on all but biom summary file:
    for i,seq in files_df.iterrows():

        if i in range(4):
            if seq['file_names'][-1] == 'q':
                identifier = '@'
            if seq['file_names'][-1]  == 'a':
                identifier = '>'
            files_df.ix[i,'seq_number'] = seq_counter(seq['file_names'], identifier)
        
        # Handle biom summary separate
        if i == 4:
            summary_lines = []
            with open(seq['file_names'], 'r') as f:
                for i in range(3):
                    summary_lines.append(f.next().strip())  
    
    files_df['fraction'] = files_df['seq_number'] / files_df.ix[0,'seq_number']
    return files_df, summary_lines



def seq_counter(file_in, identifier):
    # This function counts the header lines in a fasta or fq file. 
    # Because the wrapping and sequence length can vary, this does not use wc -l although it's much faster.
    
    seq_count = 0
    with open(file_in, 'r') as f:
        for line in f.readlines():
            if line.startswith(identifier):
                seq_count += 1
    return(seq_count)


def parse_arguments():

    # Parse Command Line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', help="Base Directory for sequence analysis workflow script")
    parser.add_argument('--log', help = "Logfile to output data to. If blank, output to terminal only.")
    
    try:
        result = parser.parse_args()
        return result
    except Exception as e: 
        parser.print_help()
        print(e)
        sys.exit(0)

def check_directory(base_dir, file_df):
    # Checks for the existence of all files and prints an error message if not found.
    for i,seq in file_df.iterrows():
        if not os.path.isfile(seq['file_names']):
            print ("{0} missing! Expected path: {1}".format(seq['descriptions'], seq['file_names']) )
            sys.exit(0)
    print ("All files present, starting!")
    return



def expected_files(base_dir):
    seq_df = pd.DataFrame(columns = ['descriptions', 'file_names', 'seq_number', 'fraction'])
    descriptions = ["Input sequences" , "Merged reads",  "Reads matching barcodes", "Reads after QC", "Reads mapped to OTUs" ]
    file_names = [base_dir  + "raw/barcodes.fastq", base_dir + "paired/fastqjoin.join.fastq" , base_dir + "split_libraries/seqs.fastq", base_dir + "vsearch/relabeled_fna.fasta", base_dir + "vsearch/biom_summary.txt"]

    seq_df['descriptions'] = descriptions
    seq_df['file_names'] = file_names
    return seq_df

if __name__ == "__main__":
    main()
