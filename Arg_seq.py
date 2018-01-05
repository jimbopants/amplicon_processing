#!/Users/jimbo/anaconda/bin/python
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


### Imports:
import pandas as pd
import Arg_seq_functions as ll
import time
from datetime import datetime
import argparse
import glob
import subprocess
import sys
import os

def main():
    # Print Version and Header
    ll.print_name()
    ll.version_info(version, author, date)
    
    # Parse arguments, initialize directory, parameters and ask for user input and store settings
    args = parse_arguments()
    params = ll.init_params(args.params)
    paths = ll.construct_directory_structure(args.base_dir)
    db_dict = ll.load_databases()
    option_check_and_set(args, params)
    logfile = ll.write_settings(paths["logging"], version, params)


### Workflow commands ###
    #Uses subprocess shell=True formatting, runs scripts in current environment. 
#1. Merge paired end reads using join_paired_ends.py.
    if any( [args.merge, args.all_thru_relabel] ):
        ll.run_command( "join_paired_ends.py -f {0} -r {1} -o {2} -b {3}".format(
                args.forward, args.reverse, paths['paired'], args.barcodes
                ), "Merging Paired Ends" )
    
#2. Relabels reads using split_libraries_fastq.py  
    if any( [args.demult, args.all_thru_relabel] ): 
        ll.run_command( "split_libraries_fastq.py -i {0} -o {1} -b {2} -m {3} --barcode_type=12 --store_demultiplexed_fastq -r 0 -q 0 -n 100".format(
        paths['paired']+"/fastqjoin.join.fastq", paths['split_libraries'], paths['paired_barcodes'], args.m
        ), "Demultiplexing" )            
    
#3. Remove primers w/ cutadapt (Not necessary if merged 515-806 reads are already ~250 bp)
    if any( [args.cutadapt, args.all_thru_relabel] ): 
        ll.run_command( "cutadapt -a {0}...{1}$ -o {2} {3}".format(
            params['fwd_primer'], params['rev_primer'], paths['no_adapter'], paths['fastq_seqs']
            ), "Trimming Forward and Reverse Primers" )
    
#4. Vsearch quality filter at specified settings:
    if any( [args.filter, args.all_thru_relabel] ):
        if os.path.isfile(paths['no_adapter']):
            seqs = paths['no_adapter']
        else:
            seqs = paths['fastq_seqs']

        ll.run_command( "vsearch --fastq_filter {0} --fastq_maxee {1} --fastq_trunclen {2} --fastaout {3}".format(
            seqs, params['max_ee'], params['trunc_len'], paths['filtered']
            ), "Vsearch Quality Filtering" )

#5. Relabel for vsearch
    if any( [args.relabel, args.all_thru_relabel] ): 
        ll.arg_seq_relabeler(paths['filtered'], paths['relabeled_fasta'])


#6a. Dereplicate, abundance sort and save a version with singletons discarded
    if any( [args.cluster, args.table_workflow]):
        ll.run_command("vsearch --derep_fulllength {0} --output {1} --sizeout --minuniquesize 2".format(
            paths["relabeled_fasta"], paths["derep_mc2"]
            ), "Dereplicating Seqs"
        )    

#6b. OTU clustering
        ll.run_command("vsearch --cluster_fast {0} --id {1} --centroids {2} --relabel OTU_ --sizeout --sizein --uc {3}".format(
            paths["derep_mc2"], params['id_thresh'], paths["otu_centroids"], paths["otu_clusters"]
            ), "Clustering Reads"
        )
#6c. Chimera filtering (de novo and reference database)
        ll.run_command("vsearch --uchime_denovo {0} --nonchimeras {1} --chimeras {2} --log {3}".format(
            paths["otu_centroids"], paths["non_chimeras_dn"], paths["chimeras_dn"], paths["dn_log"]
            ), "Filtering Chimeras"
        )
#6d. Map reads (including singletons) back to OTUs.
        ll.run_command("vsearch --usearch_global {0} -db {1} -id {2} -uc {3}".format(
            paths["relabeled_fasta"], paths["non_chimeras_dn"], params['id_thresh'], paths["otu_map"]
            ), "Mapping reads"
        )

#7. Assign taxonomy to OTUs 
    if any( [args.tax, args.table_workflow] ):
        
        # Update if I add any other options to this.
        if args.tax == 'rdp':
            tax_db = db_dict["rdp_database"]
        elif args.tax == 'silva':
            tax_db = db_dict["silva99_tax"]
            ref_db = db_dict["silva99_ref"]
        else:
            tax_db = db_dict["gg_tax_db"]


        ll.run_command( "assign_taxonomy.py -i {0} -t {1} -o {2} -r {3}".format(
            paths["non_chimeras_dn"], tax_db, paths["tax"], ref_db
            ), "Assigning taxonomy")

        ll.add_header(paths["tax"]+"/non_chimeras_dn_tax_assignments.txt")



#8. Create OTU table
    if any( [args.table, args.table_workflow] ):
        make_table = "python /drive5/uc2otutab.py {0}".format( paths["otu_map"] )
        with open(paths["otu_table"], "w") as g:
            for line in ll.run_command_and_pipe(make_table):
                if paths["otu_map"] not in line: #Don't write header line.
                    g.write(line)

        ll.run_command("biom convert -i {0} -o {1} --table-type='OTU table' --to-json".format(
            paths["otu_table"], paths["biom_table"]
            ), "Generating txt biom")

        ll.run_command('biom add-metadata -i {0} -o {1} --observation-metadata-fp {2} --sc-separated taxonomy'.format(
            paths["biom_table"], paths["biom_wtax"], paths["tax"]+"/non_chimeras_dn_tax_assignments_fixed.txt"
            ), "Adding taxonomy")

        ll.run_command('summarize_taxa.py -i {0} -o {1}'.format(
            paths["biom_wtax"], paths["gg_tax_summary"]
            ), "Summarizing Taxa")

        ll.run_command("biom summarize-table -i {0} -o {1}".format(
            paths["biom_wtax"], paths["vsearch"]+"biom_summary.txt"
            ), "Generating Biom Summary")

#9. Alignment and Tree Construction
    if any([args.align_tree, args.table_workflow] ):
        ll.run_command("align_seqs.py -i {0} -t {1} -o {2}".format(
            paths["non_chimeras_dn"], db_dict["gg_aligned"], paths["aligned"]
            ), "Aligning Seqs")

        ll.run_command( "make_phylogeny.py -i {0} -o {1} -l {2}".format(
        paths["aligned"]+"non_chimeras_dn_aligned.fasta", paths["vsearch"]+"phylo_tree.tre", paths["vsearch"]+"tree_log.txt"
        ), "Make phylogenetic tree")


# Functions:
def parse_arguments():
    parser = argparse.ArgumentParser()

## Step 1: Merge, Add barcodes to seq headers, quality filter, and relabel reads for vsearch
    parser.add_argument("--merge", help="Pair Sequences with Qiime/join_paired_ends.py", action='store_true')
    parser.add_argument("--demult", help="Demultiplex sequences with Qiime/split_libraries_fastq.py", action='store_true')
    parser.add_argument("--cutadapt", help="Demultiplex sequences with Qiime/split_libraries_fastq.py", action='store_true')
    parser.add_argument("--filter", help="Quality Filter with vsearch/fastq_filter",action='store_true')
    parser.add_argument("--relabel", help="Relabel sequences to be compatible with vsearch", action='store_true')
    parser.add_argument("--all_thru_relabel", help="Merge, Add barcodes to seq headers, quality filter, and relabel for vsearch", action='store_true')

## Step 2: Construct biom table and align OTU tables:
    parser.add_argument("--cluster", help="Derep reads, abundance sort, cluster OTUs, chimera filter and map reads",action='store_true')
    parser.add_argument("--table", help="Make OTU table and assign taxonomy",action='store_true')
    parser.add_argument("--align_tree", help="Align seqs",action='store_true')
    parser.add_argument("--table_workflow", help="Combines cluster_filter, make_table and align_tree", action='store_true')
    parser.add_argument("--tax", help="Assign taxonomy", action="store",
            choices=['gg','rdp', 'silva'])

# Files to pass
    parser.add_argument('--base_dir', help='base_directory')
    parser.add_argument("--m", help="Map file")
    parser.add_argument("--forward", help="Forward Reads in Argonne fmt")
    parser.add_argument("--reverse", help="Reverse Reads in Argonne fmt")
    parser.add_argument("--barcodes", help="Barcode file")
    parser.add_argument("--params", help="tab sep text file containing parameters in fmt: nametabvalue")
    parser.add_argument("--QC", help="Logs sequences passing QC at each step, can add significant overhead")

# Parse Command Line Arguments:
    try:
        result = parser.parse_args()
        return result
    except Exception as e: 
        parser.print_help()
        print(e)
        sys.exit(0)
    

def option_check_and_set(args, params):
    """
    After parsing arguments and reading the parameter file, 
    # this function prints relevant options and checks for any missing arguments, then asks for user input.
    """
# Error out if initial sequence files are missing and merging is seleceted. 
    if (args.merge or args.all_thru_relabel) and not(args.forward and args.reverse and args.barcodes):
        raise argparse.ArgumentTypeError('Forward, Reverse and Barcodes files required when performing merging!')

# Set OTU clustering criteria:
    if (args.filter or args.table_workflow):
        if params['id_thresh'] is None:
            params['id_thresh'] = input('[Cluster] Enter OTU clustering threshold (default = 0.97): ')
        
        if params['max_ee'] is None:
            params['max_ee'] = input('[Filter] Max Errors (~.5-1 per 100bp or 1-3 per 250bp): ')
        
        if params['trunc_len'] is None:
            params['trunc_len']= input('[Cutadapt] Enter length to trim sequences to: ')
    
    if args.cluster:
        if params['id_thresh'] is None:
            params['id_thresh'] = input('[Cluster] Enter OTU clustering threshold (default = 0.97): ')
 


# Set cutadapt options:
    if (args.cutadapt):
        if params['fwd_primer'] is None:
            params['fwd_primer'] = input('[Cutadapt] Enter forward primer: ')
        
        if params['rev_primer'] is None:
            params['rev_primer'] = input('[Cutadapt] Enter reverse primer: ')
        
    return params

if __name__ == "__main__":
    main()


