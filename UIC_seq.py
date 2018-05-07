#!/Users/jimbo/anaconda/bin/python
### Description: 
# This is a wrapper for existing command line tools with a bit of code for managing user input and performing multiple steps at once.
# It is an update to the UIC formatted 16S processing script. 
# There are slight differences in the first 4 steps due to difference in read formatting between Argonne and UIC.
### UIC Format Description:
# Argonne deliver reads in a 3 file format - 1 forward and reverse file and a separate file with barcodes. 
# By default vsearch and QIIME need concatenated single files, so this is performed after read merging and QC.


### Author and version information:
version = "0.9"
author = "Jim Griffin"
date = "1/5/18"
git_repo = 'finish'


### ChangeLog and Version History:
#   1.0 1/5/18 - Updated docs and uploaded to github. Terminal logging not fully implemented ATM.
#   0.9 10/20/17 - Main Functionality complete.
### 0.9 Changelog
# Added parameter file import option
# Merged most custom subprocess calls into a standardized function
# Reorganized upper level modules into main module and functions calls
# Added parameter logging and user input at run-time to avoid 
# Removed diversity analyses pipelines since these are not one size fits all and I want more control over these.


### Docs:

###Environment/other requirements:
# Macqiime 1.9
# conda python 3.5
# vsearch 1.11.1
# cutadapt 1.14

### Function descriptions and setup:
# arg_seq_functions.py must be in the same directory as this script or the path modified.
# reads should be located in base_dir/raw in order of: (Sample1_R1,Sample1_R2,Sample2_R1...). Can be zipped.
# parse_aguments: Reads arguments and saves them to a dictionary called args.
# init_params: Reads a tab separated text file of (parameter value) lines and returns a dictionary called params.
# db_dict: Set this for the paths on your computer in the Arg_seq_functions 
# Taxonomy databases supported: Greengenes, Silva, or RDP currently (Requires downloading an external database)
# Silva db notes (Super Useful): https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_notes.txt


### Seq processing workflow using subprocess to call individual programs:
# 1. Pair Ends with PEAR
# 2. Remove primers with cutadapt ((Primers should be listed as 5'...3'(RC))
# 3. QC using vsearch
# 4. Relabel sequences for vsearch using custom script. ##! Add label requirement options. 
# 5. Vsearch cluster reads
# 6. Assign taxonomy to OTUs using a specififed database
# 7. Create Qiime1 style BIOM table, add taxonomy, and summarize
# 8. Align seqs and construct a tree using Qiime1's default options


### TODO:
# Enable logging of all terminal output and errors
# Fix assign taxonomy with GG and RDP. Add support for the vsearch
# consensus classifier.
# Migrate code to Qiime2 if they fix vectorization issues with DADA2 performance and make better tools for interacting with artifacts.
# Write tests. LOL probably not...


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
    ll.print_name_uic()
    ll.version_info(version, author, date)
    
    # Parse arguments, initialize directory, parameters and ask for user input and store settings
    args = parse_arguments()
    params = ll.init_params(args.params)
    paths = ll.construct_directory_structure(args.base_dir)
    db_dict = ll.load_databases()
    option_check_and_set(args, params)
    logfile = ll.write_settings(paths["logging"], version, params)


### Workflow commands ###
# 1. Pair with PEAR
    if any([args.merge, args.all_thru_relabel]):

        print ("Paring reads found in {0}".format(args.base_dir+'raw/' ) )
        ll.pear_files(args.base_dir+'raw/', args.base_dir)
        print ("\n\n\nReads paired")

# 2. Remove primers w/ cutadapt (Not necessary if merged 515-806 reads are already ~250 bp)
    if any( [args.cutadapt, args.all_thru_relabel] ):
        for file_in in glob.glob(paths['assembled']+'*'):
            name = file_in.split('/')[-1].split('.')[0]
            ll.run_command( "cutadapt -a {0}...{1}$ -o {2} --discard-untrimmed {3}".format(
                params['fwd_primer'], params['rev_primer'],
                paths['split_libraries']+name+'.no_adapters.fq', file_in
            ), "Trimming Forward and Reverse Primers" )

# 3. Vsearch quality filter at specified settings:
    if any( [args.filter, args.all_thru_relabel] ):
        ll.ensure_path(paths['uic_filtered'])
        for file_in in glob.glob(paths['split_libraries']+'*'):
            print(file_in)
            name = file_in.split('/')[-1].split('.')[0]
            ll.run_command( "vsearch --fastq_filter {0} --fastq_maxee_rate {1} --minseqlength {2} --maxseqlength {3} --fastaout {4}".format(
            file_in, params['max_ee_rate'], params['min_len'], params['max_len'], paths['uic_filtered']+name+'_filtered.fa'
            ), "Vsearch Quality Filtering" )
            
#4. Relabeling sequences in fasta files. This is unique to UIC data
    if any( [args.relabel, args.all_thru_relabel] ):
        ll.uic_seq_relabeler(paths)


#5a. Dereplicate, abundance sort and save a version with singletons discarded
    if any( [args.cluster, args.table_workflow]):
        ll.run_command("vsearch --derep_fulllength {0} --output {1} --sizeout --minuniquesize 5".format(
            paths["relabeled_fasta"], paths["derep_mc2"]
            ), "Dereplicating Seqs"
        )    

#5b. OTU clustering
        ll.run_command("vsearch --cluster_fast {0} --id {1} --centroids {2} --relabel OTU_ --sizeout --sizein --uc {3}".format(
            paths["derep_mc2"], params['id_thresh'], paths["otu_centroids"], paths["otu_clusters"]
            ), "Clustering Reads"
        )
#5c. Chimera filtering (de novo and reference database)
        ll.run_command("vsearch --uchime_denovo {0} --nonchimeras {1} --chimeras {2} --log {3}".format(
            paths["otu_centroids"], paths["non_chimeras_dn"], paths["chimeras_dn"], paths["dn_log"]
            ), "Filtering Chimeras"
        )
#5d. Map reads (including singletons) back to OTUs.
        ll.run_command("vsearch --usearch_global {0} -db {1} -id {2} -uc {3}".format(
            paths["relabeled_fasta"], paths["non_chimeras_dn"], params['id_thresh'], paths["otu_map"]
            ), "Mapping reads"
        )

#6. Assign taxonomy to OTUs 
    if any( [args.tax, args.table_workflow] ):
        
        # Update if I add any other options to this.
        if args.tax == 'rdp':
            tax_db = db_dict["rdp_database"]
        elif args.tax == 'gg':
            tax_db = db_dict["gg_tax_db"]
        elif args.tax == 'euk_silva':
            tax_db = db_dict['silva_euk_tax']
            ref_db = db_dict['silva_euk_ref']
        else:
            #(Assumes args.tax == silva)
            tax_db = db_dict["silva99_tax"]
            ref_db = db_dict["silva99_ref"]


        ll.run_command( "assign_taxonomy.py -i {0} -t {1} -o {2} -r {3}".format(
            paths["non_chimeras_dn"], tax_db, paths["tax"], ref_db
            ), "Assigning taxonomy")

        ll.add_header(paths["tax"]+"/non_chimeras_dn_tax_assignments.txt")



#7. Create OTU table
    if any( [args.table, args.table_workflow] ):
        make_table = "python /drive5/uc2otutab.py {0}".format( paths["otu_map"] )
        with open(paths["otu_table"], "w") as g:
            for line in ll.run_command_and_pipe(make_table):
                line_str = line.decode('utf-8')
                if paths["otu_map"] not in line_str: #Don't write header line.
                    g.write(line_str)

        ll.run_command("biom convert -i {0} -o {1} --table-type='OTU table' --to-json".format(
            paths["otu_table"], paths["biom_table"]
            ), "Generating txt biom")

        # Check if taxonomic data exists and add to the biom table: (Skips for functional gene amplicon data)
        if os.path.isfile(paths["tax"]+"/non_chimeras_dn_tax_assignments.txt"):
            ll.run_command('biom add-metadata -i {0} -o {1} --observation-metadata-fp {2} --sc-separated taxonomy'.format(
            paths["biom_table"], paths["biom_wtax"], paths["tax"]+"/non_chimeras_dn_tax_assignments_fixed.txt"
            ), "Adding taxonomy")

            ll.run_command('summarize_taxa.py -i {0} -o {1}'.format(
            paths["biom_wtax"], paths["gg_tax_summary"]
            ), "Summarizing Taxa")

            ll.run_command("biom summarize-table -i {0} -o {1}".format(
            paths["biom_wtax"], paths["vsearch"]+"biom_summary.txt"
            ), "Generating Biom Summary")

            ll.run_command("biom convert -i {0} -o {1} --to-tsv --header-key taxonomy".format(
            paths["biom_wtax"], paths["txt_biom_w_tax"]
            ), "Generating txt biom")
        


#8. Alignment and Tree Construction
    if any([args.align_tree, args.table_workflow] ):
        ll.run_command("align_seqs.py -i {0} -m muscle -o {2}".format(
            paths["non_chimeras_dn"], db_dict["gg_aligned"], paths["aligned"]
            ), "Aligning Seqs")

        ll.run_command( "make_phylogeny.py -i {0} -o {1} -l {2}".format(
        paths["aligned"]+"non_chimeras_dn_aligned.fasta", paths["vsearch"]+"phylo_tree.tre", paths["vsearch"]+"tree_log.txt"
        ), "Make phylogenetic tree")


# Functions:
def parse_arguments():
    parser = argparse.ArgumentParser()

## Part 1: Merge PE reads, Add barcodes to seq headers, quality filter, and relabel & concatenates reads into one file for vsearch
    parser.add_argument("--merge", help="1. Pair Sequences with PEAR", action='store_true')
    parser.add_argument("--cutadapt", help="2. Trim primers. Uses parameter file values or run time specified options if not given in parameter file.", action='store_true')
    parser.add_argument("--filter", help="3. Quality Filter with vsearch/fastq_filter",action='store_true')
    parser.add_argument("--relabel", help="4. Relabel sequences to be compatible with vsearch & concatenates demultiplexed samples", action='store_true')
    parser.add_argument("--all_thru_relabel", help="Merge PE reads, trim primers, quality filter, relabel for vsearch", action='store_true')

## Step 2: Construct biom table and align OTU tables:
    parser.add_argument("--cluster", help="5. [vsearch] Derep reads, abundance sort, cluster OTUs, chimera filter and map reads",action='store_true')
    parser.add_argument("--tax", help="6. Assign taxonomy", action="store",
            choices=['gg','rdp', 'silva', 'euk_silva'])
    parser.add_argument("--table", help="7. [qiime] Make OTU table and adds taxonomy, if it exists",action='store_true')
    parser.add_argument("--align_tree", help="8. Align seqs",action='store_true')
    parser.add_argument("--table_workflow", help="Combines cluster, tax, table and align_tree", action='store_true')

# Files to pass
    parser.add_argument('--base_dir', help='base_directory')
    parser.add_argument("--m", help="Map file")
    parser.add_argument("--params", help="tab sep text file containing parameters in fmt: nametabvalue")
    parser.add_argument("--QC", help="Logs sequences passing QC at each step, can add significant overhead")


# Print help if no options given. 
    if len(sys.argv)==1:
        parser.print_help()
        print("\n\n!!!!! Read the f'n manual !!!!!\n\n")
        sys.exit(1)

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

# Set OTU clustering criteria:
    if (args.filter or args.table_workflow):
        if params['id_thresh'] is None:
            params['id_thresh'] = input('[Cluster] Enter OTU clustering threshold (default = 0.97): ')
        
        if params['max_ee_rate'] is None:
            params['max_ee_rate'] = input('[Filter] Max Errors (~.5-1 per 100bp or 1-3 per 250bp): ')
        
        if params['min_len'] is None:
            params['min_len']= input('[Cutadapt] Enter minimum sequence length: ')
    
        if params['max_len'] is None:
            params['max_len']= input('[Cutadapt] Enter maximum sequence length: ')


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


