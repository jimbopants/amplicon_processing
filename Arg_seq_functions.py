#### Low Level Functions for Sequence Processing ####
#### Jim Griffin ####
#### 10/18/17 ####

# This contains a bunch of functions for 16s read processing automation. 
#################
# 1. Shared Functions
#   a. Parameter Setting (init_params)
#   b. Logging (write_settings
#   c. Subprocess and timer (run_command)
#################
# 2. Core Facility Specific
#   a. UIC - Directory level paired end read merger workflow
#   b. Argonne - Relabel Seqs and concatenate
#################
# 3. Local User Stuff
#   a. 16s/etc. Database file locations
#   b. Directory Setup
#################
# 4. Cosmetic:
#   a. Messaging / header

##### 0. Imports #####
import glob
import subprocess
import sys
import os
import shutil
import pandas as pd
import time
from datetime import datetime

##### 1. Shared Functions #####
def init_params(params):
    """Read a tab separated parameter file and return a dict of param:value
    pairs. Includes None where value is missing.
    Skips lines starting with '#' and will ignore comments separated by tab on
    each line.
    
    """
    
    p_dict = {
            "id_thresh":None,
            "fwd_primer":None,
            "rev_primer":None,
            "min_len":None,
            "max_len":None,
            "max_ee_rate":None,
            "min_OTU_size":None,
            }

    if params is not None:
        with open(params, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line_split = line.rstrip().split('\t')
                p_dict[line_split[0]] = line_split[1]
    return p_dict

def write_settings(log_dir, version, params):
    """Creates a logging file with time, version and parameters, and returns the name of this file
    so that I can add seq counts before and after various processing steps."""
    
    now = datetime.now().strftime('%Y-%m-%d-%H-%M')
    log_file = log_dir + "logfile_{0}.txt".format(now)
    
    with open(log_file, "w") as f:
        f.write('Arg!Seq Log\n')
        f.write("time:\t{0}\n".format(now))
        f.write("version:\t{0}\n\n\n\n".format(version))
                
        for k,v in params.items():
            f.write("{0}\t{1}\n".format(k,v))
    
    return log_file


def run_command(command, description):
    """ Executes commands passed as a single string using the shell """ 
    module_start = time.time()
    
    print ("Starting: {}".format(description))
    
    # TODO:
    # Make sure that check_output is the correct method to call if I want to see error messages 
    # and stop execution of this program. 
    # Check that shell=True works with all options.
    
    subprocess.check_call(command, shell=True)
    

    module_finish = time.time()
    print ("\n\n\n{0} Finished! Time elapsed: {1}sec\n\n\n".format(
        description, str(module_finish - module_start)
            )
        )
    return 

def run_command_and_pipe(command):
    """ Runs a command at command line and returns stdout line by line"""
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, shell=True)
    return iter(p.stdout.readline, b'')

def add_header(taxonomy_file):
    # Edits a taxonomy file to fit the biom format.
    new_file = taxonomy_file.split(".txt")[0]+"_fixed.txt"
    with open (taxonomy_file) as f:
        with open(new_file, "w") as g:
            g.write("#OTUID\ttaxonomy\tconfidence\n")
            for line in f:
                g.write(line)
    return



##### 2. Platform Specific #####
def arg_seq_relabeler(post_QC, relabeled):
    """Relabels sequences labeled with Illumina information and QIIME's split_libraries_fastq headers.
    Uses the mapping file to match barcodes with sample IDs.
    """
    with open(post_QC, 'r') as f:
        with open(relabeled, 'w') as g:
            for line in f:
                if line.startswith('>'):
                    header = line.split(' ')[0]
                    bc = line.split('orig_bc=')[1].split(' ')[0]
                    header = header.rsplit('_',1)
                    new_header = header[0].replace('_', '.')
                    line = new_header + '_' + header[-1]+ ';barcodelabel='+new_header[1:] +"\n"
                g.write(line)


def uic_seq_relabeler(paths):
    """ 
    Relabels sequences for UIC style labels and concatenates output into a
    single file called combined_seqs.fna
    """
    print ("Relabeling & concatenating reads!")
    
    with open(paths['relabeled_fasta'], 'w') as g:
        for i in glob.glob(paths["uic_filtered"]+"*"):
            ind = 0
            out = i.split('/')[-1].split("_filtered")[0]
            
            with open(i, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        barcode = line.rsplit(':',1)[-1].strip()
                        line = ">barcodelabel={0};seq_{1} \
                        barcode={2}\n".format(out.replace("_","."),str(ind),barcode)
                        ind+=1
                        g.write(line)
                    else:
                        ind+=1
                        g.write(line)
    
    print ("Relabeling & concatenating reads done!")
    return

def pear_files(fwd_reverse, base_dir, threads="3"):
    """ Takes a directory of forward and reverse reads that start with identical naming conventions and are listed as:
        f1,r1,f2,r2... 
        and runs PEAR on them using the settings specified by threads, lmin and lmax. """
    # Iterate thru files and get forward and reverse, then call pear. 
    count=0
    int_assembled = base_dir+'int_assembled/'
    ensure_path(int_assembled) # Make intermediate

    # Call pear on each set of forward and reverse sequences assuming they all have the same naming convention.
    for i in glob.glob(fwd_reverse+'*'): 
            if count == 0:
                short_name = i.split('/')[-1].split('-')[0]  # Need to split this functionality out or standardize names.
                forward = i
            if count == 1:
                reverse = i
                subprocess.call(["pear", "-f", forward, "-r", reverse, "-o", int_assembled+short_name, "-j", threads])
            count+= 1
            count= count%2

    # Move paired reads to a new directory
    assembled_dir = base_dir+"assembled/"
    ensure_path(assembled_dir)
    assembled_files = glob.glob(int_assembled+"*")
    for i in range(len(assembled_files)):
        if i%4==0:
            shutil.copyfile(assembled_files[i], assembled_dir+assembled_files[i].split("/")[-1])
    return


##### 3. Directory Setup #####
def ensure_path(directory):
    """ Makes the specified directory if it does not exist."""
    if not os.path.exists(directory):
        print ("{0} Did not exist, making empty directory".format(directory))
        os.makedirs(directory)
    return

def load_databases():
    db_dict = {
            "rdp_database": "/Users/jimbo/databases/rdp_16s.fa",
            "taxa_confs" : "/Users/jimbo/databases/rdp_16s_fl.tc",
            "tax_tt" : "/Users/jimbo/databases/rdp_16s.tt",
            "gg_tax_db" : "/macqiime/greengenes/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt",
            "gg_aligned" :
            "/macqiime/greengenes/core_set_aligned.fasta.imputed",
            
            "silva99_tax" :
            "/Users/jimbo/databases/Silva_128/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt",
            
            "silva99_ref" :
            "/Users/jimbo/databases/Silva_128/rep_set/rep_set_16S_only/99/99_otus_16S.fasta",
            
            "silva_euk_ref" :
            "/Users/jimbo/databases/Silva_128/rep_set/rep_set_18S_only/99/99_otus_18S.fasta",

            "silva_euk_tax" :
            "/Users/jimbo/databases/Silva_128/taxonomy/18S_only/99/majority_taxonomy_7_levels.txt"
            }
    return db_dict

def construct_directory_structure(base_dir):
    ensure_path(base_dir)
    vsearch = base_dir+"vsearch/"
    paired = base_dir + "paired"
    split_libraries = base_dir+'split_libraries/'
    paths = {
            "vsearch" : vsearch,  # could probably do with just vsearch above?
            "logging" : base_dir + 'logging/',
            "paired" : paired,    # 
            "assembled" : base_dir + 'assembled/',
            "paired_barcodes" : paired+"/fastqjoin.join_barcodes.fastq",
            "split_libraries" :base_dir+"split_libraries/",
            "relabeled" : split_libraries + 'seqs_with_barcode_labels.fastq',
            'relabeled_fasta' : vsearch+'relabeled_fna.fasta',
            "derep_mc2" : vsearch+"derep_mc2.fa",
            "otu_centroids" : vsearch+"centroids.fa",
            "otu_clusters" : vsearch+"uc_clusters.uc",
            "fasta_seqs" : split_libraries+'seqs.fna',
            "fastq_seqs" : split_libraries+'seqs.fastq',
            "filtered" : split_libraries + 'filtered.fasta',
            "no_adapter" : split_libraries + "no_adapters.fastq",
            "uic_filtered" : base_dir + "filtered/",
            "combined_seqs" : vsearch + "combined_seqs.fna",
# Chimera detection
            "chimeras_dn" : vsearch+"chimeras_dn.fa",
            "non_chimeras_dn" : vsearch+"non_chimeras_dn.fa",
            "dn_log" : vsearch + "chimeras_dn.log",
# Mapped reads
            "otu_map" : vsearch +'otu_map.uc',
# tax files and biom files
            "tax" : vsearch+"assigned_tax",
            "otu_table" : vsearch+"otu_table.txt",
            "biom_table" : vsearch+"biom_table.biom",
            "biom_wtax" : vsearch+"biom_table_w_tax.biom",
            "txt_biom_w_tax" : vsearch + "biom_table_w_tax.txt",
            "gg_tax_summary" : vsearch+"gg_tax_summary/",
# Alignment and tree construction
            "aligned" : vsearch+'aligned/'
            }
    ensure_path(vsearch)
    ensure_path(paths["logging"])
    ensure_path(split_libraries)  # or paths['split_libraries?
    ensure_path(paths["tax"])
    ensure_path(paths["gg_tax_summary"])
    return paths

##### 4. Cosmetic #####
def print_name():
    lines=[
    "                   ___  ______  _____  _   _____              ",
    "       _~        /  _ \ | ___ \|  __ \| | /  ___|             ",
    "    _~ )_)_~    /  /_\ \| |_/ /| |  \/| | \ `--.   ___   __ _ ",
    "    )_))_))_)   |   _  ||    / | | __ | |  `--. \ / _ \ / _` |",
    "    _!__!__!_   |  | | || |\ \ | |_\ \|_| /\__/ /|  __/| (_| |",
    "    \______t/    \_| |_/\_| \_| \____/(_) \____/  \___| \__, |",
    "                                                           | |",
    "                                                           |_|"]
    for i in lines:
        print (i)
    return

def print_name_uic():
    lines=[" __    __   __       _______. _______   ______",
    "|  |  |  | |  |     /       ||   ____| /  __  \     ",
    "|  |  |  | |  |    |   (----`|  |__   |  |  |  |    ",
    "|  |  |  | |  |     \   \    |   __|  |  |  |  |    ",
    "|  `--'  | |  | .----)   |   |  |____ |  `--'  '--. ",
    " \______/  |__| |_______/    |_______| \_____\_____\ "]
    for i in lines:
        print (i)
    return

def version_info(version, author, date):
    print ( "Written by {0}\nVersion {1}\nLast Modified {2}".format(
        author, version, date
        )
    )
    return



