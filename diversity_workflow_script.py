#!/Users/jimbo/anaconda/bin/python
# This script runs the alpha and beta diversity analysis workflow scripts from QIIME1.

### Author and version information:
version = "1.0"
author = "JGrif"
date = "10/20/17"

import Arg_seq_functions as ll
import argparse

def main():
    args = parse_arguments()
    paths = ll.construct_directory_structure(args.base_dir)
    tree = paths["vsearch"]+"phylo_tree.tre"
    alpha_dir = paths["vsearch"] + "alpha"
    beta_dir = paths["vsearch"] + "beta"

    
    if args.alpha:
        ll.run_command("alpha_rarefaction.py -i {0} -m {1} -o {2} -t {3} -e {4} -O 3 -a -f".format(
        paths["biom_wtax"], args.m, alpha_dir, tree, args.rare ), "Alpha Rarefaction"
        )

    if args.beta:
        ll.run_command("beta_diversity_through_plots.py -i {0} -m {1} -o {2} -t {3} -e {4} -a -O 3 -f".format(
            paths["biom_wtax"], args.m, beta_dir, tree, args.rare), "Beta Div Workflow"
            )

# Functions
def parse_arguments():

    # Parse Command Line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--rare", help = "Max value to rarefy samples to")
    parser.add_argument("--base_dir", help="Base Directory for sequence analysis workflow script")
    parser.add_argument("--log", help = "Logfile to output data to. If blank, output to terminal only.")
    parser.add_argument("--beta", help = "Run QIIME beta diversity thru plots", action='store_true')
    parser.add_argument("--alpha", help = "Run QIIME alpha rarefaction", action='store_true')
    parser.add_argument("--m", help = "map file")
    
    try:
        result = parser.parse_args()
        return result
    except Exception as e: 
        parser.print_help()
        print(e)
        sys.exit(0)

if __name__ == "__main__":
    main()
