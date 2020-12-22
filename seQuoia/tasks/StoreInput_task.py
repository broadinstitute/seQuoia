import os
import sys
import argparse
from seQuoia.other import usefulFunctions as uF

def storeInput(fastq_sin, fastq_frw, fastq_rev, sample_name, parent_dir):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin)) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev)))
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "LSARP_Results/"
    workspace = uF.setupDirectory(parent_dir, workspace_name)
    workspace_name_input = 'Input/'
    workspace_input = uF.setupDirectory(workspace, workspace_name_input)

    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):
        res_sin_read = workspace_input + sample_name + '_R1.processed.fastq'
        if fastq_sin.endswith('.gz'): res_sin_read += '.gz'
        os.system('cp %s %s' % (fastq_sin, res_sin_read))

    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:
        res_frw_read = workspace_input + sample_name + '_R1.processed.fastq'
        res_rev_read = workspace_input + sample_name + '_R2.processed.fastq'
        if fastq_frw.endswith('.gz'): res_frw_read += '.gz'
        if fastq_rev.endswith('.gz'): res_rev_read += '.gz'
        os.system('cp %s %s' % (fastq_frw, res_frw_read))
        os.system('cp %s %s' % (fastq_rev, res_rev_read))

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "STOREINPUT.txt", 'w')
    conf_file.write("StoreInput: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for storing Illumina data for LSARP samples.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help="sample name.", required=True, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    args = parser.parse_args()

    storeInput(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir)