import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runSortMeRNA(fastq_sin, fastq_frw, fastq_rev, database_dir, sample_name, parent_dir, cores, no_gzip):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "SortMeRNA"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'SortMeRNA.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # split up ribosomal and non-ribosomal RNA data
        FastqObj.filter_ribo_rna(workspace, database_dir, cores=cores, compress=(not no_gzip))

    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # split up ribosomal and non-ribosomal RNA data
        FastqPairedObj.filter_ribo_rna(workspace, database_dir, cores=cores, compress=(not no_gzip))

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "SORTMERNA.txt", 'w')
    conf_file.write("SortMeRNA: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running removing rRNA reads. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-d', '--database_dir', help="path to sortmerna primary/root database directory.", required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    args = parser.parse_args()

    runSortMeRNA(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.database_dir, args.sample, args.sample_dir, args.cores, args.no_gzip)