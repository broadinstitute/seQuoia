import os
import sys
import argparse
from seQuoia.other import usefulFunctions as uF
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired


def runCentrifuge(fastq_sin, fastq_frw, fastq_rev, centrifuge_index, sample_name, parent_dir, cores):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Centrifuge"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Centrifuge.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # bin reads taxonomically using Centrifuge
        FastqObj.bin_taxonomically(workspace, centrifuge_index, cores=cores)

    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # bin reads taxonomically using Centrifuge
        FastqPairedObj.bin_taxonomically(workspace, centrifuge_index, cores=cores)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "CENTRIFUGE.txt", 'w')
    conf_file.write("Centrifuge: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running centrifuge to get taxonomic profile of data.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-x', '--centrifuge_index', help="Centrifuge index. Default is located at: ", required=True)
    parser.add_argument('-s', '--sample', help="sample name.", required=False, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    args = parser.parse_args()

    runCentrifuge(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.centrifuge_index, args.sample, args.sample_dir, args.cores)