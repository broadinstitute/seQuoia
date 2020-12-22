import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runBayesHammer(fastq_sin, fastq_frw, fastq_rev, sample_name, parent_dir):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "BayesHammer"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'BayesHammer.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end cutadapt operation
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)
        FastqObj.error_correction(workspace)

    ### Perform paired-end cutadapt operation
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # trim adpaters using trim-galore preset settings for nextera and return resulting FASTQ files in gzip compressed format
        FastqPairedObj.error_correction(workspace)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "BAYESHAMMER.txt", 'w')
    conf_file.write("BayesHammer: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running BayesHammer read correction module. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help="sample name.", required=False, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    args = parser.parse_args()

    runBayesHammer(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir)