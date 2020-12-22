import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runAdapterTrim(fastq_sin, fastq_frw, fastq_rev, sample_name, parent_dir, trimgalore_options, cutadapt_options):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    trimgalore_options = trimgalore_options.strip('"')
    cutadapt_options = cutadapt_options.strip('"')

    try: assert(not (trimgalore_options and cutadapt_options))
    except: sys.stderr.write("ERROR: Both filtering options with cutadapt and trim galore provided. Can only use one adapter trimmer. Exiting now ...\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "AdapterTrim"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'AdapterTrim.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end cutadapt operation
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # trim adpaters using trim-galore preset settings for nextera and return resulting FASTQ files in gzip compressed format
        if cutadapt_options:
            FastqObj.cutadapt_adapter_trim(workspace, options=cutadapt_options)
        else:
            FastqObj.trim_galore_adapter_trim(workspace, options=trimgalore_options)

    ### Perform paired-end cutadapt operation
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # trim adpaters using trim-galore preset settings for nextera and return resulting FASTQ files in gzip compressed format
        if cutadapt_options:
            FastqPairedObj.cutadapt_adapter_trim(workspace, options=cutadapt_options)
        else:
            FastqPairedObj.trim_galore_adapter_trim(workspace, options=trimgalore_options)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "ADAPTERTRIM.txt", 'w')
    conf_file.write("AdapterTrim: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running adapter trimming. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help="sample name.", required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-co', '--cutadapt_options', help="Options for cutadapt. Species to use cutadapt. [default]", required=False, default="")
    parser.add_argument('-to', '--trimgalore_options', help="Options for trim_galore which specify adapters to cut. Specifies to use trim galore", required=False, default="")
    args = parser.parse_args()

    runAdapterTrim(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.trimgalore_options, args.cutadapt_options)