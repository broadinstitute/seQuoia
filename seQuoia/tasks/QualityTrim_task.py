import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runTrimmomatic(fastq_sin, fastq_frw, fastq_rev, sample_name, parent_dir, trimmomatic_options, cores, no_gzip):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    trimmomatic_options = trimmomatic_options.strip('"')

    # set up directory structure
    workspace_name = "QualityTrim"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'QualityTrim.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end trimmomatic operation
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # trim adapters using cutadapt and return resulting FASTQ file in gzip compressed format
        FastqObj.quality_trim(workspace, options=trimmomatic_options, cores=cores, compress=(not no_gzip))

    ### Perform paired-end trimmomatic operation
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # trim adpaters using cutadapt and return resulting FASTQ files in gzip compressed format
        FastqPairedObj.quality_trim(workspace, options=trimmomatic_options, cores=cores, compress=(not no_gzip))

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "QUALITYTRIM.txt", 'w')
    conf_file.write("QualityTrim: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running Trimmomatic quality trimming. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-c', '--cores', help='select number of threads. default is 1.', required=False, default=1)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    parser.add_argument('-p', '--trimmomatic_options', help="Options for cutadapt which specify adapters to cut. Default is to use the following Trimmomatic settings: LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15",
                        default = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15", required=False)
    args = parser.parse_args()

    runTrimmomatic(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.trimmomatic_options, args.cores, args.no_gzip)