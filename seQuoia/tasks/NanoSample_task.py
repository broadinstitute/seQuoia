import os
import sys
import argparse
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.other import usefulFunctions as uF

def runNanoSample(nanopore_fastq, sample_dir, sample_name, fastqfilter_options, no_gzip):
    try: assert( nanopore_fastq and os.path.isfile(nanopore_fastq) )
    except:
        raise RuntimeError("ERROR: FASTQ input(s) were not provided. Please provide. Raising exception\n")

    fastqfilter_options = fastqfilter_options.strip('"')

    # set up directory structure
    workspace_name = "NanoSample"
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'NanoSample.log'
    logObject = uF.createLoggerObject(log_file)

    # Initialize Nanopore Object
    NanoporeObj = Nanopore(nanopore_fastq, sample_name, logObject)

    # Subsample Nanopore reads using fastqfilter by Bruce Walker
    NanoporeObj.run_fastqfilter(workspace, options=fastqfilter_options, compress=(not no_gzip))

    conf_file = open(sample_dir + "NANOSAMPLE.txt", 'w')
    conf_file.write("NanoSample: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running subsampling of Nanopore reads. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-f', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-p', '--fastqfilter_parameters', help='filtering parameters for FASTQ filter.', required=False, default='-l 3000 -b 300000000')
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    args = parser.parse_args()

    runNanoSample(args.nanopore_fastq, args.sample_dir, args.sample, args.fastqfilter_parameters, args.no_gzip)