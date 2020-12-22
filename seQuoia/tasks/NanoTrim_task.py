import os
import sys
import argparse
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.other import usefulFunctions as uF

def runNanoTrim(nanopore_fastq, sample_dir, sample_name, no_gzip):
    try: assert( nanopore_fastq and (os.path.isfile(nanopore_fastq) or os.path.isdir(nanopore_fastq)) )
    except:
        sys.stderr.write("ERROR: FASTQ input(s) were not provided. Please provide. Raising exception\n")
        raise RuntimeError

    # set up directory structure
    workspace_name = "NanoTrim"
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'NanoTrim.log'
    logObject = uF.createLoggerObject(log_file)

    # Initialize Nanopore Object
    NanoporeObj = Nanopore(nanopore_fastq, sample_name, logObject)

    # Trim any adapters with PoreChop
    NanoporeObj.run_nanotrim(workspace, compress=(not no_gzip))

    conf_file = open(sample_dir + "NANOTRIM.txt", 'w')
    conf_file.write("NanoTrim: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for trimming adapters from nanopore fastq using PoreChop. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-f', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    args = parser.parse_args()

    runNanoTrim(args.nanopore_fastq, args.sample_dir, args.sample, args.no_gzip)