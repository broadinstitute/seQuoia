import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runMLST(fastq_frw, fastq_rev, sample_name, parent_dir, ariba_db_dir):
    try: assert(fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    try: assert(os.path.isdir(ariba_db_dir))
    except: sys.stderr.write("ERROR: Some issue occurred with provided databases/options. Please check the input and retry.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "MLST"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'MLST.log'
    logObject = uF.createLoggerObject(log_file)

    # create FastqPaired Object
    FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

    # Run ARIBA analysis to detect STs in raw reads
    FastqPairedObj.ariba(workspace, name='ariba_mlst', ariba_db=ariba_db_dir)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "MLST.txt", 'w')
    conf_file.write("MLST: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running multi-locus subtyping analysis. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help="sample name.", required=False, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-ad', '--ariba_db_dir', help='MLST database for organism of interest for ARIBA analysis.', required=False)
    args = parser.parse_args()

    runMLST(args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.ariba_db_dir)