import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runAMRP(fastq_frw, fastq_rev, sample_name, parent_dir, shortbred_markers, ariba_database, ariba_names):
    try: assert( fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    try: assert(shortbred_markers or ariba_database)
    except: sys.stderr.write("ERROR: Some issue occurred with provided databases/options. Please check the input and retry.\n"); raise RuntimeError

    try:
        if ariba_database or ariba_names:
            assert(ariba_database and ariba_names)
    except: sys.stderr.write("ERROR: ARIBA database provided without names or visa versa, either way please check the input and retry!.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "AMRP"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'AMPR.log'
    logObject = uF.createLoggerObject(log_file)

    # create FastqPaired Object
    FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

    # Run AMR Prediction analysis
    if ariba_database:
        for i, adb in enumerate(ariba_database):
            adb_name = ariba_names[i]
            FastqPairedObj.ariba(workspace, name=adb_name, ariba_db=adb)
    if shortbred_markers:
        FastqPairedObj.shortbred_amrp(workspace, shortbred_markers)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "AMRP.txt", 'w')
    conf_file.write("AMRP: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for identifying AMR genes in samples from reads. Uses FastqAnalyzer.py class OOP framework. Currently, implements the ability to find AMR genes using 
    ARIBA and/or shortBRED.
    """)
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-sb', '--shortbred_markers', help="Protein markers prepared using shortBRED identify.", required=False)
    parser.add_argument('-ad', '--ariba_database', nargs='+', help="Prepared reference for ARIBA.", required=False)
    parser.add_argument('-an', '--ariba_names', nargs='+', help="Reference identifiers for ARIBA", required=False)

    args = parser.parse_args()

    runAMRP(args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.shortbred_markers,
            args.ariba_database, args.ariba_names)