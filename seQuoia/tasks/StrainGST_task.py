import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.classes.KmerAnalyzer import Kmer
from seQuoia.other import usefulFunctions as uF

def runStrainGST(fastq_sin, fastq_frw, fastq_rev, db, sample_name, parent_dir, options_kmerize, options_straingst):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev))))
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    try: assert(os.path.isfile(db))
    except: sys.stderr.write("ERROR: StrainGST pangenome database file is not available."); raise RuntimeError()

    # set up directory structure
    workspace_name = "StrainGST"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'StrainGST.log'
    logObject = uF.createLoggerObject(log_file)

    kmer_file = None
    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # run strainge kmerize
        kmer_file = FastqObj.kmerize(workspace, options=options_kmerize)

    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:

        # create Fastq object
        FastqObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # run strainge kmerize
        kmer_file = FastqObj.kmerize(workspace, options=options_kmerize)

    # create Kmer object
    KmerObj = Kmer(kmer_file, sample_name, logObject)

    # run straingst
    KmerObj.run_straingst(workspace, db, options=options_straingst)

    # produce kmer histogram - in progress - issues running.
    # KmerObj.create_histogram(workspace)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "STRAINGST.txt", 'w')
    conf_file.write("StrainGST: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running StrainGST analysis to find closest reference strain and performing k-merization of reads.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-d', '--db', help="location of target StrainGST pangenome database.", required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-pk', '--options_kmerize', help='options for running StrainGE kmerize.', required=False, default="-k 23")
    parser.add_argument('-ps', '--options_straingst', help='options for running StrainGST', required=False, default="")
    args = parser.parse_args()

    runStrainGST(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.db, args.sample, args.sample_dir, args.options_kmerize, args.options_straingst)