import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runAssembly(fastq_frw, fastq_rev, sample_name, parent_dir, read_length, unicycler, cores):
    try: assert( os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Assembly"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Assembly.log'
    logObject = uF.createLoggerObject(log_file)

    # create FastqPaired Object
    FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

    if unicycler: FastqPairedObj.run_unicycler(workspace, cores=cores)
    else: FastqPairedObj.run_spades(workspace, read_length=read_length, cores=cores)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "ASSEMBLY.txt", 'w')
    conf_file.write("Assembly: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running Illumina-only assembly. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=True)
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-rl', '--read_length', type=int, help='read lengths in sequencing data.', required=False, default=150)
    parser.add_argument('-u', '--unicycler', action='store_true', help='use Unicycler for extra perks and pilon polishing.', default=False, required=False)
    parser.add_argument('-c', '--cores', type=int, help='specify number of cores.', required=False, default=1)
    args = parser.parse_args()

    runAssembly(args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.read_length, args.unicycler, args.cores)