import os
import sys
import argparse
from seQuoia.other import usefulFunctions as uF
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired

def runSymlinkInput(fastq_sin, fastq_frw, fastq_rev, parent_dir, sample_name):
        try:
            assert ((fastq_sin and os.path.isfile(fastq_sin)) or (
            fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev)))
        except:
            sys.stderr.write(
                "ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

        # set up directory structure
        workspace_name = "Symlink_Input"
        workspace = uF.setupDirectory(parent_dir, workspace_name)

        # create logging object
        log_file = workspace + 'Symlink.log'
        logObject = uF.createLoggerObject(log_file)

        ### Perform single end QC analysis
        if fastq_sin and os.path.isfile(fastq_sin):

            # create Fastq object
            FastqObj = Fastq(fastq_sin, sample_name, logObject)

            # create symlink
            FastqObj.create_symlink(workspace)

        ### Perform paired-end QC analysis
        elif fastq_frw and fastq_rev:

            # create FastqPaired Object
            FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

            # create symlink
            FastqPairedObj.create_symlink(workspace)
        conf_file = open(parent_dir + "SYMLINK_INPUT.txt", 'w')
        conf_file.write("SymlinkInput: Module Completed Succesfully!")
        conf_file.close()


if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for just symlinking input FASTQ files with better naming for pipeline support.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    args = parser.parse_args()

    runSymlinkInput(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample_dir, args.sample)