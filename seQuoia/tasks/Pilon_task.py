import os
import sys
import argparse
from seQuoia.classes.AlignmentAnalyzer import Alignment
from seQuoia.other import usefulFunctions as uF

def runPilon(bam_file, reference_fasta, sample_name, parent_dir, pilon_options, cores):
    try: assert( os.path.isfile(bam_file) )
    except: sys.stderr.write("ERROR: BAM input was not provided. Please provide. Exiting now ...\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Pilon"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Pilon.log'
    logObject = uF.createLoggerObject(log_file)

    #### Start Pilon workflow

    # create Alignment object
    AlignmentObj = Alignment(bam_file, sample_name, logObject)

    # run Pilon
    AlignmentObj.run_pilon(workspace, reference_fasta, options=pilon_options)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "PILON.txt", 'w')
    conf_file.write("Pilon: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Run Pilon for variant calling.
    """)

    parser.add_argument('-b', '--bam_file', help="location of BAM file. Should be pre-processed already.", required=True)
    parser.add_argument('-r', '--reference_fasta', help="location of reference FASTA assembly.", required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-p', '--pilon_options', help='BWA options.', required=False, default="--variant")
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    args = parser.parse_args()

    runPilon(args.bam_file, args.reference_fasta, args.sample, args.sample_dir, args.pilon_options, args.cores)