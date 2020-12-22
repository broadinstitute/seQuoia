import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.classes.AlignmentAnalyzer import Alignment
from seQuoia.other import usefulFunctions as uF

def runRefAlignment(fastq_sin, fastq_frw, fastq_rev, reference_fasta, sample_name, parent_dir, bwa_options, cores):
    try: assert( (os.path.isfile(fastq_sin)) or (os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev)) )
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    try: assert(os.path.isfile(reference_fasta))
    except: sys.stderr.write("ERROR: Reference FASTA file does not have the correct format.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "ReferenceAlignment"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'ReferenceAlignment.log'
    logObject = uF.createLoggerObject(log_file)

    sam_file = None
    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # Align reads to reference genome
        sam_file = FastqObj.align_to_reference(workspace, reference_fasta, options=bwa_options, cores=cores)


    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # Align reads to reference genome
        sam_file = FastqPairedObj.align_to_reference(workspace, reference_fasta, options=bwa_options, cores=cores)

    # create Alignment object
    AlignmentObj = Alignment(sam_file, sample_name, logObject)

    # compress SAM to BAM
    AlignmentObj.compress_sam(workspace, clean=True)

    # sort BAM file
    AlignmentObj.sort_bam(workspace, clean=True)

    # index BAM file
    AlignmentObj.index_bam(workspace)

    # mark duplicates
    AlignmentObj.mark_dups(workspace, clean=True)

    # index BAM file
    AlignmentObj.index_bam(workspace)

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "REFALIGNMENT.txt", 'w')
    conf_file.write("Reference Alignment: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Run reference alignment of reads and perform some downstream processing.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False)
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False)
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False)
    parser.add_argument('-r', '--reference_fasta', help="location of reference FASTA assembly.", required=False, default="")
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-b', '--bwa_options', help='BWA options.', required=False, default="-M")
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    args = parser.parse_args()

    runRefAlignment(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.reference_fasta, args.sample, args.sample_dir, args.bwa_options, args.cores)