import os
import sys
import argparse
from seQuoia.other import usefulFunctions as uF
from seQuoia.other import nanoporeFunctions as npF

def runNanoQC(nanopore_fastq, nanopore_seqsum, nanopore_barcode, sample_dir, cores):
    try: assert( nanopore_fastq and (os.path.isfile(nanopore_fastq)) )
    except:
        sys.stderr.write("ERROR: FASTQ input(s) were not provided. Please provide. Raising exception\n")
        raise RuntimeError

    # set up directory structure
    workspace_name = "NanoQC"
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'NanoQC.log'
    logObject = uF.createLoggerObject(log_file)

    npF.run_nanoplot_qc(nanopore_fastq, workspace, logObject, cores=cores)
    if os.path.isfile(nanopore_seqsum):
        sample_seqsum = npF.filter_sequence_summary(nanopore_seqsum, nanopore_barcode, workspace, logObject)
        npF.run_minion_qc(sample_seqsum, workspace, logObject)

    # create successful completion file if steps completed!
    conf_file = open(sample_dir + "NANOQC.txt", 'w')
    conf_file.write("NanoQC: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running NanoQC/minion_qc quality control assessment. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-fq', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True, default="")
    parser.add_argument('-ss', '--nanopore_seqsum', help="location of sequencing_summary.txt file produced by basecallers Guppy/Albacore.", required=False, default="")
    parser.add_argument('-bc', '--nanopore_barcode', help="sample barcode as referenced in sequencing summary file.", required=False, default="")
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    args = parser.parse_args()

    runNanoQC(args.nanopore_fastq, args.nanopore_seqsum, args.nanopore_barcode, args.sample_dir, args.cores)