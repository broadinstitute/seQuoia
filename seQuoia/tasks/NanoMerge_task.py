import os
import sys
import argparse
from seQuoia.other import nanoporeFunctions as npF
from seQuoia.other import usefulFunctions as uF

def runNanoMerge(nanopore_fastq, sample_dir, sample_name, barcode, no_gzip):
    try: assert( nanopore_fastq and (os.path.isfile(nanopore_fastq) or os.path.isdir(nanopore_fastq)) )
    except:
        sys.stderr.write("ERROR: FASTQ input(s) were not provided. Please provide. Raising exception\n")
        raise RuntimeError

    # set up directory structure
    workspace_name = "NanoMerge"
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'NanoMerge.log'
    logObject = uF.createLoggerObject(log_file)

    if os.path.isdir(nanopore_fastq):
        npF.concat_fastqs(nanopore_fastq, workspace, sample_name, logObject, barcode=barcode, compress=(not no_gzip))
    elif os.path.isfile(nanopore_fastq):
        if nanopore_fastq.endswith('.gz'):
            os.system('cp %s %s' % (nanopore_fastq, workspace + sample_name + '.fastq.gz'))
            if no_gzip:
                os.system('gunzip %s' % workspace + sample_name + '.fastq.gz')
        else:
            os.system('cp %s %s' % (nanopore_fastq, workspace + sample_name + '.fastq'))
            if not no_gzip:
                os.system('gzip %s' % workspace + sample_name + '.fastq')

    conf_file = open(sample_dir + "NANOMERGE.txt", 'w')
    conf_file.write("NanoMerge: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running merging of multiple nanopore FASTQ files. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-f', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-b', '--barcode', help='provide barcode. This will make it so that /barcode_id/ is needed in the file path to be merged and regarded as belonging to the sample.', required=False)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    args = parser.parse_args()

    runNanoMerge(args.nanopore_fastq, args.sample_dir, args.sample, args.barcode, args.no_gzip)