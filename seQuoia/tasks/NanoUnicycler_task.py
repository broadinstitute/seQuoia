import os
import argparse
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.other import usefulFunctions as uF

def runNanoUnicycler(nanopore_fastq, illumina_forward, illumina_reverse, sample_dir, sample_name, unicycler_options, identifier, cores):
    try: assert( nanopore_fastq and os.path.isfile(nanopore_fastq) )
    except:
        raise RuntimeError("ERROR: FASTQ input(s) were not provided properly. Please fix. Raising exception\n")

    try: assert( illumina_forward and illumina_reverse and os.path.isfile(illumina_forward) and os.path.isfile(illumina_reverse) )
    except:
        raise RuntimeError("ERROR: Optional Illumina FASTQ input(s) / assembly were not provided properly. Please fix. Raising exception\n")

    unicycler_options = unicycler_options.strip('"')

    # set up directory structure
    workspace_name = "Unicycler_Assembly"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Unicycler_Assembly.log'
    logObject = uF.createLoggerObject(log_file)

    # initialize Nanopore object
    NanoporeObj = Nanopore(nanopore_fastq, sample_name, logObject)

    # Run Unicycler for assembly
    NanoporeObj.run_unicycler(illumina_forward, illumina_reverse, workspace, options=unicycler_options, cores=cores)

    conf_file_name = sample_dir + "UNICYCLER_ASSEMBLY"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("Unicycler Assembly: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running Unicycler for nanopore + illumina hybrid assembly. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-f', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-1', '--illumina_forward', help='location of forward illumina reads for sample.', required=True)
    parser.add_argument('-2', '--illumina_reverse', help='location of reverse illumina reads for sample.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-p', '--unicycler_parameters', help='parameters for Unicycler.', required=False, default="--mode normal --verbosity 2")
    parser.add_argument('-i', '--identifier', help='identifier, in case multiple runs in workflow.', required=False, default=None)
    parser.add_argument('-c', '--cores', type=int, help='number of threads.', required=False, default=1)

    args = parser.parse_args()

    runNanoUnicycler(args.nanopore_fastq, args.illumina_forward, args.illumina_reverse, args.sample_dir, args.sample, args.unicycler_parameters, args.identifier, args.cores)