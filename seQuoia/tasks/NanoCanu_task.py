import os
import argparse
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.classes.FastqAnalyzer import Fastq
from seQuoia.other import usefulFunctions as uF

def runNanoCanu(nanopore_fastq, sample_dir, sample_name, canu_options, memory, cores):
    try: assert( nanopore_fastq and os.path.isfile(nanopore_fastq) )
    except:
        raise RuntimeError("ERROR: FASTQ input(s) were not provided properly. Please fix. Raising exception\n")

    unicycler_options = canu_options.strip('"')

    # set up directory structure
    workspace_name = "Canu_Assembly"
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Canu_Assembly.log'
    logObject = uF.createLoggerObject(log_file)

    if nanopore_fastq.endswith('.gz'):
        FastqObj = Fastq(nanopore_fastq, sample_name, logObject)
        FastqObj.create_new_instance(workspace, compress=False, change_reference=True)

        # Initialize Nanopore Object
        NanoporeObj = Nanopore(FastqObj.fastq, sample_name, logObject)

        # Run Canu for assembly
        NanoporeObj.run_canu(workspace, options=canu_options, memory=memory, cores=cores)

        # Clean up temporary FASTQ instance
        os.system('rm -f %s' % FastqObj.fastq)

    else:
        # Initialize Nanopore Object
        NanoporeObj = Nanopore(nanopore_fastq, sample_name, logObject)

        # Run Canu for assembly
        NanoporeObj.run_canu(workspace, options=canu_options, memory=memory, cores=cores)


    conf_file = open(sample_dir + "CANU_ASSEMBLY.txt", 'w')
    conf_file.write("Canu Assembly: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running Canu for nanopore only assembly. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-f', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-p', '--canu_parameters', help='parameters for Canu.', required=False, default="stopOnReadQuality=false genomeSize=5m")
    parser.add_argument('-m', '--memory', type=int, help='Gb of memory needed.', required=False, default=24)
    parser.add_argument('-c', '--cores', type=int, help='number of threads.', required=False, default=1)

    args = parser.parse_args()

    runNanoCanu(args.nanopore_fastq, args.sample_dir, args.sample, args.canu_parameters, args.memory, args.cores)