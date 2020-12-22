import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.other import usefulFunctions as uF

def runSubsample(fastq_sin, fastq_frw, fastq_rev, sample_name, parent_dir, reads, bases, no_gzip):
    try: assert( (fastq_sin and os.path.isfile(fastq_sin)) or (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev)))
    except: sys.stderr.write("ERROR: FASTQ inputs were not provided. Please provide either a pair (frw and rev) or a single FASTQ file.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Subsample"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Subsample.log'
    logObject = uF.createLoggerObject(log_file)

    ### Perform single end QC analysis
    if fastq_sin and os.path.isfile(fastq_sin):

        # create Fastq object
        FastqObj = Fastq(fastq_sin, sample_name, logObject)

        # run FastQC and parse results.
        if reads:
            FastqObj.subsample(workspace, reads=reads, compress=(not no_gzip))
        elif bases:
            FastqObj.downsample(workspace, bases=bases, compress=(not no_gzip))
        else:
            logObject.error("No subsampling quantity specified defaulting to 100K reads being subsampled!")
            FastqObj.subsample(workspace, reads=reads, compress=(not no_gzip))

    ### Perform paired-end QC analysis
    elif fastq_frw and fastq_rev:

        # create FastqPaired Object
        FastqPairedObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)

        # run FastQC and parse results.
        if reads:
            FastqPairedObj.subsample(workspace, reads=reads, compress=(not no_gzip))
        elif bases:
            FastqPairedObj.downsample(workspace, bases=bases, compress=(not no_gzip))
        else:
            logObject.error("No subsampling quantity specified defaulting to 100K ")
            FastqPairedObj.subsample(workspace, reads=reads, compress=(not no_gzip))

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "SUBSAMPLE.txt", 'w')
    conf_file.write("Subsample: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    seQc module for subsampling FASTQ file to have reduced number of reads. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-f', '--fastq_sin', help='location of input single end FASTQ file.', required=False, default="")
    parser.add_argument('-frw', '--fastq_frw', help="location of forward FASTQ file.", required=False, default="")
    parser.add_argument('-rev', '--fastq_rev', help="location of reverse FASTQ file.", required=False, default="")
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-n', '--reads', type=int, help='specify the number of reads to sample. Deafult is 100,000', required=False, default=None)
    parser.add_argument('-m', '--bases', type=int, help='specify the number of bases to sample.', required=False, default=None)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses or leaves compressed).', required=False, default=False)
    args = parser.parse_args()

    runSubsample(args.fastq_sin, args.fastq_frw, args.fastq_rev, args.sample, args.sample_dir, args.reads, args.bases, args.no_gzip)