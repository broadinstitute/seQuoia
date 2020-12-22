import os
import sys
import argparse
import shutil
from seQuoia.other import usefulFunctions as uF
from seQuoia.classes.AlignmentAnalyzer import Alignment

def processGpDirectory(bam_location, sample_name, parent_dir, no_gzip):
    try: assert(os.path.isdir(bam_location) or os.path.isfile(bam_location))
    except: sys.stderr.write("ERROR: BAM/GP directory does not exist! Exiting now ...\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "ProcessGPDirectory"
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'ProcessGPDirectory.log'
    logObject = uF.createLoggerObject(log_file)

    bam_location = os.path.abspath(bam_location)

    logObject.info('*' * 70)
    logObject.info("Beginning to convert BAM location %s to FASTQ(s)" % (os.path.abspath(bam_location) + '/'))

    input_bam = bam_location
    if os.path.isdir(bam_location):
        gp_directory = bam_location + '/'

        # copy over all txt files
        metric_files = [gp_directory + f for f in os.listdir(gp_directory) if not f.endswith('.pdf') and not f.endswith('.bam') and not f.endswith('.bai') and not os.path.isdir(gp_directory + f)]
        for mf in metric_files:
            mf_basename = mf.split('/')[-1]
            try:
                 shutil.copy(mf, workspace + mf_basename)
            except:
                 pass

        bam_files = [gp_directory + f for f in os.listdir(gp_directory) if f.endswith('.bam')]
        try: assert(len(bam_files) == 1)
        except:
            logObject.error()
            raise RuntimeError

        input_bam = bam_files[0]

    # create Alignment Object
    AlignmentObj = Alignment(input_bam, sample_name, logObject)

    # extract reads from BAM
    AlignmentObj.extract_fastqs(workspace, compress=(not no_gzip))

    # create successful completion file if steps completed!
    conf_file = open(parent_dir + "GPPROCESS.txt", 'w')
    conf_file.write("ProcessGPDirectory: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Extracts FASTQ files from BAM file and gather picard metrics.
    """)

    parser.add_argument('-i', '--bam_location', help="GP directory for sample.", required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-z', '--no_gzip', action='store_true', help='leave resulting files uncompressed. Default is False (compresses results).', required=False, default=False)
    args = parser.parse_args()

    processGpDirectory(args.bam_location, args.sample, args.sample_dir, args.no_gzip)
