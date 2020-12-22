import os
import argparse
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.classes.AlignmentAnalyzer import Alignment
from seQuoia.classes.AssemblyAnalyzer import Assembly
from seQuoia.other import usefulFunctions as uF

def runNanoPolish(assembly, nanopore_fastq, nanopore_fast5, nanopore_barcode, sample_dir, sample_name, nanopolish_variants_options, identifier, cores):
    try: assert(cores >= 2 and cores%2 == 0)
    except:
        raise RuntimeError("ERROR: Number of cores required for Nanopolish need to be larger than or equal to 2 and a multiple of 2.\n")
    try: assert( assembly and os.path.isfile(assembly) )
    except:
        raise RuntimeError("ERROR: Assembly was not provided properly. Please fix. Raising exception\n")

    try: assert( nanopore_fastq and os.path.isfile(nanopore_fastq) )
    except:
        raise RuntimeError("ERROR: FASTQ input(s) were not provided properly. Please fix. Raising exception\n")

    try: assert( nanopore_fast5 and os.path.isdir(nanopore_fast5))
    except:
        raise RuntimeError("ERROR: fast5 input(s) were not provided properly. Please fix. Raising exception\n")

    # set up directory structure
    workspace_name = "Nanopolish"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'NanoPolish_Assembly.log'
    logObject = uF.createLoggerObject(log_file)

    # Initialize Nanopore Object
    NanoporeObj = Nanopore(nanopore_fastq, sample_name, logObject, fast5=nanopore_fast5, barcode=nanopore_barcode)

    workspace_a = uF.setupDirectory(workspace, "indexing")

    # Index Nanopore FASTQ with Signal Data from Fast5
    NanoporeObj.index_fast5(workspace_a)

    workspace_b = uF.setupDirectory(workspace, "alignment")

    # Run minimap alignment of indexed reads to assembly
    sam_alignment = NanoporeObj.minimap_alignment(assembly, workspace_b, cores=cores)

    # Initialize Alignment Object
    NanoporeAlignmentObj = Alignment(sam_alignment, sample_name, logObject)

    # compress SAM to BAM
    NanoporeAlignmentObj.compress_sam(workspace_b, clean=True)

    # sort BAM file
    NanoporeAlignmentObj.sort_bam(workspace_b, clean=True)

    # index BAM file
    NanoporeAlignmentObj.index_bam(workspace_b)

    # Initialize Assembly Object
    AssemblyObj = Assembly(assembly, sample_name, logObject)

    workspace_c = uF.setupDirectory(workspace, "polishing")

    # Run Nanopolish
    AssemblyObj.run_nanopolish(NanoporeObj.fastq, NanoporeAlignmentObj.alignment_file,
                               workspace_c, nanopolish_variants_options=nanopolish_variants_options, cores=cores)

    # Remove / clean-up alignment file
    os.system('rm -f %s' % NanoporeAlignmentObj.alignment_file)

    # Write checkpoint file for successful completion of module!
    conf_file_name = sample_dir + "NANOPOLISH"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("Nanopolish: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""Module for performing assembly polishing with Nanopolish.""")

    parser.add_argument('-a', '--assembly', help='location of canu assembly', required=True)
    parser.add_argument('-q', '--nanopore_fastq', help='location of nanopore FASTQ file / directory with nanopore FASTQ files.', required=True)
    parser.add_argument('-5', '--nanopore_fast5', help='directory of nanopore FAST5 files.', required=True)
    parser.add_argument('-b', '--nanopore_barcode', help='barcode identifier for sample.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-p', '--nanopolish_variants_options', help='Provide options for Nanopolish variants --consensus step.', required=False, default="")
    parser.add_argument('-i', '--identifier', help='identifier, in case multiple runs in workflow.', required=False, default=None)
    parser.add_argument('-c', '--cores', type=int, help='number of cores', required=False, default=16)
    args = parser.parse_args()

    runNanoPolish(args.assembly, args.nanopore_fastq, args.nanopore_fast5, args.nanopore_barcode, args.sample_dir, args.sample, args.nanopolish_variants_options, args.identifier, args.cores)