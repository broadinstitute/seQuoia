import os
import sys
import argparse
from seQuoia.classes.FastqAnalyzer import Fastq, FastqPaired
from seQuoia.classes.AlignmentAnalyzer import Alignment
from seQuoia.classes.NanoporeAnalyzer import Nanopore
from seQuoia.other import usefulFunctions as uF
from seQuoia.other import nanoporeFunctions as npF

def runPilon(fastq_frw, fastq_rev, nanopore_fastq, assembly, sample_name, sample_dir, pilon_options, max_iterations, identifier, cores):

    # set up directory structure
    workspace_name = "Pilon_Polishing"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Pilon_Polishing.log'
    logObject = uF.createLoggerObject(log_file)

    best_assembly = workspace + sample_name + '.pilon-refined.fasta'

    IlluminaObj = None
    NanoporeObj = None
    try:
        assert( (fastq_frw and fastq_rev and os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev) ) or (nanopore_fastq and os.path.isfile(nanopore_fastq)))
        if os.path.isfile(fastq_frw) and os.path.isfile(fastq_rev):
            IlluminaObj = FastqPaired(fastq_frw, fastq_rev, sample_name, logObject)
            IlluminaObj.validate()
    except:
        sys.stderr.write("ERROR: Sequencing data inputs are missing or not properly provided. Exiting now ...\n")
        raise RuntimeError

    try:
        assert( os.path.isfile(assembly))
        if os.path.isfile(nanopore_fastq):
            NanoporeObj = Nanopore(nanopore_fastq, sample_name + '_ont', logObject)
            NanoporeObj.validate()
    except:
        sys.stderr.write("ERROR: Input assembly to be refined has a faulty path. Exiting now ...\n")
        raise RuntimeError

    logObject.info("Beginning pilon polishing process!")

    iteration = 1
    while iteration < (max_iterations+1):

        logObject.info("<<< starting iteration %d >>>" % iteration)

        iteration_workspace_name = 'iteration_' + str(iteration)
        iteration_workspace = uF.setupDirectory(workspace, iteration_workspace_name, panic_if_exists=False)

        iteration_log_file = workspace + 'iteration_' + str(iteration) + '.log'
        iterationLogObject = uF.createLoggerObject(iteration_log_file)

        IlluminaAlignmentObj = None
        NanoporeAlignmentObj = None
        # Align Illumina Reads to Assembly
        if IlluminaObj:

            # Reinitialize for logging
            iterationIlluminaObj = FastqPaired(fastq_frw, fastq_rev, sample_name, iterationLogObject)

            # Align reads to reference genome
            illumina_sam_file = iterationIlluminaObj.align_to_reference(iteration_workspace, assembly, cores=cores)

            # create Alignment object
            IlluminaAlignmentObj = Alignment(illumina_sam_file, sample_name, iterationLogObject)

            # compress SAM to BAM
            IlluminaAlignmentObj.compress_sam(iteration_workspace, clean=True)

            # sort BAM file
            IlluminaAlignmentObj.sort_bam(iteration_workspace, cores=cores, clean=True)

            # mark duplicates
            IlluminaAlignmentObj.mark_dups(iteration_workspace, clean=True)

            # index BAM file
            IlluminaAlignmentObj.index_bam(iteration_workspace)

        if NanoporeObj:
            # Reinitialize for logging
            iterationNanoporeObj = Nanopore(nanopore_fastq, sample_name + '_ont', iterationLogObject)

            # Align nanopore reads to reference genome
            nanopore_sam_file = iterationNanoporeObj.minimap_alignment(assembly, iteration_workspace, cores=cores)

            # create Alignment Object
            NanoporeAlignmentObj = Alignment(nanopore_sam_file, sample_name + '_ont', iterationLogObject)

            # compress SAM to BAM
            NanoporeAlignmentObj.compress_sam(iteration_workspace, clean=True)

            # sort BAM file
            NanoporeAlignmentObj.sort_bam(iteration_workspace, clean=True, cores=cores)

            # index BAM file
            NanoporeAlignmentObj.index_bam(iteration_workspace)

        # Run Pilon Assembly Refinement
        polished_assembly = None
        if IlluminaObj and NanoporeObj:
            polished_assembly, changes = npF.run_pilon(assembly, sample_name, iteration_workspace, iterationLogObject,
                                              illumina_bam=IlluminaAlignmentObj.alignment_file,
                                              nanopore_bam=NanoporeAlignmentObj.alignment_file, options=pilon_options, cores=cores)
        elif NanoporeObj:
            polished_assembly, changes = npF.run_pilon(assembly, sample_name, iteration_workspace, iterationLogObject,
                                              nanopore_bam=NanoporeAlignmentObj.alignment_file, options=pilon_options, cores=cores)
        elif IlluminaObj:
            polished_assembly, changes = npF.run_pilon(assembly, sample_name, iteration_workspace, iterationLogObject,
                                              illumina_bam=IlluminaAlignmentObj.alignment_file, options=pilon_options, cores=cores)

        os.system('mv %s %s' % (iteration_workspace + 'pilon.log', workspace + 'iteration_' + str(iteration) + '.pilon.log'))
        os.system('mv %s %s' % (iteration_workspace + 'pilon.err', workspace + 'iteration_' + str(iteration) + '.pilon.err'))

        convergence = uF.compare_assemblies(assembly, polished_assembly)

        if polished_assembly:
            if os.path.isfile(best_assembly): os.system('rm -f %s*' % best_assembly)
            os.system('mv %s %s' % (polished_assembly, best_assembly))
            os.system('mv %s %s' % (changes, workspace + 'pilon_iteration_' + str(iteration) + '.changes'))

        assembly = best_assembly
        os.system('rm -rf %s' % iteration_workspace)

        if convergence:
            logObject.info("iteration %d finished with convergence." % iteration)
            break
        else:
            logObject.info("iteration %d finished without converging." % iteration)

        uF.closeLoggerObject(iterationLogObject)
        iteration += 1

    if max_iterations == iteration:
        logObject.info("No convergence could be reached after %d iterations" % max_iterations)
        logObject.info("Using last assembly as final polished version.")

    conf_file_name = sample_dir + "PILON_POLISHING"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("Pilon Assembly Polishing: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""Run Pilon Polishing for Assemblies. needs either paired end illumina data or nanopore fastq data.""")

    parser.add_argument('-1', '--fastq_frw', help='provide path to forward illumina FASTQ file.', required=False, default=None)
    parser.add_argument('-2', '--fastq_rev', help='provide path to reverse illumina FASTQ file.', required=False, default=None)
    parser.add_argument('-n', '--nanopore_fastq', help='provide path to nanopore FASTQ file.', required=False, default=None)
    parser.add_argument('-a', '--assembly', help='provide path to unrefined assembly.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-p', '--pilon_options', help='Pilon options.', required=False, default="")
    parser.add_argument('-m', '--max_iterations', type=int, help='maximum iterations to do before calling it quits. default is 10.', required=False, default=10)
    parser.add_argument('-i', '--identifier', help='identifier, in case multiple runs in workflow.', required=False, default=None)
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    args = parser.parse_args()

    runPilon(args.fastq_frw, args.fastq_rev, args.nanopore_fastq, args.assembly, args.sample, args.sample_dir, args.pilon_options, args.max_iterations, args.identifier, args.cores)