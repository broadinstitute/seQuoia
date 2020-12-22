import os
import sys

true_set = set(['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh'])
class Workflow():
    def __init__(self, sample, workspace, metadata, grid_dir, code_dir,
                 singleFastq="", forwardFastq="", reverseFastq="", bamLocation="",
                 nanoporeFastq="", nanoporeFast5="", nanoporeSeqSummary="", nanoporeBarcode="",
                 cluster="Local", project=None, config_parameters=None):
        self.workspace = os.path.abspath(workspace) + '/'
        self.sample = sample
        self.bamLocation = bamLocation
        self.grid_dir = grid_dir
        self.singleFastq = singleFastq
        self.forwardFastq = forwardFastq
        self.reverseFastq = reverseFastq
        self.nanoporeFastq = nanoporeFastq
        self.nanoporeFast5 = nanoporeFast5
        self.nanoporeSeqSummary = nanoporeSeqSummary
        self.nanoporeBarcode = nanoporeBarcode
        self.metadata = metadata
        self.cluster = cluster
        self.project = project
        self.config_parameters = config_parameters
        self.seqc_dir = code_dir
        self.code_dir = code_dir + '/tasks/'
        self.sample_workspace = self.workspace + sample + '/'
        if not os.path.isfile(self.metadata):
            self.metadata = "'\"" + self.metadata + "\"'"


    """ 
    This is what the USER should alter as they see fit. E.g. the options for running external tools or the arrangment 
    of modules / flow of the pipeline. 
    """
    def runWorkflow(self):
        from seQuoia.other import usefulFunctions as uF

        ##################################
        ### Open for User Modification ###
        ##################################

        ####
        ## Start the pipeline!
        ####

        """ Parameters """
        reference_fasta = None
        pilon_options = "--variant"
        run_cleanup = False

        arguments = {'reference_fasta': reference_fasta,
                     'pilon_options': pilon_options,
                     'run_cleanup': run_cleanup}

        # overwrite variable values with config provided values
        if self.config_parameters:
            for var in self.config_parameters.items():
                if len(var[1].split()) > 1: arguments[var[0]] = uF.wrapAppropriately(var[1])
                else: arguments[var[0]] = var[1]

        if arguments['run_cleanup'].lower() in uF.true_set:
            arguments['run_cleanup'] = True

        # Run Setup_task module
        ###################################################################################################################################
        setup_checkpoint_file = self.sample_workspace + "SETUP.txt"
        if not os.path.isfile(setup_checkpoint_file):
            setup_cmd = ["sh", self.code_dir + "Setup_task.sh", "--sample", self.sample,
                                                                "--parent_dir", self.workspace,
                                                                "--meta", self.metadata]
            uF.runSingleJob(setup_cmd, "setup", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
            uF.check_completion(setup_checkpoint_file)

        if self.bamLocation:
            # Convert GP BAM file to FASTQ(s) and store any available Picard metric reports
            ###################################################################################################################################
            gp_process_checkpoint_file = self.sample_workspace + "GPPROCESS.txt"
            if not os.path.isfile(gp_process_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'ProcessGPDirectory/')
                gp_process_cmd = ["sh", self.code_dir + "ProcessGPDirectory_task.sh", "--bam_location=" + self.bamLocation,
                                                                                      "--sample=" + self.sample,
                                                                                      "--sample_dir=" + self.sample_workspace]
                uF.runSingleJob(gp_process_cmd, "gp_process", self.sample, self.grid_dir, mem=32, cluster=self.cluster, project=self.project)
                uF.check_completion(gp_process_checkpoint_file)

            if not (self.singleFastq or (self.forwardFastq and self.reverseFastq)):
                self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("gp_process",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)


        # Step 1.75: Run Input Symlinking (Dumb module to just symlink illumina sequencing reads to have sample names in files)
        ###################################################################################################################################
        symlink_checkpoint_file = self.sample_workspace + "SYMLINK_INPUT.txt"
        if not os.path.isfile(symlink_checkpoint_file):
            uF.cleanUp(self.sample_workspace + 'Symlink_Input/')
            symlink_cmd = ["sh", self.code_dir + "SymlinkInput_task.sh", "--fastq_sin=" + self.singleFastq,
                          "--fastq_frw=" + self.forwardFastq,
                          "--fastq_rev=" + self.reverseFastq,
                          "--sample=" + self.sample,
                          "--sample_dir=" + self.sample_workspace]
            uF.runSingleJob(symlink_cmd, "symlink_input", self.sample, self.grid_dir, cluster=self.cluster,
                            project=self.project)
            uF.check_completion(symlink_checkpoint_file)

        self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("symlink",
                                                                                           self.sample_workspace,
                                                                                           self.singleFastq,
                                                                                           self.forwardFastq,
                                                                                           self.reverseFastq)


        # Run Adapter Trimming
        ###################################################################################################################################
        adaptertrim_checkpoint_file = self.sample_workspace + "ADAPTERTRIM.txt"
        if not os.path.isfile(adaptertrim_checkpoint_file):
            uF.cleanUp(self.sample_workspace + 'AdapterTrim/')
            cutadapt_cmd = ["sh", self.code_dir + "AdapterTrim_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                         "--fastq_frw=" + self.forwardFastq,
                                                                         "--fastq_rev=" + self.reverseFastq,
                                                                         "--sample=" + self.sample,
                                                                         "--sample_dir", self.sample_workspace]
            uF.runSingleJob(cutadapt_cmd, "adapter_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=10,
                            threads=1, project=self.project)
            uF.check_completion(adaptertrim_checkpoint_file)

        self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("adapter_trim",
                                                                                           self.sample_workspace,
                                                                                           self.singleFastq,
                                                                                           self.forwardFastq,
                                                                                           self.reverseFastq)

        if arguments['reference_fasta']:
            # Step 7: Run Reference Alignment (any BWA reference index should already be built,
            #         this is not the objective of seQc), includes post processing steps
            ###################################################################################################################################
            reference_align_checkpoint_file = self.sample_workspace + "REFALIGNMENT.txt"
            if not os.path.isfile(reference_align_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'ReferenceAlignment/')
                refalign_cmd = ["sh", self.code_dir + "RefAlignment_task.sh",
                                "--fastq_sin=" + self.singleFastq,
                                "--fastq_frw=" + self.forwardFastq,
                                "--fastq_rev=" + self.reverseFastq,
                                "--sample=" + self.sample,
                                "--sample_dir=" + self.sample_workspace,
                                "--reference_fasta", arguments['reference_fasta']]
                uF.runSingleJob(refalign_cmd, "ref_alignment", self.sample, self.grid_dir, cluster=self.cluster, mem=40,
                                threads=1, project=self.project)
                try:
                    assert (os.path.isfile(reference_align_checkpoint_file))
                except:
                    raise RuntimeWarning(); sys.exit(1)
            bam_file = uF.extractResultingAlignment("reference_aligner", self.sample_workspace)

            # Step 8: Run Pilon
            ###################################################################################################################################
            pilon_checkpoint_file = self.sample_workspace + "PILON.txt"
            if not os.path.isfile(pilon_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Pilon/')
                pilon_cmd = ["sh", self.code_dir + "Pilon_task.sh",
                             "--bam_file=" + bam_file,
                             "--sample="+  self.sample,
                             "--sample_dir=" + self.sample_workspace,
                             "--pilon_options=" + arguments['pilon_options'],
                             "--reference_fasta=" + os.path.abspath(arguments['reference_fasta'])]
                uF.runSingleJob(pilon_cmd, "pilon", self.sample, self.grid_dir, cluster=self.cluster, mem=40, threads=1,
                                project=self.project)
                try:
                    assert (os.path.isfile(pilon_checkpoint_file))
                except:
                    raise RuntimeWarning(); sys.exit(1)


        else:
            sys.stderr.write("You forgot the reference! Try again ...\n")
            raise RuntimeWarning
            sys.exit(1)

        if arguments['run_cleanup']:
            # Run Cleanup
            ###################################################################################################################################
            cleanup_checkpoint_file = self.sample_workspace + "CLEANUP.txt"
            if not os.path.isfile(cleanup_checkpoint_file):
                cleanup_cmd = ["sh", self.code_dir + "Cleanup_task.sh", "--sample_dir", self.sample_workspace]
                uF.runSingleJob(cleanup_cmd, "cleanup", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                uF.check_completion(cleanup_checkpoint_file)