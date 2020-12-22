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
        run_subsample = False
        read_subsampling = 1000000
        gaemr_threads = 4
        gaemr_memory = 24
        gaemr_timelimit = "48:00:00"
        gaemr_formatter_options = uF.wrapAppropriately("-g 1 -c 100 -r")
        gaemr_qc_options = uF.wrapAppropriately("--force --analyze_rna")
        run_cleanup = False
        arguments = {'reference_fasta': reference_fasta,
                     'run_subsample': run_subsample,
                     'read_subsampling': read_subsampling,
                     'gaemr_formatter_options': gaemr_formatter_options,
                     'gaemr_qc_options': gaemr_qc_options,
                     'gaemr_threads': gaemr_threads,
                     'gaemr_memory': gaemr_memory,
                     'gaemr_timelimit': gaemr_timelimit,
                     'run_cleanup': run_cleanup}

        # overwrite variable values with config provided values
        if self.config_parameters:
            for var in self.config_parameters.items():
                if len(var[1].split()) > 1: arguments[var[0]] = uF.wrapAppropriately(var[1])
                else: arguments[var[0]] = var[1]

        boolean_arguments = ['run_subsample', 'run_cleanup']
        for ba in boolean_arguments:
            if str(arguments[ba]).lower() in uF.true_set:
                arguments[ba] = True
            else:
                arguments[ba] = False

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
            # Step 1.5: Convert GP BAM file to FASTQ(s) and store any available Picard metric reports
            ###################################################################################################################################
            gp_process_checkpoint_file = self.sample_workspace + "GPPROCESS.txt"
            if not os.path.isfile(gp_process_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'ProcessGPDirectory/')
                gp_process_cmd = ["sh", self.code_dir + "ProcessGPDirectory_task.sh", "--bam_location=" + self.bamLocation,
                                                                                      "--sample=" + self.sample,
                                                                                      "--sample_dir=" + self.sample_workspace]
                uF.runSingleJob(gp_process_cmd, "gp_process", self.sample, self.grid_dir, mem=32, timelimit="10:00:00", cluster=self.cluster, project=self.project)
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
                          "--sample_dir=" + self.sample_workspace,
                          "--sample=" + self.sample]
            uF.runSingleJob(symlink_cmd, "symlink_input", self.sample, self.grid_dir, cluster=self.cluster,
                            project=self.project)
            uF.check_completion(symlink_checkpoint_file)

        self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("symlink",
                                                                                           self.sample_workspace,
                                                                                           self.singleFastq,
                                                                                           self.forwardFastq,
                                                                                           self.reverseFastq)

        if arguments['run_subsample']:
            # Step 2: Run Subsample
            ###################################################################################################################################
            subsample_checkpoint_file = self.sample_workspace + "SUBSAMPLE.txt"
            if not os.path.isfile(subsample_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Subsample/')
                subsample_cmd = ["sh", self.code_dir + "Subsample_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                         "--fastq_frw=" + self.forwardFastq,
                                                                         "--fastq_rev=" + self.reverseFastq,
                                                                         "--sample=" + self.sample,
                                                                         "--sample_dir", self.sample_workspace,
                                                                         "--reads", str(arguments['read_subsampling']),
                                                                         "--no_gzip"]
                uF.runSingleJob(subsample_cmd, "subsample", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                uF.check_completion(subsample_checkpoint_file)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("subsample", self.sample_workspace,
                                                                                                self.singleFastq,
                                                                                                self.forwardFastq, self.reverseFastq)

        if arguments['reference_fasta']:
            # Run GAEMR Assembly QC
            ###################################################################################################################################
            gaemr_checkpoint_file = self.sample_workspace + "GAEMR.txt"
            if not os.path.isfile(gaemr_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'GAEMR/')

                gaemr_cmd = ["sh", self.code_dir + "GAEMR_task.sh",
                             "--assembly", arguments['reference_fasta'],
                             "--sample_dir", self.sample_workspace,
                             "--sample", self.sample,
                             "--format_options=" + str(arguments['gaemr_formatter_options']),
                             "--qc_options=" + str(arguments['gaemr_qc_options']),
                             "--illumina_frw=" + self.forwardFastq,
                             "--illumina_rev=" + self.reverseFastq,
                            "--cores", str(arguments['gaemr_threads'])]

                uF.runSingleJob(gaemr_cmd, "gaemr", self.sample, self.grid_dir, cluster=self.cluster,
                                timelimit=arguments['gaemr_timelimit'],
                                threads=int(arguments['gaemr_threads']), mem=int(arguments['gaemr_memory']),
                                project=self.project)
                uF.check_completion(gaemr_checkpoint_file)
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
                uF.runSingleJob(cleanup_cmd, "cleanup", self.sample, self.grid_dir, cluster=self.cluster,
                                project=self.project)
                uF.check_completion(cleanup_checkpoint_file)