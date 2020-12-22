import os
import sys

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

        run_subsample = False
        read_subsampling = 100000
        trimgalore_options = ""
        run_centrifuge = False
        centrifuge_memory= 32
        centrifuge_threads = 1
        centrifuge_index = None # recommended to set this when installing seQuoia
        run_kneaddata = True
        kneaddata_options = "" # recommended to set this when installing seQuoia (for instance to set path to human reference genome)
        kneaddata_threads = 1
        kneaddata_memory = 32
        run_metaphlan2 = True
        metaphlan2_threads = 1
        run_sortmerna = False
        sortmerna_db = "" # recommended to set this when installing seQuoia
        sortmerna_threads = 4
        sortmerna_memory = 16
        run_straingst = False
        straingst_db = None # recommended to set this when installing seQuoia
        run_shortbred = False
        amr_shortbred_markers = None # recommended to set this when installing seQuoia
        set_timelimit_sevendays = False
        run_cleanup = True

        arguments = {'run_subsample': run_subsample,
                     'read_subsampling': read_subsampling,
                     'trimgalore_options': trimgalore_options,
                     'run_centrifuge': run_centrifuge,
                     'centrifuge_memory': centrifuge_memory,
                     'centrifuge_threads': centrifuge_threads,
                     'centrifuge_index': centrifuge_index,
                     'run_kneaddata': run_kneaddata,
                     'kneaddata_threads': kneaddata_threads,
                     'kneaddata_memory': kneaddata_memory,
                     'kneaddata_options': kneaddata_options,
                     'run_metaphlan2': run_metaphlan2,
                     'metaphlan2_threads': metaphlan2_threads,
                     'run_sortmerna': run_sortmerna,
                     'sortmerna_threads': sortmerna_threads,
                     'sortmerna_memory': sortmerna_memory,
                     'sortmerna_db': sortmerna_db,
                     'run_straingst': run_straingst,
                     'straingst_db': straingst_db,
                     'run_shortbred': run_shortbred,
                     'amr_shortbred_markers': amr_shortbred_markers,
                     'set_timelimit_sevendays': set_timelimit_sevendays,
                     'run_cleanup': run_cleanup}

        # overwrite variable values with config provided values
        if self.config_parameters:
            for var in self.config_parameters.items():
                if len(var[1].split()) > 1: arguments[var[0]] = uF.wrapAppropriately(var[1])
                else: arguments[var[0]] = var[1]

        boolean_arguments = ['run_subsample', 'run_centrifuge', 'run_kneaddata', 'run_centrifuge', 'run_sortmerna', 'run_shortbred',
                             'run_straingst', 'run_cleanup', 'run_metaphlan2', 'set_timelimit_sevendays']
        for ba in boolean_arguments:
            if str(arguments[ba]).lower() in uF.true_set:
                arguments[ba] = True
            else:
                arguments[ba] = False

        ####
        ## Start the pipeline!
        ####

        # Step 1: Run Setup_task module
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
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(gp_process_cmd, "gp_process", self.sample, self.grid_dir, cluster=self.cluster, timelimit="7:00:00:00", project=self.project)
                else:
                    uF.runSingleJob(gp_process_cmd, "gp_process", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)

                uF.check_completion(gp_process_checkpoint_file)

            if not (self.singleFastq or (self.forwardFastq and self.reverseFastq)):
                self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("gp_process",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)

        # Run Input Symlinking (Dumb module to just symlink illumina sequencing reads to have sample names in files)
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

        if arguments['run_subsample']:
            # Run Subsample
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
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(subsample_cmd, "subsample", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit='7:00:00:00',  project=self.project)
                else:
                    uF.runSingleJob(subsample_cmd, "subsample", self.sample, self.grid_dir, cluster=self.cluster, mem=48, project=self.project)
                uF.check_completion(subsample_checkpoint_file)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("subsample",
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
                                                                         "--sample_dir", self.sample_workspace,
                                                                         "--trimgalore_options=" + arguments['trimgalore_options']]
            if arguments['set_timelimit_sevendays']:
                uF.runSingleJob(cutadapt_cmd, "adapter_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit="7:00:00:00", threads=1, project=self.project)
            else:
                uF.runSingleJob(cutadapt_cmd, "adapter_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit="12:00:00", threads=1, project=self.project)
            uF.check_completion(adaptertrim_checkpoint_file)

        self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("adapter_trim",
                                                                                           self.sample_workspace,
                                                                                           self.singleFastq,
                                                                                           self.forwardFastq,
                                                                                           self.reverseFastq)


        if arguments['run_centrifuge'] and arguments['centrifuge_index']:
            # Run Centrifuge
            ###################################################################################################################################
            centrifuge_checkpoint_file = self.sample_workspace + "CENTRIFUGE.txt"
            if not os.path.isfile(centrifuge_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Centrifuge/')
                centrifuge_cmd = ["sh", self.code_dir + "Centrifuge_task.sh", "--fastq_sin=" + self.singleFastq,
                                  "--fastq_frw=" + self.forwardFastq,
                                  "--fastq_rev=" + self.reverseFastq,
                                  "--sample=" + self.sample,
                                  "--sample_dir", self.sample_workspace,
                                  "--centrifuge_index=" + arguments['centrifuge_index'],
                                  "--cores", str(arguments['centrifuge_threads'])]

                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(centrifuge_cmd, "centrifuge", self.sample, self.grid_dir, cluster=self.cluster, mem=int(arguments['centrifuge_memory']), threads=int(arguments['centrifuge_threads']), timelimit="7:00:00:00", project=self.project)
                else:
                    uF.runSingleJob(centrifuge_cmd, "centrifuge", self.sample, self.grid_dir, cluster=self.cluster, mem=int(arguments['centrifuge_memory']), threads=int(arguments['centrifuge_threads']), timelimit="12:00:00", project=self.project)

                uF.check_completion(centrifuge_checkpoint_file)

        if arguments['run_sortmerna']:
            # Run SortMeRNA
            ###################################################################################################################################
            sortmerna_checkpoint_file = self.sample_workspace + "SORTMERNA.txt"
            if not os.path.isfile(sortmerna_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'SortMeRNA/')
                sortmerna_cmd = ["sh", self.code_dir + "SortMeRNA_task.sh", "--fastq_sin=" + self.singleFastq,
                                 "--fastq_frw=" + self.forwardFastq,
                                 "--fastq_rev=" + self.reverseFastq,
                                 "--sample=" + self.sample,
                                 "--sample_dir=" + self.sample_workspace,
                                 "--database_dir=" + arguments['sortmerna_db'],
                                 "--cores", str(arguments['sortmerna_threads'])]
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(sortmerna_cmd, "sortmerna", self.sample, self.grid_dir, cluster=self.cluster, project=self.project, timelimit="7:00:00:00", mem=int(arguments['sortmerna_memory']), threads=int(arguments['sortmerna_threads']))
                else:
                    uF.runSingleJob(sortmerna_cmd, "sortmerna", self.sample, self.grid_dir, cluster=self.cluster, project=self.project, timelimit="12:00:00", mem=int(arguments['sortmerna_memory']), threads=int(arguments['sortmerna_threads']))
                try:
                    assert (os.path.isfile(sortmerna_checkpoint_file))
                except:
                    raise RuntimeWarning(); sys.exit(1)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("sortmerna",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)

        if arguments['run_kneaddata']:
            # Run KneadData
            ###################################################################################################################################
            kneaddata_checkpoint_file = self.sample_workspace + "KNEADDATA.txt"
            if not os.path.isfile(kneaddata_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'KneadData/')
                kneaddata_cmd = ["sh", self.code_dir + "Kneaddata_task.sh", "--fastq_sin=" + self.singleFastq,
                                 "--fastq_frw=" + self.forwardFastq,
                                 "--fastq_rev=" + self.reverseFastq,
                                 "--sample=" + self.sample,
                                 "--sample_dir", self.sample_workspace,
                                 "--kneaddata_options=" + arguments['kneaddata_options'],
                                 "--cores", str(arguments['kneaddata_threads'])]
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(kneaddata_cmd, "kneaddata", self.sample, self.grid_dir, cluster=self.cluster, project=self.project, timelimit="7:00:00:00", mem=int(arguments['kneaddata_memory']), threads=int(arguments['kneaddata_threads']))
                else:
                    uF.runSingleJob(kneaddata_cmd, "kneaddata", self.sample, self.grid_dir, cluster=self.cluster, project=self.project, timelimit="12:00:00", mem=int(arguments['kneaddata_memory']), threads=int(arguments['kneaddata_threads']))

                try:
                    assert (os.path.isfile(kneaddata_checkpoint_file))
                except:
                    raise RuntimeWarning();
                    sys.exit(1)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("kneaddata",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)

        # Run FastQC
        ###################################################################################################################################
        fastqc_checkpoint_file = self.sample_workspace + "FASTQC.txt"
        if not os.path.isfile(fastqc_checkpoint_file):
            uF.cleanUp(self.sample_workspace + 'FastQC/')
            fastqc_cmd = ["sh", self.code_dir + "FastQC_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                  "--fastq_frw=" + self.forwardFastq,
                                                                  "--fastq_rev=" + self.reverseFastq,
                                                                  "--sample=" + self.sample,
                                                                  "--sample_dir", self.sample_workspace,
                                                                  "--cores", "1"]
            if arguments['set_timelimit_sevendays']:
                uF.runSingleJob(fastqc_cmd, "fastqc", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit="7:00:00:00", project=self.project)
            else:
                uF.runSingleJob(fastqc_cmd, "fastqc", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit="12:00:00", project=self.project)
            uF.check_completion(fastqc_checkpoint_file)

        if arguments['run_straingst'] and arguments['straingst_db']:
            # Run StrainGST
            ###################################################################################################################################
            straingst_checkpoint_file = self.sample_workspace + "STRAINGST.txt"
            if not os.path.isfile(straingst_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'StrainGST/')
                straingst_cmd = ["sh", self.code_dir + "StrainGST_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                              "--fastq_frw=" + self.forwardFastq,
                                                                              "--fastq_rev=" + self.reverseFastq,
                                                                              "--sample=" + self.sample,
                                                                              "--sample_dir", self.sample_workspace,
                                                                              "--db=" + arguments['straingst_db']]
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(straingst_cmd, "straingst", self.sample, self.grid_dir, cluster=self.cluster, timelimit='7:00:00:00', mem=20, threads=1, project=self.project)
                else:
                    uF.runSingleJob(straingst_cmd, "straingst", self.sample, self.grid_dir, cluster=self.cluster, timelimit='12:00:00', mem=20, threads=1, project=self.project)
                uF.check_completion(straingst_checkpoint_file)

        if arguments['run_metaphlan2']:
            # Run Metaphlan Taxonomic Profiling
            ###################################################################################################################################
            metaphlan_checkpioing_file = self.sample_workspace + "METAPHLAN.txt"
            if not os.path.isfile(metaphlan_checkpioing_file):
                uF.cleanUp(self.sample_workspace + 'MetaPhlAn2/')
                refalign_cmd = ["sh", self.code_dir + "Metaphlan_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                    "--fastq_frw=" + self.forwardFastq,
                                                                    "--fastq_rev=" + self.reverseFastq,
                                                                    "--cores=" + str(arguments['metaphlan2_threads']),
                                                                    "--sample=" + self.sample,
                                                                    "--sample_dir", self.sample_workspace]
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(refalign_cmd, "metaphlan", self.sample, self.grid_dir, cluster=self.cluster, mem=32, timelimit="7:00:00:00", threads=arguments['metaphlan2_threads'], project=self.project)
                else:
                    uF.runSingleJob(refalign_cmd, "metaphlan", self.sample, self.grid_dir, cluster=self.cluster, mem=32, timelimit="12:00:00", threads=arguments['metaphlan2_threads'], project=self.project)

                uF.check_completion(metaphlan_checkpioing_file)

        if arguments['run_shortbred'] and arguments['amr_shortbred_markers']:
            # Run AMR Identification
            ###################################################################################################################################
            amrp_checkpoint_file = self.sample_workspace + "AMRP.txt"
            if not os.path.isfile(amrp_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'AMRP/')
                refalign_cmd = ["sh", self.code_dir + "AMRP_task.sh", "--fastq_frw=" + self.forwardFastq,
                                                                    "--fastq_rev=" + self.reverseFastq,
                                                                    "--sample=" + self.sample,
                                                                    "--sample_dir=" + self.sample_workspace,
                                                                    "--shortbred_markers", arguments['amr_shortbred_markers']]
                if arguments['set_timelimit_sevendays']:
                    uF.runSingleJob(refalign_cmd, "amr_prediction", self.sample, self.grid_dir, cluster=self.cluster, timelimit='7:00:00:00', mem=48, threads=1, project=self.project)
                else:
                    uF.runSingleJob(refalign_cmd, "amr_prediction", self.sample, self.grid_dir, cluster=self.cluster, mem=48, timelimit='12:00:00', threads=1, project=self.project)

                uF.check_completion(amrp_checkpoint_file)

        if arguments['run_cleanup']:
            # Step 5: Run Cleanup
            ###################################################################################################################################
            cleanup_checkpoint_file = self.sample_workspace + "CLEANUP.txt"
            if not os.path.isfile(cleanup_checkpoint_file):
                cleanup_cmd = ["sh", self.code_dir + "Cleanup_task.sh", "--sample_dir", self.sample_workspace]
                uF.runSingleJob(cleanup_cmd, "cleanup", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                uF.check_completion(cleanup_checkpoint_file)