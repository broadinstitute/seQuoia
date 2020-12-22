import os
import sys

true_set = set(['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh'])
class Workflow():
    def __init__(self, sample, workspace, metadata, grid_dir, code_dir,
                 singleFastq="", forwardFastq="", reverseFastq="", bamLocation="",
                 nanoporeFastq="", nanoporeFast5="", nanoporeSeqSummary="", nanoporeBarcode="",
                 cluster="UGER", project=None, config_parameters=None):
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
        sys.path.append(self.code_dir)
        from seQuoia.other import usefulFunctions as uF

        ##################################
        ### Open for User Modification ###
        ##################################

        ####
        ## Start the pipeline!
        ####

        """ Parameters """

        """ Default Taxa-Non-Specfic Database files """


        run_adaptertrim = True
        trimgalore_options = ""
        run_qualitytrim = False
        trimmomatic_options = uF.wrapAppropriately("LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15")
        run_store_input = True
        run_centrifuge = True
        centrifuge_index = None
        run_mlst_ariba = False
        mlst_ariba_db = None
        amr_ariba_card_db = None
        other_ariba_db_paths = ""
        other_ariba_db_names = ""
        run_straingst = True
        straingst_db = None
        run_pilon = False
        reference_fasta = None
        run_subsample_for_assembly = False
        read_subsampling = 1000000
        run_assembly = True
        spades_read_length = 150
        assembly_threads = 4
        assembly_memory = 16
        assembly_timelimit = "48:00:00"
        unicycler_flag = True
        run_gaemr = False
        gaemr_formatter_options = uF.wrapAppropriately("-g 1 -c 100 -r")
        gaemr_qc_options = uF.wrapAppropriately("--force --analyze_rna")
        run_cleanup = True

        arguments = {'run_adaptertrim': run_adaptertrim,
                     'trimgalore_options': trimgalore_options,
                     'run_qualitytrim': run_qualitytrim,
                     'trimmomatic_options': trimmomatic_options,
                     'run_store_input': run_store_input,
                     'run_centrifuge': run_centrifuge,
                     'centrifuge_index': centrifuge_index,
                     'run_mlst_ariba': run_mlst_ariba,
                     'mlst_ariba_db': mlst_ariba_db,
                     'amr_ariba_card_db': amr_ariba_card_db,
                     'other_ariba_db_paths': other_ariba_db_paths,
                     'other_ariba_db_names': other_ariba_db_names,
                     'run_straingst': run_straingst,
                     'straingst_db': straingst_db,
                     'run_pilon': run_pilon,
                     'reference_fasta': reference_fasta,
                     'run_subsample_for_assembly': run_subsample_for_assembly,
                     'read_subsampling': read_subsampling,
                     'run_assembly': run_assembly,
                     'unicycler_flag': unicycler_flag,
                     'spades_read_length': spades_read_length,
                     'assembly_threads': assembly_threads,
                     'assembly_memory': assembly_memory,
                     'assembly_timelimit': assembly_timelimit,
                     'run_gaemr': run_gaemr,
                     'gaemr_formatter_options': gaemr_formatter_options,
                     'gaemr_qc_options': gaemr_qc_options,
                     'run_cleanup': run_cleanup}

        # overwrite variable values with config provided values
        if self.config_parameters:
            for var in self.config_parameters.items():
                if len(var[1].split()) > 1 and not var[0] == "other_ariba_db_names" and not var[0] == "other_ariba_db_paths":
                    arguments[var[0]] = uF.wrapAppropriately(var[1])
                else: arguments[var[0]] = var[1]

        boolean_arguments = ['run_store_input', 'run_adaptertrim', 'run_qualitytrim', 'run_mlst_ariba', 'run_assembly',
                             'unicycler_flag', 'run_subsample_for_assembly', 'run_straingst', 'run_pilon', 'run_gaemr', 'run_cleanup']
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
            # Convert GP BAM file to FASTQ(s) and store any available Picard metric reports
            ###################################################################################################################################
            gp_process_checkpoint_file = self.sample_workspace + "GPPROCESS.txt"
            if not os.path.isfile(gp_process_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'ProcessGPDirectory/')
                gp_process_cmd = ["sh", self.code_dir + "ProcessGPDirectory_task.sh", "--bam_location=" + self.bamLocation,
                                                                                      "--sample=" + self.sample,
                                                                                      "--sample_dir=" + self.sample_workspace]
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


        # Run Adapter Trimming
        ###################################################################################################################################
        if arguments['run_adaptertrim']:
            adaptertrim_checkpoint_file = self.sample_workspace + "ADAPTERTRIM.txt"
            if not os.path.isfile(adaptertrim_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'AdapterTrim/')
                cutadapt_cmd = ["sh", self.code_dir + "AdapterTrim_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                             "--fastq_frw=" + self.forwardFastq,
                                                                             "--fastq_rev=" + self.reverseFastq,
                                                                             "--sample=" + self.sample,
                                                                             "--sample_dir=" + self.sample_workspace,
                                                                             "--trimgalore_options=" + arguments['trimgalore_options']]
                uF.runSingleJob(cutadapt_cmd, "adapter_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=5,
                                threads=1, project=self.project)
                uF.check_completion(adaptertrim_checkpoint_file)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("adapter_trim",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)


        # Run Quality Trimming
        ###################################################################################################################################
        if arguments['run_qualitytrim']:
            qualitytrim_checkpoint_file = self.sample_workspace + "QUALITYTRIM.txt"
            if not os.path.isfile(qualitytrim_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'QualityTrim/')
                trimmomatic_cmd = ["sh", self.code_dir + "QualityTrim_task.sh", "--fastq_sin=" + self.singleFastq,
                                   "--fastq_frw=" + self.forwardFastq,
                                   "--fastq_rev=" + self.reverseFastq,
                                   "--sample=" + self.sample,
                                   "--sample_dir=" + self.sample_workspace,
                                   "--trimmomatic_options=" + arguments['trimmomatic_options'],
                                   "--cores", "1"]
                uF.runSingleJob(trimmomatic_cmd, "quality_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=16,
                                threads=1, project=self.project)
                uF.check_completion(qualitytrim_checkpoint_file)

            self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("quality_trim",
                                                                                               self.sample_workspace,
                                                                                               self.singleFastq,
                                                                                               self.forwardFastq,
                                                                                               self.reverseFastq)

        # Store Input
        ###################################################################################################################################
        if arguments['run_store_input']:
            input_checkpoint_file = self.sample_workspace + "STOREINPUT.txt"
            if not os.path.isfile(input_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Input/')
                fastqc_cmd = ["sh", self.code_dir + "StoreInput_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                          "--fastq_frw=" + self.forwardFastq,
                                                                          "--fastq_rev=" + self.reverseFastq,
                                                                          "--sample=" + self.sample,
                                                                          "--sample_dir", self.sample_workspace]
                uF.runSingleJob(fastqc_cmd, "store_input", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                uF.check_completion(input_checkpoint_file)

        # Run FastQC
        ###################################################################################################################################
        fastqc_checkpoint_file = self.sample_workspace + "FASTQC.txt"
        if not os.path.isfile(fastqc_checkpoint_file):
            uF.cleanUp(self.sample_workspace + 'FastQC/')
            fastqc_cmd = ["sh", self.code_dir + "FastQC_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                  "--fastq_frw=" + self.forwardFastq,
                                                                  "--fastq_rev=" + self.reverseFastq,
                                                                  "--sample=" + self.sample,
                                                                  "--sample_dir=" + self.sample_workspace,
                                                                  "--cores", "1"]
            uF.runSingleJob(fastqc_cmd, "fastqc", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
            uF.check_completion(fastqc_checkpoint_file)

        # Run Centrifuge
        ###################################################################################################################################
        if arguments['run_centrifuge'] and arguments['centrifuge_index']:
            centrifuge_checkpoint_file = self.sample_workspace + "CENTRIFUGE.txt"
            if not os.path.isfile(centrifuge_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Centrifuge/')
                centrifuge_cmd = ["sh", self.code_dir + "Centrifuge_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                              "--fastq_frw=" + self.forwardFastq,
                                                                              "--fastq_rev=" + self.reverseFastq,
                                                                              "--sample=" + self.sample,
                                                                              "--sample_dir", self.sample_workspace,
                                                                              "--centrifuge_index=" + arguments['centrifuge_index'],
                                                                              "--cores", "1"]
                uF.runSingleJob(centrifuge_cmd, "centrifuge", self.sample, self.grid_dir, cluster=self.cluster, mem=32, threads=1, timelimit="05:00:00", project=self.project)
                uF.check_completion(centrifuge_checkpoint_file)

        # Run MLST
        ###################################################################################################################################
        if arguments['run_mlst_ariba'] and arguments['mlst_ariba_db']:
            mlst_checkpoint_file = self.sample_workspace + "MLST.txt"
            if not os.path.isfile(mlst_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'MLST/')
                mlst_cmd = ["sh", self.code_dir + "MLST_task.sh",  "--fastq_frw=" + self.forwardFastq,
                                                                   "--fastq_rev=" + self.reverseFastq,
                                                                   "--sample=" + self.sample,
                                                                   "--sample_dir=" + self.sample_workspace,
                                                                   "--ariba_db_dir=" + arguments['mlst_ariba_db']]
                uF.runSingleJob(mlst_cmd, "mlst", self.sample, self.grid_dir, cluster=self.cluster, mem=20, threads=1, project=self.project)
                uF.check_completion(mlst_checkpoint_file)

        # Run AMR Identification
        ####################################################################################################################################
        amrp_checkpoint_file = self.sample_workspace + "AMRP.txt"
        if not os.path.isfile(amrp_checkpoint_file):
            uF.cleanUp(self.sample_workspace + 'AMRP/')
            refalign_cmd = ["sh", self.code_dir + "AMRP_task.sh", "--fastq_frw=" + self.forwardFastq,
                                                                  "--fastq_rev=" + self.reverseFastq,
                                                                  "--sample=" + self.sample,
                                                                  "--sample_dir", self.sample_workspace,
                                                                  "--ariba_database", arguments['amr_ariba_card_db'], arguments['other_ariba_db_paths'],
                                                                  "--ariba_names", 'CARD', arguments['other_ariba_db_names']]
            uF.runSingleJob(refalign_cmd, "amr_prediction", self.sample, self.grid_dir, cluster=self.cluster, mem=20, threads=1, project=self.project)
            uF.check_completion(amrp_checkpoint_file)


        # Run StrainGST
        ###################################################################################################################################
        if arguments['run_straingst'] and arguments['straingst_db']:
            straingst_checkpoint_file = self.sample_workspace + "STRAINGST.txt"
            if not os.path.isfile(straingst_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'StrainGST/')
                straingst_cmd = ["sh", self.code_dir + "StrainGST_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                            "--fastq_frw=" + self.forwardFastq,
                                                                            "--fastq_rev=" + self.reverseFastq,
                                                                            "--sample=" + self.sample,
                                                                            "--sample_dir", self.sample_workspace,
                                                                            "--db=" + arguments['straingst_db']]
                uF.runSingleJob(straingst_cmd, "straingst", self.sample, self.grid_dir, cluster=self.cluster, mem=20, threads=1, project=self.project)
                uF.check_completion(straingst_checkpoint_file)

        if arguments['run_pilon']:
            # Run Reference Alignment (any BWA reference index should already be built, this is not the objective of seQc), includes post processing steps
            ###################################################################################################################################
            reference_align_checkpoint_file = self.sample_workspace + "REFALIGNMENT.txt"
            if not os.path.isfile(reference_align_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'ReferenceAlignment/')
                refalign_cmd = ["sh", self.code_dir + "RefAlignment_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                              "--fastq_frw=" + self.forwardFastq,
                                                                              "--fastq_rev=" + self.reverseFastq,
                                                                              "--sample=" + self.sample,
                                                                              "--sample_dir=" + self.sample_workspace,
                                                                              "--reference_fasta=" + arguments['reference_fasta']] # gathers the fasta file for the best reference identified by StrainGST
                uF.runSingleJob(refalign_cmd, "ref_alignment", self.sample, self.grid_dir, cluster=self.cluster, mem=20, threads=1, project=self.project)
                uF.check_completion(reference_align_checkpoint_file)

            bam_file = uF.extractResultingAlignment("reference_aligner", self.sample_workspace)

            # Run Pilon
            ###################################################################################################################################
            pilon_checkpoint_file = self.sample_workspace + "PILON.txt"
            if not os.path.isfile(pilon_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Pilon/')
                pilon_cmd = ["sh", self.code_dir + "Pilon_task.sh", "--bam_file=" + bam_file,
                                                                    "--sample=" + self.sample,
                                                                    "--sample_dir=" + self.sample_workspace,
                                                                    "--reference_fasta=" + arguments['reference_fasta']]  # gathers the fasta file for the best reference identified by StrainGST
                uF.runSingleJob(pilon_cmd, "pilon", self.sample, self.grid_dir, cluster=self.cluster, mem=24, threads=1, project=self.project)
                uF.check_completion(pilon_checkpoint_file)


        if arguments['run_assembly']:
            if arguments['run_subsample_for_assembly']:
                # Run Subsample
                ###################################################################################################################################
                subsample_checkpoint_file = self.sample_workspace + "SUBSAMPLE.txt"
                if not os.path.isfile(subsample_checkpoint_file):
                    uF.cleanUp(self.sample_workspace + 'Subsample/')
                    subsample_cmd = ["sh", self.code_dir + "Subsample_task.sh", "--fastq_sin=" + self.singleFastq,
                                     "--fastq_frw=" + self.forwardFastq,
                                     "--fastq_rev=" + self.reverseFastq,
                                     "--sample=" + self.sample,
                                     "--sample_dir=" + self.sample_workspace,
                                     "--reads", str(arguments['read_subsampling'])]
                    uF.runSingleJob(subsample_cmd, "subsample", self.sample, self.grid_dir, cluster=self.cluster,
                                    project=self.project)
                    uF.check_completion(subsample_checkpoint_file)

                self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("subsample",
                                                                                                   self.sample_workspace,
                                                                                                   self.singleFastq,
                                                                                                   self.forwardFastq,
                                                                                                   self.reverseFastq)

            # Run SPAdes Assembly
            ###################################################################################################################################
            assembly_checkpoint_file = self.sample_workspace + "ASSEMBLY.txt"
            if not os.path.isfile(assembly_checkpoint_file):
                uF.cleanUp(self.sample_workspace + 'Assembly/')
                assembly_cmd = ["sh", self.code_dir + "Assembly_task.sh", "--fastq_frw=" + self.forwardFastq,
                                                                      "--fastq_rev=" + self.reverseFastq,
                                                                      "--sample=" + self.sample,
                                                                      "--sample_dir", self.sample_workspace,
                                                                      "--read_length", str(arguments['spades_read_length']),
                                                                      "--cores", str(arguments['assembly_threads'])]
                if arguments['unicycler_flag']: assembly_cmd += ['--unicycler']
                uF.runSingleJob(assembly_cmd, "assembly", self.sample, self.grid_dir, cluster=self.cluster, mem=int(arguments['assembly_memory']), threads=int(arguments['assembly_threads']), timelimit=arguments['assembly_timelimit'], project=self.project)
                uF.check_completion(assembly_checkpoint_file)

            assembly = self.sample_workspace + '/Assembly/assembly.fasta'

            # Run GAEMR Assembly QC
            ###################################################################################################################################
            if arguments['run_gaemr']:
                gaemr_checkpoint_file = self.sample_workspace + "GAEMR.txt"
                if not os.path.isfile(gaemr_checkpoint_file):
                    uF.cleanUp(self.sample_workspace + 'GAEMR/')

                    gaemr_cmd = ["sh", self.code_dir + "GAEMR_task.sh",  "--assembly=" + assembly,
                                                                         "--sample_dir=" + self.sample_workspace,
                                                                         "--sample=" + self.sample,
                                                                         "--fastq_frw=" + self.forwardFastq,
                                                                         "--fastq_rev=" + self.reverseFastq,
                                                                         "--format_options=" + str(arguments['gaemr_formatter_options']),
                                                                         "--qc_options=" + str(arguments['gaemr_qc_options']),
                                                                         "--cores", '4']

                    uF.runSingleJob(gaemr_cmd, "gaemr", self.sample, self.grid_dir, cluster=self.cluster,
                                    timelimit="48:00:00",
                                    threads=4, mem=25,
                                    project=self.project)
                    uF.check_completion(gaemr_checkpoint_file)

        # Create LSARP Data Tables
        ###################################################################################################################################
        lsarp_checkpoint_file = self.sample_workspace + "LSARP.txt"
        if not os.path.isfile(lsarp_checkpoint_file):
            #uF.cleanUp(self.sample_workspace + 'LSARP_Data_Tables/')
            lsarp_cmd = ["sh", self.code_dir + "TablesForLSARP_task.sh", "--sample_dir", self.sample_workspace]
            uF.runSingleJob(lsarp_cmd, "lsarp_tables", self.sample, self.grid_dir, cluster=self.cluster, mem=5, threads=1, project=self.project)
            uF.check_completion(lsarp_checkpoint_file)



        # Run Cleanup
        ###################################################################################################################################
        if arguments["run_cleanup"]:
            cleanup_checkpoint_file = self.sample_workspace + "CLEANUP.txt"
            if not os.path.isfile(cleanup_checkpoint_file):
                cleanup_cmd = ["sh", self.code_dir + "Cleanup_task.sh", "--sample_dir", self.sample_workspace]
                uF.runSingleJob(cleanup_cmd, "cleanup", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                uF.check_completion(cleanup_checkpoint_file)
