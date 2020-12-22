import os
import sys
import subprocess
import copy
import time
import datetime

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
        self.progress = {}

    def runWorkflow(self):
        from seQuoia.other import usefulFunctions as uF

        ##################################
        ### Open for User Modification ###
        ##################################

        illumina_reads_to_subsample = 2500000
        nanomerge_memory = 24
        nanomerge_timelimit = "02:00:00"
        fastqfilter_options = uF.wrapAppropriately("-l 3000 -L 20000 -b 300000000")
        nanotrim_memory = 24
        nanotrim_timelimit = "10:00:00"
        nanosubsample_memory = 32
        nanosubsample_timelimit = "10:00:00"
        unicycler_options = uF.wrapAppropriately("--mode normal --verbosity 2")
        unicycler_threads = 4
        unicycler_memory = 16 # per cpu
        unicycler_timelimit = "48:00:00"
        pilonpolishing_max_iterations = 10
        pilonpolishing_threads = 1 #4
        pilonpolishing_memory = 64 #8
        pilonpolishing_timelimit = "24:00:00"
        run_guinan = False
        contig_size_filter = 0
        run_gaemr_ont = False
        gaemr_formatter_options = uF.wrapAppropriately("-g 1 -c 100 -r")
        gaemr_qc_options = uF.wrapAppropriately("--force --analyze_rna")
        gaemr_threads= 2 #4
        gaemr_memory = 32 #16 # per cpu
        gaemr_timelimit = "24:00:00"
        run_canu = False
        canu_options = uF.wrapAppropriately("stopOnReadQuality=false genomeSize=5m")
        canu_threads = 4
        canu_memory = 24
        canu_timelimit = "07:00:00:00"
        run_nanopolish = False
        nanopolish_options = ""
        nanopolish_threads = 16
        nanopolish_memory = 10
        nanopolish_timelimit = "72:00:00"
        run_pilon_on_canu = True
        run_reorganize = False


        arguments = {'illumina_reads_to_subsample': illumina_reads_to_subsample,
                     'nanomerge_memory': nanomerge_memory,
                     'nanomerge_timelimit': nanomerge_timelimit,
                     'fastqfilter_options': fastqfilter_options,
                     'nanotrim_memory': nanotrim_memory,
                     'nanotrim_timelimit': nanotrim_timelimit,
                     'nanosubsample_memory': nanosubsample_memory,
                     'nanosubsample_timelimit': nanosubsample_timelimit,
                     'unicycler_threads': unicycler_threads,
                     'unicycler_memory': unicycler_memory,
                     'unicycler_options': unicycler_options,
                     'unicycler_timelimit': unicycler_timelimit,
                     'pilonpolishing_max_iterations': pilonpolishing_max_iterations,
                     'pilonpolishing_threads': pilonpolishing_threads,
                     'pilonpolishing_memory': pilonpolishing_memory,
                     'pilonpolishing_timelimit': pilonpolishing_timelimit,
                     'run_guinan': run_guinan,
                     'contig_size_filter': contig_size_filter,
                     'run_gaemr_ont': run_gaemr_ont,
                     'gaemr_formatter_options': gaemr_formatter_options,
                     'gaemr_qc_options' : gaemr_qc_options,
                     'gaemr_threads' : gaemr_threads,
                     'gaemr_memory' : gaemr_memory,
                     'gaemr_timelimit': gaemr_timelimit,
                     'canu_options': canu_options,
                     'canu_threads': canu_threads,
                     'canu_memory': canu_memory,
                     'canu_timelimit': canu_timelimit,
                     'run_canu': run_canu,
                     'run_nanopolish': run_nanopolish,
                     'nanopolish_options': nanopolish_options,
                     'nanopolish_threads': nanopolish_threads,
                     'nanopolish_memory': nanopolish_memory,
                     'nanopolish_timelimit': nanopolish_timelimit,
                     'run_pilon_on_canu': run_pilon_on_canu,
                     'run_reorganize': run_reorganize}

        # overwrite variable values with config provided values
        if self.config_parameters:
            for var in self.config_parameters.items():
                if len(var[1].split()) > 1: arguments[var[0]] = uF.wrapAppropriately(var[1])
                else: arguments[var[0]] = var[1]

        boolean_arguments = ['run_canu', 'run_reorganize', 'run_nanopolish', 'run_pilon_on_canu', 'run_gaemr_ont', 'run_guinan']
        for ba in boolean_arguments:
            if str(arguments[ba]).lower() in uF.true_set: arguments[ba] = True
            else: arguments[ba] = False

        # Necessities for plumbing features!
        workflow_steps = ['setup', 'gp_to_fastq', 'symlink', 'il_adapter_trim', 'il_subsample', 'ont_merge',
                          'ont_adapter_trim', 'ont_subsample', 'uni_assembly_1', 'pipol_1', 'adapt_screen_1',
                          'gaemr_1', 'mlst_1', 'uni_assembly_2', 'pipol_2', 'adapt_screen_2', 'gaemr_2', 'mlst_2',
                          'canu', 'nanopolish', 'pipol_3', 'adapt_screen_3', 'gaemr_3', 'mlst_3', 'completion']

        checkpoint_files = ['SETUP', 'GPPROCESS', 'SYMLINK_INPUT', 'ADAPTERTRIM', 'SUBSAMPLE', 'NANOMERGE',
                            'NANOTRIM', 'NANOSAMPLE', 'UNICYCLER_ASSEMBLY_full-np', 'PILON_POLISHING_full-np',
                            'ASSEMBLY_ADAPTER_REMOVAL_full-np', 'GAEMR_full-np', 'ASSEMBLY_MLST_full-np',
                            'UNICYCLER_ASSEMBLY_sub-np', 'PILON_POLISHING_sub-np', 'ASSEMBLY_ADAPTER_REMOVAL_sub-np',
                            'GAEMR_sub-np', 'ASSEMBLY_MLST_sub-np', 'CANU_ASSEMBLY', 'NANOPOLISH', 'PILON_POLISHING_canu',
                            'ASSEMBLY_ADAPTER_REMOVAL_canu', 'GAEMR_canu', 'ASSEMBLY_MLST_canu', 'COMPLETION']

        dependent_steps = [[], ['setup'], ['setup'], ['symlink'], ['il_adapter_trim'], ['setup'],
                           ['ont_merge'], ['ont_adapter_trim'], ['ont_adapter_trim'], ['uni_assembly_1'], ['pipol_1'],
                           ['adapt_screen_1'], ['adapt_screen_1'], ['ont_subsample'], ['uni_assembly_2'], ['pipol_2'],
                           ['adapt_screen_2'], ['adapt_screen_2'], ['ont_adapter_trim'], ['canu'], ['canu'], ['canu'],
                           ['adapt_screen_3'], ['adapt_screen_3'], ['gaemr_1', 'gaemr_2', 'mlst_1', 'mlst_2']]

        for i, step in enumerate(workflow_steps):
            dependency_status = 'not_ready'
            if i == 0: dependency_status = 'ready'
            checkpoint_file = self.sample_workspace + checkpoint_files[i] + '.txt'
            self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i]), dependency_status, None]
            if step == 'symlink' and self.bamLocation:
                self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['gp_to_fastq']), dependency_status, None]
            if (step == 'uni_assembly_1' or step == 'uni_assembly_2') and (self.bamLocation or (self.forwardFastq and self.reverseFastq)):
                self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['il_subsample']), dependency_status, None]
            if step == 'pipol_3':
                if (self.bamLocation or (self.forwardFastq and self.reverseFastq)):
                    self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['il_subsample']), dependency_status, None]
                if arguments['run_nanopolish']:
                    self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['nanopolish']), dependency_status, None]
            if step == 'adapt_screen_3':
                if arguments['run_nanopolish']:
                    self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['nanopolish']), dependency_status, None]
                if arguments['run_pilon_on_canu']:
                    self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['pipol_3']), dependency_status, None]
            if step == 'completion' and arguments['run_canu']:
                self.progress[step] = ['not_started', checkpoint_file, set(dependent_steps[i] + ['gaemr_3', 'mlst_3']), dependency_status, None]

        # loop to watch for job completion and update progress

        if not os.path.isfile(self.sample_workspace + 'REORGANIZATION.txt'):
            in_progress = True
            while in_progress:
                self.updateProgress()
                in_progress = not (self.progress['completion'][0] == 'completed')
                try: self.workflowInstance(arguments)
                except: sys.stderr.write("Error in workflow for sample %s!\n" % self.sample); break
                time.sleep(30)

    def updateProgress(self):
        from seQuoia.other import usefulFunctions as uF

        try:
            iter = 0
            while iter < 2:
                qstat_cmd = ['qstat']
                if self.cluster == 'SLURM':
                    qstat_cmd = ['squeue']
                proc = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                user_running_jobs = set([])
                for i, line in enumerate(out.decode('utf-8').split('\n')):
                    line = line.strip()
                    ls = line.split()
                    if (i > 1 or (i > 0 and self.cluster == 'SLURM'))and len(ls) > 4:
                        if uF.is_number(ls[0]) and not 'E' in ls[4]:
                            user_running_jobs.add(ls[0])

                progress_copy = copy.deepcopy(self.progress)
                for s in self.progress.items():
                    # first update job status
                    s_values = s[1][0]
                    job_id = s[1][-1]
                    checkpoint_file = s[1][1]
                    if s_values == 'completed': continue
                    elif uF.checkIfFileExists(checkpoint_file): progress_copy[s[0]][0] = 'completed'; continue
                    elif job_id:
                        if job_id in user_running_jobs: progress_copy[s[0]][0] = 'in_progress'; continue
                        else:
                            if uF.checkIfFileExists(checkpoint_file): progress_copy[s[0]][0] = 'completed'; continue
                            elif iter == 0: progress_copy[s[0]][0] = 'warning'; time.sleep(15); continue
                            elif iter == 1 and s_values != 'warning':
                                iter -= 1
                                progress_copy[s[0]][0] = 'in_progress'
                                continue
                            else:
                                time.sleep(300)
                                if uF.checkIfFileExists(checkpoint_file): progress_copy[s[0]][0] = 'completed'; continue
                                # give up hope and kill any remaining workflow jobs
                                all_workflow_jobs = set([x[-1] for x in self.progress.values()])
                                ongoing_jobs = user_running_jobs.intersection(all_workflow_jobs)
                                if len(ongoing_jobs) > 0:
                                    qdel_cmd = ['qdel'] + list(ongoing_jobs)
                                    if self.cluster == 'SLURM':
                                        qdel_cmd = ['scancel'] + list(ongoing_jobs)
                                    proc = subprocess.Popen(qdel_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                    proc.wait()

                                sys.stderr.write('Warning with sample %s on workflow step %s\n' % (self.sample, s[0]))
                                sys.exit(1)
                user_running_jobs = None
                self.progress = progress_copy
                iter += 1

            progress_copy = copy.deepcopy(self.progress)
            for s in self.progress.items():
                # second update job dependency status
                if s[1][3] == 'ready': continue
                dependencies = s[1][2]
                dependency_stati = 'ready'
                for ds in dependencies:
                    ds_status = self.progress[ds][0]
                    if not ds_status == 'completed': dependency_stati = 'not_ready'
                progress_copy[s[0]][3] = dependency_stati

            self.progress = progress_copy
        except:
            raise RuntimeWarning()
            sys.exit(1)

    def workflowInstance(self, arguments):
        from seQuoia.other import usefulFunctions as uF

        ####
        ## Start the pipeline!
        ####

        # Step 1: Run Setup_task module
        ###################################################################################################################################
        if self.progress['setup'][0] == 'not_started' and self.progress['setup'][3] == 'ready':
            setup_cmd = ["sh", self.code_dir + "Setup_task.sh", "--sample", self.sample,
                                                                "--parent_dir", self.workspace,
                                                                "--meta", self.metadata]
            job_id = uF.runSingleJobBackground(setup_cmd, "setup", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
            self.progress['setup'][-1] = job_id

        if self.bamLocation or self.singleFastq or (self.forwardFastq and self.reverseFastq):
            """ STARTING ILLUMINA SPECIFIC QC/PROCCESSING MODULES"""
            if self.bamLocation and self.progress['gp_to_fastq'][0] == 'not_started' and self.progress['gp_to_fastq'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'ProcessGPDirectory/')
                gp_process_cmd = ["sh", self.code_dir + "ProcessGPDirectory_task.sh", "--bam_location=" + self.bamLocation,
                                                                                      "--sample=" + self.sample,
                                                                                      "--sample_dir=" + self.sample_workspace]
                job_id = uF.runSingleJobBackground(gp_process_cmd, "gp_process", self.sample, self.grid_dir, mem=32, timelimit="12:00:00", cluster=self.cluster, project=self.project)
                self.progress['gp_to_fastq'][-1] = job_id

            if self.progress['gp_to_fastq'][0] == 'completed': self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("gp_process", self.sample_workspace, self.singleFastq, self.forwardFastq, self.reverseFastq)

            # Step 1.5: Run Input Symlinking (Dumb module to just symlink illumina sequencing reads to have sample names in files)
            ###################################################################################################################################
            if self.progress['symlink'][0] == 'not_started' and self.progress['symlink'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Symlink_Input/')
                symlink_cmd = ["sh", self.code_dir + "SymlinkInput_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                             "--fastq_frw=" + self.forwardFastq,
                                                                             "--fastq_rev=" + self.reverseFastq,
                                                                             "--sample_dir=" + self.sample_workspace,
                                                                             "--sample=" + self.sample]
                job_id = uF.runSingleJobBackground(symlink_cmd, "symlink_input", self.sample, self.grid_dir, cluster=self.cluster, project=self.project)
                self.progress['symlink'][-1] = job_id

            if self.progress['symlink'][0] == 'completed': self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("symlink", self.sample_workspace, self.singleFastq, self.forwardFastq, self.reverseFastq)

            # Step 2: Run Adapter Trimming
            ###################################################################################################################################
            if self.progress['il_adapter_trim'][0] == 'not_started' and self.progress['il_adapter_trim'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'AdapterTrim/')
                cutadapt_cmd = ["sh", self.code_dir + "AdapterTrim_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                             "--fastq_frw=" + self.forwardFastq,
                                                                             "--fastq_rev=" + self.reverseFastq,
                                                                             "--sample=" + self.sample,
                                                                             "--sample_dir=" + self.sample_workspace]
                job_id = uF.runSingleJobBackground(cutadapt_cmd, "adapter_trim", self.sample, self.grid_dir, cluster=self.cluster, mem=5, threads=1, project=self.project)
                self.progress['il_adapter_trim'][-1] = job_id

            if self.progress['il_adapter_trim'][0] == 'completed': self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("adapter_trim", self.sample_workspace, self.singleFastq, self.forwardFastq, self.reverseFastq)

            # Step 3: Run Subsample
            ###################################################################################################################################
            if self.progress['il_subsample'][0] == 'not_started' and self.progress['il_subsample'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Subsample/')
                subsample_cmd = ["sh", self.code_dir + "Subsample_task.sh", "--fastq_sin=" + self.singleFastq,
                                                                            "--fastq_frw=" + self.forwardFastq,
                                                                            "--fastq_rev=" + self.reverseFastq,
                                                                            "--sample=" + self.sample,
                                                                            "--sample_dir=" + self.sample_workspace,
                                                                            "--reads=" + str(arguments['illumina_reads_to_subsample'])]
                job_id = uF.runSingleJobBackground(subsample_cmd, "subsample", self.sample, self.grid_dir, mem=24, cluster=self.cluster, project=self.project)
                self.progress['il_subsample'][-1] = job_id

            if self.progress['il_subsample'][0] == 'completed': self.singleFastq, self.forwardFastq, self.reverseFastq = uF.extractResultingFastqs("subsample", self.sample_workspace, self.singleFastq, self.forwardFastq, self.reverseFastq)


        # Step 4: Run Nanopore Merging
        ###################################################################################################################################
        if self.progress['ont_merge'][0] == 'not_started' and self.progress['ont_merge'][3] == 'ready':
            uF.cleanUp(self.sample_workspace + 'NanoMerge/')
            nanomerge_cmd = ["sh", self.code_dir + "NanoMerge_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                        "--sample_dir=" + self.sample_workspace,
                                                                        "--sample=" + self.sample,
                                                                        "--barcode=" + self.nanoporeBarcode]
            job_id = uF.runSingleJobBackground(nanomerge_cmd, "nanomerge", self.sample, self.grid_dir, mem=int(arguments['nanomerge_memory']), timelimit=arguments['nanomerge_timelimit'], cluster=self.cluster, project=self.project)
            self.progress['ont_merge'][-1] = job_id
            
        if self.progress['ont_merge'][0] == 'completed': 
            self.nanoporeFastq = self.sample_workspace + '/NanoMerge/' + self.sample + '.fastq'
            if not os.path.isfile(self.nanoporeFastq): self.nanoporeFastq += '.gz'

        # Step 5: Run Nanopore Adapter Trimmming
        ###################################################################################################################################
        if self.progress['ont_adapter_trim'][0] == 'not_started' and self.progress['ont_adapter_trim'][3] == 'ready':
            uF.cleanUp(self.sample_workspace + 'NanoTrim/')
            nanomerge_cmd = ["sh", self.code_dir + "NanoTrim_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                       "--sample_dir=" + self.sample_workspace,
                                                                       "--sample=" + self.sample]
            job_id = uF.runSingleJobBackground(nanomerge_cmd, "nanotrim", self.sample, self.grid_dir, mem=int(arguments['nanotrim_memory']), timelimit=str(arguments['nanotrim_timelimit']), cluster=self.cluster, project=self.project)
            self.progress['ont_adapter_trim'][-1] = job_id
            
        if self.progress['ont_adapter_trim'][0] == 'completed':
            self.nanoporeFastq = self.sample_workspace + '/NanoTrim/' + self.sample + '.fastq'
            if not os.path.isfile(self.nanoporeFastq): self.nanoporeFastq += '.gz'

        # Step 6: Run Nanopore Subsample
        ###################################################################################################################################
        if self.progress['ont_subsample'][0] == 'not_started' and self.progress['ont_subsample'][3] == 'ready':
            uF.cleanUp(self.sample_workspace + 'NanoSample/')
            nanosample_cmd = ["sh", self.code_dir + "NanoSample_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                          "--fastqfilter_parameters=" + arguments['fastqfilter_options'],
                                                                          "--sample_dir="+ self.sample_workspace,
                                                                          "--sample=" + self.sample]
            job_id = uF.runSingleJobBackground(nanosample_cmd, "nanosample", self.sample, self.grid_dir, mem=int(arguments['nanosubsample_memory']), timelimit=str(arguments['nanosubsample_timelimit']), cluster=self.cluster, project=self.project)
            self.progress['ont_subsample'][-1] = job_id
            
        nanoporeFastq_subsampled = ''
        if self.progress['ont_subsample'][0] == 'completed':
            nanoporeFastq_subsampled = self.sample_workspace + 'NanoSample/' + self.sample + '.fastq'
            if not os.path.isfile(nanoporeFastq_subsampled): nanoporeFastq_subsampled += '.gz'

        if self.forwardFastq and self.reverseFastq:
            """ Assembly 1 - All ONT Reads Based Assembly """


            # Run 1 Identifier
            run1_id = "full-np"

            # Step 7-1: Run Unicycler
            ###################################################################################################################################
            if self.progress['uni_assembly_1'][0] == 'not_started' and self.progress['uni_assembly_1'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Unicycler_Assembly_' + run1_id + '/')
                nanounicyc_cmd = ["sh", self.code_dir + "NanoUnicycler_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                                  "--illumina_forward=" + self.forwardFastq,
                                                                                  "--illumina_reverse=" + self.reverseFastq,
                                                                                  "--unicycler_parameters=" + str(arguments['unicycler_options']),
                                                                                  "--sample_dir=" + self.sample_workspace,
                                                                                  "--sample=" + self.sample,
                                                                                  "--identifier=" + run1_id,
                                                                                  "--cores=" + str(arguments['unicycler_threads'])]
                job_id = uF.runSingleJobBackground(nanounicyc_cmd, "unicycler-" + run1_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=str(arguments['unicycler_timelimit']), threads=int(arguments['unicycler_threads']), mem=int(arguments['unicycler_memory']), project=self.project)
                self.progress['uni_assembly_1'][-1] = job_id

            assembly_1 = self.sample_workspace + 'Unicycler_Assembly_' + run1_id + '/' + self.sample + '.fasta'

            # Step 8-1: Run Iterative Pilon Polishing
            if self.progress['pipol_1'][0] == 'not_started' and self.progress['pipol_1'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Pilon_Polishing_' + run1_id + '/')
                pipol_cmd = ["sh", self.code_dir + "PilonPolishing_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                              "--fastq_frw=" + self.forwardFastq,
                                                                              "--fastq_rev=" + self.reverseFastq,
                                                                              "--assembly=" + assembly_1,
                                                                              "--sample=" + self.sample,
                                                                              "--sample_dir=" + self.sample_workspace,
                                                                              "--max_iterations=" + str(arguments['pilonpolishing_max_iterations']),
                                                                              "--identifier=" + run1_id,
                                                                              "--cores=" + str(arguments['pilonpolishing_threads'])]
                job_id = uF.runSingleJobBackground(pipol_cmd, "pilonpolishing-" + run1_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['pilonpolishing_timelimit'], threads=int(arguments['pilonpolishing_threads']), mem=int(arguments['pilonpolishing_memory']), project=self.project)
                self.progress['pipol_1'][-1] = job_id

            assembly_1 = self.sample_workspace + '/Pilon_Polishing_' + run1_id + '/' + self.sample + '.pilon-refined.fasta'

            # Step 9-1: Run Adapter Contamination Screen
            ###################################################################################################################################
            if self.progress['adapt_screen_1'][0] == 'not_started' and self.progress['adapt_screen_1'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_Adapter_Removal_' + run1_id + '/')
                aar_cmd = ["sh", self.code_dir + "AssemblyAdapterRemoval_task.sh", "--assembly=" + assembly_1,
                                                                                     "--sample_dir=" + self.sample_workspace,
                                                                                     "--sample=" + self.sample,
                                                                                     "--size_filter=" + str(arguments['contig_size_filter']),
                                                                                     "--identifier=" + run1_id]
                if arguments['run_guinan']: aar_cmd += ['--run_guinan']

                job_id = uF.runSingleJobBackground(aar_cmd, "adapter-contam-screen-" + run1_id, self.sample, self.grid_dir,
                                                   cluster=self.cluster, timelimit="01:59:00",
                                                   threads=1,
                                                   mem=10, project=self.project)
                self.progress['adapt_screen_1'][-1] = job_id

            assembly_1 = self.sample_workspace + '/Assembly_Adapter_Removal_' + run1_id + '/' + self.sample + '.filtered.fasta'

            # Step 10-1: Run GAEMR
            ###################################################################################################################################
            if self.progress['gaemr_1'][0] == 'not_started' and self.progress['gaemr_1'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'GAEMR_' + run1_id + '/')
                gaemr_cmd = ["sh", self.code_dir + "GAEMR_task.sh", "--assembly=" + assembly_1,
                                                                      "--sample_dir=" + self.sample_workspace,
                                                                      "--sample=" + self.sample,
                                                                      "--format_options=" + str(arguments['gaemr_formatter_options']),
                                                                      "--qc_options=" + str(arguments['gaemr_qc_options']),
                                                                      "--identifier=" + run1_id,
                                                                      "--nanopore_fastq=" + self.nanoporeFastq,
                                                                      "--cores=" + str(arguments['gaemr_threads'])]
                if arguments['run_gaemr_ont']: gaemr_cmd += ['--gaemr_ont']
                job_id = uF.runSingleJobBackground(gaemr_cmd, "gaemr-" + run1_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['gaemr_timelimit'], threads= int(arguments['gaemr_threads']), mem= int(arguments['gaemr_memory']), project=self.project)
                self.progress['gaemr_1'][-1] = job_id

            # Step 11-1: Run MLST on Assembly
            ###################################################################################################################################
            if self.progress['mlst_1'][0] == 'not_started' and self.progress['mlst_1'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_MLST_' + run1_id + '/')
                mlst_cmd = ["sh", self.code_dir + "AssemblyMLST_task.sh", "--assembly=" + assembly_1,
                                                                          "--sample_dir=" + self.sample_workspace,
                                                                          "--sample=" + self.sample,
                                                                          "--identifier=" + run1_id]

                job_id = uF.runSingleJobBackground(mlst_cmd, "assembly-mlst-" + run1_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit="01:59:00", threads=1, mem=10, project=self.project)
                self.progress['mlst_1'][-1] = job_id

            """ Assembly 2 - Subsampled ONT Reads Based Assembly """

            # Run 2 Identifier
            run2_id = "sub-np"

            # Step 7-2: Run Unicycler
            ###################################################################################################################################
            if self.progress['uni_assembly_2'][0] == 'not_started' and self.progress['uni_assembly_2'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Unicycler_Assembly_' + run2_id + '/')
                nanounicyc_cmd = ["sh", self.code_dir + "NanoUnicycler_task.sh", "--nanopore_fastq=" + nanoporeFastq_subsampled,
                                                                                  "--illumina_forward=" + self.forwardFastq,
                                                                                  "--illumina_reverse=" + self.reverseFastq,
                                                                                  "--unicycler_parameters=" + str(arguments['unicycler_options']),
                                                                                  "--sample_dir=" + self.sample_workspace,
                                                                                  "--sample=" + self.sample,
                                                                                  "--identifier=" + run2_id,
                                                                                  "--cores=" + str(arguments['unicycler_threads']) ]
                job_id = uF.runSingleJobBackground(nanounicyc_cmd, "unicycler-" + run2_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=str(arguments['unicycler_timelimit']), threads=int(arguments['unicycler_threads']), mem=int(arguments['unicycler_memory']), project=self.project)
                self.progress['uni_assembly_2'][-1] = job_id

            assembly_2 = self.sample_workspace + '/Unicycler_Assembly_' + run2_id + '/' + self.sample + '.fasta'

            # Step 8-2: Run Iterative Pilon Polishing
            if self.progress['pipol_2'][0] == 'not_started' and self.progress['pipol_2'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Pilon_Polishing_' + run2_id + '/')
                pipol_cmd = ["sh", self.code_dir + "PilonPolishing_task.sh", "--nanopore_fastq=" + nanoporeFastq_subsampled,
                                                                              "--fastq_frw=" + self.forwardFastq,
                                                                              "--fastq_rev=" + self.reverseFastq,
                                                                              "--assembly=" + assembly_2,
                                                                              "--sample=" + self.sample,
                                                                              "--sample_dir=" + self.sample_workspace,
                                                                              "--max_iterations=" + str(arguments['pilonpolishing_max_iterations']),
                                                                              "--identifier=" + run2_id,
                                                                              "--cores=" + str(arguments['pilonpolishing_threads'])]
                job_id = uF.runSingleJobBackground(pipol_cmd, "pilonpolishing-" + run2_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['pilonpolishing_timelimit'], threads=int(arguments['pilonpolishing_threads']), mem=int(arguments['pilonpolishing_memory']), project=self.project)
                self.progress['pipol_2'][-1] = job_id

            assembly_2 = self.sample_workspace + '/Pilon_Polishing_' + run2_id + '/' + self.sample + '.pilon-refined.fasta'

            # Step 9-2: Run Adapter Contamination Screen
            ###################################################################################################################################
            if self.progress['adapt_screen_2'][0] == 'not_started' and self.progress['adapt_screen_2'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_Adapter_Removal_' + run2_id + '/')
                aar_cmd = ["sh", self.code_dir + "AssemblyAdapterRemoval_task.sh", "--assembly=" + assembly_2,
                                                                                 "--sample_dir=" + self.sample_workspace,
                                                                                 "--sample=" + self.sample,
                                                                                 "--size_filter="+  str(arguments['contig_size_filter']),
                                                                                 "--identifier=" + run2_id]
                if arguments['run_guinan']: aar_cmd += ['--run_guinan']
                job_id = uF.runSingleJobBackground(aar_cmd, "adapter-contam-screen-" + run2_id, self.sample, self.grid_dir,
                                                   cluster=self.cluster, timelimit="01:59:00",
                                                   threads=1,
                                                   mem=10, project=self.project)
                self.progress['adapt_screen_2'][-1] = job_id

            assembly_2 = self.sample_workspace + '/Assembly_Adapter_Removal_' + run2_id + '/' + self.sample + '.filtered.fasta'

            # Step 10-2: Run GAEMR
            ###################################################################################################################################
            if self.progress['gaemr_2'][0] == 'not_started' and self.progress['gaemr_2'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'GAEMR_' + run2_id + '/')
                gaemr_cmd = ["sh", self.code_dir + "GAEMR_task.sh", "--assembly=" + assembly_2,
                                                                      "--sample_dir=" + self.sample_workspace,
                                                                      "--sample=" + self.sample,
                                                                      "--format_options=" + str(arguments['gaemr_formatter_options']),
                                                                      "--qc_options=" + str(arguments['gaemr_qc_options']),
                                                                      "--identifier=" + run2_id,
                                                                      "--nanopore_fastq=" + nanoporeFastq_subsampled,
                                                                      "--cores=" + str(arguments['gaemr_threads'])]
                if arguments['run_gaemr_ont']: gaemr_cmd += ['--gaemr_ont']
                job_id = uF.runSingleJobBackground(gaemr_cmd, "gaemr-" + run2_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['gaemr_timelimit'], threads= int(arguments['gaemr_threads']), mem= int(arguments['gaemr_memory']), project=self.project)
                self.progress['gaemr_2'][-1] = job_id

            # Step 11-2: Run MLST on Assembly
            ###################################################################################################################################
            if self.progress['mlst_2'][0] == 'not_started' and self.progress['mlst_2'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_MLST_' + run2_id + '/')
                mlst_cmd = ["sh", self.code_dir + "AssemblyMLST_task.sh", "--assembly=" + assembly_2,
                                                                          "--sample_dir=" + self.sample_workspace,
                                                                          "--sample=" + self.sample,
                                                                          "--identifier=" + run2_id]

                job_id = uF.runSingleJobBackground(mlst_cmd, "assembly-mlst-" + run2_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit="01:59:00", threads=1, mem=10, project=self.project)
                self.progress['mlst_2'][-1] = job_id

        if arguments['run_canu']:
            run3_id = "canu"
            # Step 7-3: Run Canu
            ###################################################################################################################################
            if self.progress['canu'][0] == 'not_started' and self.progress['canu'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Canu_Assembly/')
                nanocanu_cmd = ["sh", self.code_dir + "NanoCanu_task.sh", "--nanopore_fastq=" + self.nanoporeFastq,
                                                                          "--canu_parameters=" + str(arguments['canu_options']),
                                                                          "--sample_dir=" + self.sample_workspace,
                                                                          "--sample=" + self.sample,
                                                                          "--memory=" + str(arguments['canu_memory']),
                                                                          "--cores=" + str(arguments['canu_threads'])]
                job_id = uF.runSingleJobBackground(nanocanu_cmd, "canu", self.sample, self.grid_dir, cluster=self.cluster, timelimit=str(arguments['canu_timelimit']), threads=int(arguments['canu_threads']), mem=int(arguments['canu_memory']), project=self.project)
                self.progress['canu'][-1] = job_id

            assembly_3 = self.sample_workspace + '/Canu_Assembly/' + self.sample + '.fasta'

            if arguments['run_nanopolish']:
                # Step 8-3: Run Nanopolish
                ###################################################################################################################################
                if self.progress['nanopolish'][0] == 'not_started' and self.progress['nanopolish'][3] == 'ready':
                    uF.cleanUp(self.sample_workspace + 'Nanopolish/')
                    nanocanu_cmd = ["sh", self.code_dir + "NanoPolish_task.sh", "--assembly=" + assembly_3,
                                                                                  "--nanopore_fastq=" + self.nanoporeFastq,
                                                                                  "--nanopore_fast5=" + self.nanoporeFast5,
                                                                                  "--nanopore_barcode=" + self.nanoporeBarcode,
                                                                                  "--sample_dir=" + self.sample_workspace,
                                                                                  "--sample=" + self.sample,
                                                                                  "--nanopolish_variants_options=" + str(arguments['nanopolish_options']),
                                                                                  "--cores=" + str(arguments['nanopolish_threads'])]
                    job_id = uF.runSingleJobBackground(nanocanu_cmd, "nanopolish", self.sample, self.grid_dir, cluster=self.cluster, timelimit=str(arguments['nanopolish_timelimit']), threads=int(arguments['nanopolish_threads']), mem=int(arguments['nanopolish_memory']), project=self.project)
                    self.progress['nanopolish'][-1] = job_id

                assembly_3 = self.sample_workspace + '/Nanopolish/polishing/' + self.sample + '.fasta'

            # Step 9-3: Run Iterative Pilon Polishing
            if arguments['run_pilon_on_canu']:
                if self.progress['pipol_3'][0] == 'not_started' and self.progress['pipol_3'][3] == 'ready':
                    uF.cleanUp(self.sample_workspace + 'Pilon_Polishing_' + run3_id + '/')
                    pipol_cmd = ["sh", self.code_dir + "PilonPolishing_task.sh", "--nanopore_fastq=" + nanoporeFastq_subsampled,
                                                                                  "--fastq_frw=" + self.forwardFastq,
                                                                                  "--fastq_rev=" + self.reverseFastq,
                                                                                  "--sample=" + self.sample,
                                                                                  "--assembly=" + assembly_3,
                                                                                  "--sample_dir=" + self.sample_workspace,
                                                                                  "--max_iterations=" + str(arguments['pilonpolishing_max_iterations']),
                                                                                  "--identifier=" + run3_id,
                                                                                  "--cores=" + str(arguments['pilonpolishing_threads'])]
                    job_id = uF.runSingleJobBackground(pipol_cmd, "pilonpolishing-" + run3_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['pilonpolishing_timelimit'], threads=int(arguments['pilonpolishing_threads']), mem=int(arguments['pilonpolishing_memory']), project=self.project)
                    self.progress['pipol_3'][-1] = job_id

                assembly_3 = self.sample_workspace + '/Pilon_Polishing_' + run3_id + '/' + self.sample + '.pilon-refined.fasta'

            # Step 10-3: Run Adapter Contamination Screen
            ###################################################################################################################################
            if self.progress['adapt_screen_3'][0] == 'not_started' and self.progress['adapt_screen_3'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_Adapter_Removal_' + run3_id + '/')
                aar_cmd = ["sh", self.code_dir + "AssemblyAdapterRemoval_task.sh", "--assembly=" + assembly_3,
                                                                                     "--sample_dir=" + self.sample_workspace,
                                                                                     "--sample=" + self.sample,
                                                                                     "--size_filter=" + str(arguments['contig_size_filter']),
                                                                                     "--identifier=" + run3_id]
                if arguments['run_guinan']: aar_cmd += ['--run_guinan']
                job_id = uF.runSingleJobBackground(aar_cmd, "adapter-contam-screen-" + run3_id, self.sample, self.grid_dir,
                                                   cluster=self.cluster, timelimit="01:59:00",
                                                   threads=1,
                                                   mem=10, project=self.project)
                self.progress['adapt_screen_3'][-1] = job_id

            assembly_3 = self.sample_workspace + '/Assembly_Adapter_Removal_' + run3_id + '/' + self.sample + '.filtered.fasta'

            # Step 11-3: Run GAEMR
            ###################################################################################################################################
            if self.progress['gaemr_3'][0] == 'not_started' and self.progress['gaemr_3'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'GAEMR_' + run3_id + '/')
                gaemr_cmd = ["sh", self.code_dir + "GAEMR_task.sh", "--assembly=" + assembly_3,
                                                                     "--sample_dir=" + self.sample_workspace,
                                                                     "--sample=" + self.sample,
                                                                     "--format_options=" + str(arguments['gaemr_formatter_options']),
                                                                     "--qc_options=" + str(arguments['gaemr_qc_options']),
                                                                     "--identifier=" + run3_id,
                                                                     "--nanopore_fastq=" + self.nanoporeFastq,
                                                                     "--cores=" + str(arguments['gaemr_threads'])]
                if arguments['run_gaemr_ont']: gaemr_cmd += ['--gaemr_ont']
                job_id = uF.runSingleJobBackground(gaemr_cmd, "gaemr-" + run3_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit=arguments['gaemr_timelimit'], threads=int(arguments['gaemr_threads']), mem=int(arguments['gaemr_memory']), project=self.project)
                self.progress['gaemr_3'][-1] = job_id

            # Step 12-3: Run MLST on Assembly
            ###################################################################################################################################
            if self.progress['mlst_3'][0] == 'not_started' and self.progress['mlst_3'][3] == 'ready':
                uF.cleanUp(self.sample_workspace + 'Assembly_MLST_' + run3_id + '/')
                mlst_cmd = ["sh", self.code_dir + "AssemblyMLST_task.sh", "--assembly="+ assembly_3,
                                                                          "--sample_dir="+ self.sample_workspace,
                                                                          "--sample="+ self.sample,
                                                                          "--identifier=" + run3_id]

                job_id = uF.runSingleJobBackground(mlst_cmd, "assembly-mlst-" + run3_id, self.sample, self.grid_dir, cluster=self.cluster, timelimit="01:59:00", threads=1, mem=10, project=self.project)
                self.progress['mlst_3'][-1] = job_id

        # Step 10: Reorganize the Sample Repository
        ###################################################################################################################################
        if (self.progress['completion'][0] == 'not_started' and self.progress['completion'][3] == 'ready') or arguments['run_reorganize']:
            os.system('rm %s' % self.sample_workspace + 'COMPLETION.txt')
            completion_cmd = ["sh", self.code_dir + "BacONTCompletionCheck_task.sh", "--sample_dir=" + self.sample_workspace,
                                                                                      "--sample=" + self.sample]
            if arguments['run_reorganize']:
                completion_cmd += ['--reorganize']
            job_id = uF.runSingleJobBackground(completion_cmd, "completion", self.sample, self.grid_dir, cluster=self.cluster,
                                               project=self.project)
            self.progress['completion'][-1] = job_id
