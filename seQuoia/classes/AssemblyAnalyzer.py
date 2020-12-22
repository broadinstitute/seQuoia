import os
import sys
import subprocess
from Bio import SeqIO
import multiprocessing

external_wrappers_dir = '/'.join(
    os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

def run_nanopolish_variants(data):
    workspace, window, cmd = data
    stdout_file_1 = workspace + window.replace(':', '_') + '.log'
    stderr_file_1 = workspace + window.replace(':', '_') + '.err'
    stdout_handle_1 = open(stdout_file_1, 'w')
    stderr_handle_1 = open(stderr_file_1, 'w')
    proc = subprocess.Popen(cmd, stdout=stdout_handle_1, stderr=stderr_handle_1)
    proc.wait()
    stdout_handle_1.close();
    stderr_handle_1.close();


class Assembly():
    def __init__(self, assembly_fasta, sample_name, logObject, contig_fasta=None, assembly_graph=None):
        self.assembly_fasta = os.path.abspath(assembly_fasta)
        self.sample_name = sample_name
        self.logFile = logObject
        self.contig_fasta = contig_fasta
        self.assembly_graph = assembly_graph

    def filter_contigs_by_size(self, working_dir, size_filter=0, change_reference=True):
        """ Filter assembly contigs by size threshold. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            assert(os.path.isdir(working_dir))
            self.logFile.info('-'*70)
            self.logFile.info("Filtering assembly %s for contigs shorter than %d." % (self.assembly_fasta, size_filter))
            resulting_assembly = working_dir + self.sample_name + '.filtered.fasta'
            assembly_handle = open(resulting_assembly, 'w')
            with open(self.assembly_fasta) as oaf:
                for rec in SeqIO.parse(oaf, 'fasta'):
                    if len(str(rec.seq)) >= size_filter:
                        assembly_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
            assembly_handle.close()
            assert(os.path.isfile(resulting_assembly))
            self.logFile.info("Successfully filtered assembly. Resulting assembly file can be found at: %s" % resulting_assembly)
            if change_reference:
                self.assembly_fasta = resulting_assembly
                self.assembly_graph = None
                self.contig_fasta = None
            return resulting_assembly
        except RuntimeError as e:
            self.logFile.error("Something went wrong when filtering assembly for short contigs. Raising exception ..." % e)
            raise RuntimeError

    def detect_adapters(self, working_dir, options='-r -e -f Directive -v --no_protection --adaptorsdb', cores=1):
        """ Use GAEMR to detect adapters in assembly. Returns commands for guinan analysis and adapter removal. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            assert(os.path.isdir(working_dir) and os.path.isfile(self.contig_fasta))
            gaemr_version_info = open(external_wrappers_dir + "run_contamination_screen.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Screening for adapters with GAEMR version %s" % gaemr_version_info)
            result_txt = working_dir + self.sample_name + '.screen.out'
            gaemr_cmd = ['/bin/bash', external_wrappers_dir + 'run_contamination_screen.sh', '-t', str(cores), options, '-o',
                        result_txt, '-i', self.contig_fasta]
            self.logFile.info("Executing the command: %s" % ' '.join(gaemr_cmd))
            stdout_file = working_dir + self.sample_name + ".contamination_screen.log"
            stderr_file = working_dir + self.sample_name + ".contamination_screen.err"
            stdout_handle = open(result_txt, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(gaemr_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            commands_result_file = working_dir + self.sample_name + '.screen.out.directives.txt'
            assert(os.path.isfile(result_txt) and os.path.isfile(commands_result_file))
            self.logFile.info("GAEMR adapter detection finished successfully!")
            return commands_result_file
        except RuntimeError as e:
            self.logFile.error("GAEMR adapter detection was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def remove_adapters(self, command_file, working_dir, options='--sort_by_length --rename', change_reference=True):
        #  Use guinan suite to remove adapters in assembly. 
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            assert(os.path.isdir(working_dir) and os.path.isfile(command_file) and os.path.isfile(self.contig_fasta) and os.path.isfile(self.assembly_graph))
            guinan_version_info = open(external_wrappers_dir + "modify_assembly.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Removing adapters with guinan version %s" % guinan_version_info)
            result_header = working_dir +  self.sample_name + '.adapter-free'
            guinan_cmd = ['/bin/bash', external_wrappers_dir + 'modify_assembly.sh', options, '--agp', self.assembly_graph,
                          '--command_file', command_file, '--output_header', result_header, self.contig_fasta]
            self.logFile.info("Executing the command: %s" % ' '.join(guinan_cmd))
            stdout_file = working_dir + self.sample_name + ".contamination_screen.log"
            stderr_file = working_dir + self.sample_name + ".contamination_screen.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(guinan_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            # TODO: improve checksafe
            result_file = result_header + '.scaffolds.fasta'
            assert(os.path.isfile(result_file))
            self.logFile.info("Successfully filtered out any adapters with guianan, filtered assembly can be found at: %s" % result_file)
            if change_reference:
                self.assembly_fasta = result_file
                self.assembly_graph = None
                self.contig_fasta = None
            self.logFile.info("guinan adapter removal finished successfully!")
            return result_file
        except RuntimeError as e:
            self.logFile.error("guinan adapter removal was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_mlst(self, working_dir):
        """ Performing mlst on assembly. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            assert(os.path.isdir(working_dir))
            mlst_version_info = open(external_wrappers_dir + "mlst.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Performing MLST of %s" % mlst_version_info)
            result_txt = working_dir + self.sample_name + '.txt'
            mlst_cmd = ['/bin/bash', external_wrappers_dir + 'mlst.sh', self.assembly_fasta]
            self.logFile.info("Executing the command: %s" % ' '.join(mlst_cmd))
            stderr_file = working_dir + self.sample_name + ".mlst.err"
            stdout_handle = open(result_txt, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(mlst_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            # TODO: improve checksafe
            assert(os.path.isfile(result_txt) and os.path.getsize(result_txt) > 100)
            self.logFile.info("mlst finished successfully!")
        except RuntimeError as e:
            self.logFile.error("mlst was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_gaemr_formatter(self, working_dir, options='-g 1 -c 100 -r', reference_change=True):
        try:
            assert (os.path.isfile(self.assembly_fasta))
            workspace = os.path.abspath(working_dir) + '/'

            gaemr_version = open(external_wrappers_dir + 'make_standard_assembly_files.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Formatting assemblies with GAEMR for QC and Vesper2, using %s" % gaemr_version)

            output_results = workspace + '/' + self.sample_name

            gaemr_cmd = ["sh", external_wrappers_dir + "make_standard_assembly_files.sh", options, '-S', self.assembly_fasta,
                         '-o', output_results]

            self.logFile.info("Executing the command: %s" % ' '.join(gaemr_cmd))
            stdout_file = workspace + "gaemr_formatting.log"
            stderr_file = workspace + "gaemr_formatting.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(gaemr_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            resulting_graph = output_results + '.agp'
            resulting_contig = output_results + '.contigs.fasta'
            resulting_scaffold = output_results + '.scaffolds.fasta'

            assert (os.path.isfile(resulting_graph) and os.path.isfile(resulting_contig) and os.path.isfile(resulting_scaffold))
            self.logFile.info("Results are available at %s" % workspace)

            if reference_change:
                self.assembly_fasta = resulting_scaffold
                self.contig_fasta = resulting_contig
                self.assembly_graph = resulting_graph
            return ([resulting_scaffold, resulting_contig, resulting_graph])
        except RuntimeError as e:
            self.logFile.error("Unable to run GAEMR formatting of assembly for QC. Raising exception ..." % e)
            raise RuntimeError

    def run_gaemr_qc(self, read_list, working_dir, options='--force --analyze_rna', cores=1):
        try:
            assert (os.path.isfile(self.assembly_fasta) and os.path.isfile(self.contig_fasta) and os.path.isfile(self.assembly_graph))
            workspace = os.path.abspath(working_dir) + '/'

            gaemr_version = open(external_wrappers_dir + 'gaemr.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running GAEMR for QC, using %s" % gaemr_version)

            os.chdir(workspace)

            gaemr_cmd = ["sh", external_wrappers_dir + "gaemr.sh", options, '-t', str(cores), '-a',
                         self.assembly_graph, '-c', self.contig_fasta, '-s', self.assembly_fasta, '--read_list=' + read_list]

            self.logFile.info("Executing the command: %s" % ' '.join(gaemr_cmd))
            stdout_file = workspace + "gaemr_qc.log"
            stderr_file = workspace + "gaemr_qc.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(gaemr_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            # Check GAEMR QC ran successfully
            assert(os.path.isdir(workspace + '/gaemr/'))
            result_files = ['assembly.basic_assembly_stats.table.txt', 'assembly.blast_hit_taxonomy.table.txt',
                            'assembly.contig_detail.table.txt', 'assembly.gap_analysis.table.txt',
                            'assembly.gap_ss_analysis.table.txt', 'assembly.scaffold_detail.table.txt',
                            'assembly.scaffolds.bam_phys_cvg_stats.table.txt',
                            'assembly.scaffolds.bam_seq_cvg_stats.table.txt',
                            'assembly.scaffolds.simple_bam_stats.table.txt', 'fragment.scaffolds.insert_size.table.txt']
            for f in result_files:
                fpath = workspace + '/gaemr/table/' + f
                assert(fpath)
                with open(fpath) as ofp:
                    file_line_count = len(ofp.readlines())
                    assert(file_line_count > 2)

            self.logFile.info("Results are available at %s" % workspace + '/gaemr/')
            return
        except RuntimeError as e:
            self.logFile.error("Unable to run GAEMR for QC. Raising exception ..." % e)
            raise RuntimeError

    def run_gaemr_qc_ont(self, read_list, working_dir, options='--force --analyze_rna', cores=1):
        try:
            assert (os.path.isfile(self.assembly_fasta) and os.path.isfile(self.contig_fasta) and os.path.isfile(self.assembly_graph))
            workspace = os.path.abspath(working_dir) + '/'

            gaemr_version = open(external_wrappers_dir + 'gaemr_ont.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running GAEMR for QC, using %s" % gaemr_version)

            os.chdir(workspace)

            gaemr_cmd = ["sh", external_wrappers_dir + "gaemr_ont.sh", options, '-t', str(cores), '-a',
                         self.assembly_graph, '-c', self.contig_fasta, '-s', self.assembly_fasta, '--read_list=' + read_list]

            self.logFile.info("Executing the command: %s" % ' '.join(gaemr_cmd))
            stdout_file = workspace + "gaemr_qc.log"
            stderr_file = workspace + "gaemr_qc.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(gaemr_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            # Check GAEMR QC ran successfully
            assert(os.path.isdir(workspace + '/gaemr/'))
            result_files = ['assembly.basic_assembly_stats.table.txt', 'assembly.blast_hit_taxonomy.table.txt',
                            'assembly.contig_detail.table.txt', 'assembly.gap_analysis.table.txt',
                            'assembly.gap_ss_analysis.table.txt', 'assembly.scaffold_detail.table.txt',
                            'assembly.scaffolds.bam_phys_cvg_stats.table.txt',
                            'assembly.scaffolds.bam_seq_cvg_stats.table.txt',
                            'assembly.scaffolds.simple_bam_stats.table.txt', 'fragment.scaffolds.insert_size.table.txt']
            for f in result_files:
                fpath = workspace + '/gaemr/table/' + f
                assert(fpath)
                with open(fpath) as ofp:
                    file_line_count = len(ofp.readlines())
                    assert(file_line_count > 2)

            self.logFile.info("Results are available at %s" % workspace + '/gaemr/')
            return
        except RuntimeError as e:
            self.logFile.error("Unable to run GAEMR for QC. Raising exception ..." % e)
            raise RuntimeError

    def run_pilon(self, working_dir, nanopore_bam=None, illumina_bam=None, options='', cores=1, reference_change=False):
        try:
            assert (os.path.isfile(illumina_bam) or os.path.isfile(nanopore_bam))
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isdir(working_dir))
            pilon_version = open(external_wrappers_dir + 'pilon.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running pilon for illumina and/or nanopore to refine assembly, using version %s" % pilon_version)

            result_assembly = workspace + self.sample_name + '.fasta'
            result_changes = workspace + self.sample_name + '.changes'
            pilon_cmd = ["/bin/bash", external_wrappers_dir + "pilon.sh", options, '--genome', self.assembly_fasta, '--changes',
                         '--threads', str(cores), '--output', self.sample_name, '--outdir', workspace]
            if nanopore_bam:
                pilon_cmd += ['--nanopore', nanopore_bam]
            if illumina_bam:
                pilon_cmd += ['--frags', illumina_bam]

            self.logFile.info("Executing the command: %s" % ' '.join(pilon_cmd))
            stdout_file = workspace + 'pilon.log'
            stderr_file = workspace + "pilon.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(pilon_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            assert (os.path.isfile(result_assembly))
            self.logFile.info("Results are available at %s" % result_assembly)

            if reference_change:
                self.assembly_fasta = result_assembly

            return [result_assembly, result_changes]
        except RuntimeError as e:
            self.logFile.error("Unable to run Pilon. Raising exception ..." % e)
            raise RuntimeError

    def run_nanopolish(self, nanopore_fastq, nanopore_bam, working_dir, nanopolish_variants_options='', cores=8, reference_change=True):
        try:
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isfile(nanopore_fastq) and os.path.isfile(nanopore_bam) and os.path.isdir(working_dir))

            workspace_vcf_dir = working_dir + 'vcfs/'
            workspace_log_dir = working_dir + 'logs/'

            os.system('mkdir %s' % workspace_vcf_dir)
            os.system('mkdir %s' % workspace_log_dir)
            nanopolish_version = open(external_wrappers_dir + 'nanopolish.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running nanopolish to refine assembly, using version %s" % nanopolish_version)
            result_assembly = workspace + self.sample_name + '.fasta'

            # run polishing.
            nanopolish_variants_cmd_1 = ["/bin/bash", external_wrappers_dir + "nanopolish_makerange.sh",
                                               self.assembly_fasta]

            self.logFile.info("Executing the command: %s" % ' '.join(nanopolish_variants_cmd_1))
            stdout_file_1 = workspace + 'nanopolish_variants_1.log'
            stderr_file_1 = workspace + 'nanopolish_variants_1.err'
            stdout_handle_1 = open(stdout_file_1, 'w')
            stderr_handle_1 = open(stderr_file_1, 'w')
            proc1 = subprocess.Popen(nanopolish_variants_cmd_1, stdout=stdout_handle_1, stderr=stderr_handle_1)
            proc1.wait()
            stdout_handle_1.close(); stderr_handle_1.close();

            cmds = []
            with open(stdout_file_1) as osf1:
                for line in osf1:
                    line = line.strip()
                    nanopolish_cmd = ['/bin/bash', external_wrappers_dir + "nanopolish.sh", "variants", "--consensus",
                                      nanopolish_variants_options,
                                      '-g', self.assembly_fasta, "-o",
                                      workspace_vcf_dir + 'polished.' + line.replace(':', '_') + '.vcf', '-w', line,
                                      '-r', nanopore_fastq, '-b', nanopore_bam, '-t', '2', '--min-candidate-frequency',
                                      '0.1']
                    self.logFile.info("Will be executing the command: %s" % ' '.join(nanopolish_cmd))
                    cmds.append([workspace_log_dir, line, nanopolish_cmd])

            p = multiprocessing.Pool(int(cores/2))
            try:
                p.map(run_nanopolish_variants, cmds)
            except:
                self.logFile.error("Problem with running nanopolish variants --consensus. Raising exception ...")
                p.close()
                raise RuntimeError
            else:
                p.close()

            vcflist_file = workspace + 'vcfs.list'
            vcflist_handle = open(vcflist_file, 'w')
            for f in os.listdir(workspace_vcf_dir):
                vcflist_handle.write(workspace_vcf_dir + f + '\n')
            vcflist_handle.close()

            # combine results across windows and convert to FASTA format.
            nanopolish_vcf2fasta_cmd = ["/bin/bash", external_wrappers_dir + "nanopolish.sh", "vcf2fasta",
                                        "--skip-checks", '-g', self.assembly_fasta, '-f', vcflist_file]
            self.logFile.info("-" * 70)
            self.logFile.info("Executing the command: %s" % ' '.join(nanopolish_vcf2fasta_cmd))
            stderr_file = workspace + 'nanopolish_vcf2fasta.err'
            stdout_handle = open(result_assembly, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(nanopolish_vcf2fasta_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            assert (os.path.isfile(result_assembly) and os.path.getsize(result_assembly) > 100)
            self.logFile.info("Polished assembly available at %s" % result_assembly)

            if reference_change:
                self.assembly_fasta = result_assembly
            return result_assembly
        except RuntimeError as e:
            self.logFile.error("Unable to run Nanopolish. Raising exception ..." % e)
            raise RuntimeError