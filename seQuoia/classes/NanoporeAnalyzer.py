import os
import sys
import subprocess
from operator import itemgetter
sys.path.append('/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]))
from other import usefulFunctions as uF

external_wrappers_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

class Nanopore():
    def __init__(self, fastq, sample_name, logObject, fast5=None, seqsum=None, barcode=None):
        self.fastq = os.path.abspath(fastq)
        self.sample_name = sample_name
        self.logFile = logObject
        self.fast5 = fast5
        self.seqsum = seqsum
        self.barcode = barcode

    def validate(self, fastq=None):
        """ Validate file is indeed a FASTQ using fqtools. """
        try:
            result = True
            if not fastq: fastq = self.fastq
            assert (fastq.endswith(".fastq.gz") or fastq.endswith(".fastq"))
            fqtools_version_info = open(external_wrappers_dir + 'fqtools.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Validating FASTQ file %s using fqtools installation %s" % (fastq, fqtools_version_info))
            validate_cmd = ["/bin/bash", external_wrappers_dir + 'fqtools.sh', 'validate', fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(validate_cmd))
            proc = subprocess.Popen(validate_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            if out.decode('utf-8').strip() != "OK" or not out.decode('utf-8').strip():
                result = False
                self.logFile.warning("Unable to validate file %s as a FASTQ" % fastq)
            else:
                self.logFile.info("Was able to successfully validate file %s as a FASTQ" % fastq)
            return result
        except RuntimeError as e:
            self.logFile.error("fqtools was not able to be run! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_nanotrim(self, working_dir, compress=True, change_reference=True):
        """ Remove adapters using PoreChop from Fastq file. """
        try:
            assert (os.path.isfile(self.fastq))
            self.logFile.info(
                "Trimming adapter sequences from reads in nanopore FASTQ file %s.\nCurrent working directory:%s" % (
                self.fastq, working_dir))

            adapter_free_fastq_file = working_dir + self.sample_name + '.fastq'

            porechop_version_info = open(external_wrappers_dir + 'porechop.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Trimming adapters with PoreChop version %s" % (porechop_version_info))
            nanoplot_cmd = ["/bin/bash", external_wrappers_dir + 'porechop.sh', '-i', self.fastq]

            self.logFile.info("Executing the command: %s" % ' '.join(nanoplot_cmd))

            stderr_file = working_dir + "porechop.err"
            stdout_handle = open(adapter_free_fastq_file, 'w');
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(nanoplot_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()
            assert (os.path.getsize(adapter_free_fastq_file) >= 100)
            assert (self.validate(adapter_free_fastq_file))
            if compress:
                try:
                    os.system('gzip ' + adapter_free_fastq_file)
                    adapter_free_fastq_file += ".gz"
                    assert (os.path.isfile(adapter_free_fastq_file))
                    self.logFile.info("Successfully compressed resulting FASTQ file.")
                except:
                    self.logFile.error("Unable to compress output FASTQ file! Raising exception ...")
                    raise RuntimeError
            if change_reference:
                self.fastq = adapter_free_fastq_file
            return adapter_free_fastq_file
        except RuntimeError as e:
            self.logFile.error("Unable to trim adapters using porechop, error message %s" % e)
            raise RuntimeError

    def run_filtlong(self, working_dir, options='--min_length 1000 --keep_percent 90 --target_bases 500000000',
                     illumina_forward=False, illumina_reverse=False,
                     reference_assembly=False, compress=True, change_reference=True):
        """ Subset Nanopore reads using filtlong for improved assembly. Options only used if neither illumina sequencing
         data nor reference fasta is provided. """
        try:
            assert (os.path.isfile(self.fastq) or os.path.isdir(working_dir))
            fitlong_version_info = open(external_wrappers_dir + 'fitlong.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Subsetting Nanopore reads using fitlong version %s" % fitlong_version_info)

            filtlong_cmd = []
            subset_fastq_results = self.working_dir + self.sample_name + '.fastq'
            if os.path.isfile(illumina_forward) and os.path.isfile(illumina_reverse):
                filtlong_cmd = ["/bin/bash", external_wrappers_dir + "fitlong.sh", '-1', illumina_forward, '-2',
                                illumina_reverse, '--trim', '--split 500', self.fastq, '>', subset_fastq_results]
            elif os.path.isfile(reference_assembly):
                filtlong_cmd = ["/bin/bash", external_wrappers_dir + "fitlong.sh", '-a', reference_assembly, '--trim',
                                '--split 500', self.fastq, '>', subset_fastq_results]
            else:
                filtlong_cmd = ["/bin/bash", external_wrappers_dir + "fitlong.sh", options, self.fastq, '>',
                                subset_fastq_results]

            self.logFile.info("Executing the command: %s" % ' '.join(filtlong_cmd))
            stdout_file = working_dir + "filtlong.log"
            stderr_file = working_dir + "filtlong.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(filtlong_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            try:
                assert (self.validate(subset_fastq_results))
            except:
                raise RuntimeError("Output doesn't exist or is invalid!")

            if compress:
                os.system('gzip %s' % subset_fastq_results)
                subset_fastq_results += '.gz'

            assert (os.path.isfile(subset_fastq_results))

            self.logFile.info("Results are available at %s" % subset_fastq_results)

            if change_reference:
                self.fastq = subset_fastq_results

            return subset_fastq_results
        except RuntimeError as e:
            self.logFile.error("Unable to run_filtlong for Nanopore read subsetting. Raising exception ..." % e)
            raise RuntimeError

    def run_fastqfilter(self, working_dir, options='-l 3000 -b 300000000', compress=True, change_reference=True):
        try:
            assert (os.path.isfile(self.fastq) and os.path.isdir(working_dir))
            fastqfilter_version = open(external_wrappers_dir + 'fastqfilter.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Subsetting Nanopore reads using fastqfilter version %s" % fastqfilter_version)

            subset_fastq_results = os.path.abspath(working_dir) + '/' + self.sample_name + '.fastq'

            fastqfilter_cmd = ["/bin/bash", external_wrappers_dir + "fastqfilter.sh", options, '-o', subset_fastq_results,
                               self.fastq]

            self.logFile.info("Executing the command: %s" % ' '.join(fastqfilter_cmd))
            stdout_file = working_dir + "fastqfilter.log"
            stderr_file = working_dir + "fastqfilter.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(fastqfilter_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            if compress:
                os.system('gzip %s' % subset_fastq_results)
                subset_fastq_results += '.gz'

            assert (os.path.isfile(subset_fastq_results))
            self.logFile.info("Results are available at %s" % subset_fastq_results)

            if change_reference:
                self.fastq = subset_fastq_results

            return subset_fastq_results

        except RuntimeError as e:
            self.logFile.error("Unable to fastqfilter for Nanopore read subsetting. Raising exception ..." % e)
            raise RuntimeError

    def run_unicycler(self, illumina_f_fastq, illumina_r_fastq, working_dir, options="--mode normal --verbosity 2", cores=1):
        try:
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isfile(illumina_f_fastq) and os.path.isfile(illumina_r_fastq) and os.path.isfile(
                self.fastq) and os.path.isdir(workspace))
            workspace = os.path.abspath(workspace) + '/'

            unicycler_version = open(external_wrappers_dir + 'unicycler.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Generating hybrid assembly with Unicyler version %s" % unicycler_version)

            unicycler_cmd = ["/bin/bash", external_wrappers_dir + "unicycler.sh", options, '-1', illumina_f_fastq,
                             '-2', illumina_r_fastq, '-l', self.fastq, '-o', workspace, '-t', str(cores)]

            self.logFile.info("Executing the command: %s" % ' '.join(unicycler_cmd))
            stdout_file = workspace + "unicycler.log"
            stderr_file = workspace + "unicycler.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(unicycler_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            result_fasta = workspace + 'assembly.fasta'
            result_gfa = workspace + 'assembly.gfa'
            assert (os.path.isfile(result_fasta) and os.path.isfile(result_gfa))
            final_fasta = workspace + self.sample_name + '.fasta'
            os.system('mv %s %s' % (result_fasta, final_fasta))
            self.logFile.info("Results are available at %s" % workspace)
            self.logFile.info("Assembly is available at %s" % final_fasta)
            return final_fasta
        except RuntimeError as e:
            self.logFile.error("Unable to run Unicycler for Nanopore/Illumina Hybrid Bacterial Assembly. Raising exception ..." % e)
            raise RuntimeError

    def run_canu(self, working_dir, options="stopOnReadQuality=false genomeSize=5m", memory=24, cores=1):
        try:
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isfile(self.fastq) and os.path.isdir(workspace))
            canu_version = open(external_wrappers_dir + 'canu.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Generating nanopore raw assembly with Canu version %s" % canu_version)

            canu_cmd = ["/bin/bash", external_wrappers_dir + "canu.sh", '-p', self.sample_name, '-d', workspace,
                        '-nanopore-raw', self.fastq, options, 'maxMemory=' + str(int(memory*cores)),
                        'maxThreads=' + str(cores), 'minReadLength=500', 'minOverlapLength=500',
                        'useGrid=false']

            self.logFile.info("Executing the command: %s" % ' '.join(canu_cmd))
            stdout_file = workspace + "canu.log"
            stderr_file = workspace + "canu.err"
            stdout_handle = open(stdout_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(canu_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            result_fasta = workspace + self.sample_name + '.contigs.fasta'
            assert (os.path.isfile(result_fasta))
            final_fasta = workspace + self.sample_name + '.fasta'
            os.system('cp %s %s' % (result_fasta, final_fasta))
            self.logFile.info("Results are available at %s" % workspace)
            self.logFile.info("Assembly is available at %s" % final_fasta)
            return final_fasta
        except RuntimeError as e:
            self.logFile.error("Unable to run Canu for Nanopore Raw Assembly. Raising exception ..." % e)
            raise RuntimeError

    def minimap_alignment(self, assembly_fasta, working_dir, options='-ax map-ont', cores=1):
        try:
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isfile(self.fastq) and os.path.isfile(assembly_fasta))

            minimap_version = open(external_wrappers_dir + 'minimap2.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running minimap2 for nanopore to assembly alignment, using version %s" % minimap_version)

            minimap_cmd = ["/bin/bash", external_wrappers_dir + "minimap2.sh", options, '-t', str(cores), assembly_fasta, self.fastq]

            self.logFile.info("Executing the command: %s" % (' '.join(minimap_cmd)))
            sam_file = workspace + self.sample_name + '.sam'
            stderr_file = workspace + "minimap.err"
            sam_file_handle = open(sam_file, 'w')
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(minimap_cmd, stdout=sam_file_handle, stderr=stderr_handle)
            proc.wait()
            sam_file_handle.close();
            stderr_handle.close()

            assert (os.path.isfile(sam_file))
            self.logFile.info("Results are available at %s" % sam_file)
            return sam_file
        except RuntimeError as e:
            self.logFile.error("Unable to run minimap alignment. Raising exception ..." % e)
            raise RuntimeError

    def recursiveFast5DirectoryIdentifier(self, directory, searching_directories):
        """ function for getting every"""
        if directory.strip('/').split('/')[-1] == self.barcode:
            searching_directories.add(directory)
            return searching_directories

        for sd in os.listdir(directory):
            subdir = os.path.abspath(directory + sd) + '/'
            if os.path.isdir(subdir):
                if (not self.barcode or self.barcode in subdir) and not '/workspace/fail/' in subdir:
                    searching_directories.add(subdir)
                elif not '/workspace/fail/' in subdir and not uF.is_number(sd):
                    searching_directories = self.recursiveFast5DirectoryIdentifier(subdir, searching_directories)
        return searching_directories

    def index_fast5(self, working_dir, change_reference=True):
        """
        Perform Nanopolish indexing. This requires the Nanopore object to have
        """
        try:
            workspace = os.path.abspath(working_dir) + '/'
            assert (os.path.isfile(self.fastq) and os.path.isdir(self.fast5))

            nanopolish_version = open(external_wrappers_dir + 'nanopolish.txt').readlines()[0]
            self.logFile.info("-" * 70)
            self.logFile.info("Running nanopolish index using version %s" % nanopolish_version)

            symlink_fastq = workspace + self.fastq.split('/')[-1]
            os.system('ln -s %s %s' % (self.fastq, symlink_fastq))

            directories_to_search = self.recursiveFast5DirectoryIdentifier(self.fast5, set([]))

            self.logFile.info("Will be looking for FAST5 files in the following directories and subdirectories recursively:" + '\n'.join(directories_to_search))
            nanopolish_cmd = ["/bin/bash", external_wrappers_dir + "nanopolish.sh", "index"] + \
                             ['-d ' + x for x in sorted(directories_to_search)] +\
                             [symlink_fastq]

            self.logFile.info("Executing the command: %s" % ' '.join(nanopolish_cmd))

            stdout_file = workspace + "nanopolish.log"
            stderr_file = workspace + "nanopolish.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(nanopolish_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            fastq_basename = self.fastq.split('/')[-1]
            check_files = ['.index', '.index.gzi', '.index.fai', '.index.readdb']

            for cf in check_files:
                cf_path = workspace + fastq_basename + cf
                assert(os.path.isfile(cf_path))

            if change_reference:
                self.fastq = symlink_fastq
            self.logFile.info("Indexing results are available at %s" % workspace)
        except RuntimeError as e:
            self.logFile.error("Unable to run Nanopolish indexing. Raising exception ..." % e)
            raise RuntimeError