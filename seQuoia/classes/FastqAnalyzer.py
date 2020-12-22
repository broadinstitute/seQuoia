import os
import subprocess
import shutil
from seQuoia.other import usefulFunctions as uF

external_wrappers_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

class Fastq():
    def __init__(self, fastq, sample_name, logObject):
        """ The Fastq object consists of just two elements, a current path to some FASTQ file and a logging object.
        A name can optionally be provided """
        self.fastq = fastq
        self.sample_name = sample_name
        self.logFile = logObject

    def create_symlink(self, working_dir, change_reference=True):
        working_dir = os.path.abspath(working_dir) + '/'
        symlink_fastq = working_dir + self.sample_name + '.fastq'
        if self.fastq.endswith('.gz'):
            symlink_fastq += '.gz'
        os.system('ln -s %s %s' % (self.fastq, symlink_fastq))
        if change_reference:
            self.fastq = symlink_fastq
        return symlink_fastq

    def create_new_instance(self, working_dir, compress=True, change_reference=True):
        """ Create new FASTQ instance """
        self.logFile.info("Creating new instance of FASTQ file %s in working directory %s" % (self.fastq, working_dir))
        result = None
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            if not os.path.isfile(working_dir + os.path.basename(self.fastq)):
                result = working_dir + os.path.basename(self.fastq)
                shutil.copy(self.fastq, result)
            else:
                compressed = False
                if self.fastq.endswith('.gz'): compressed = True
                result = working_dir + self.sample_name + '.copy.fastq'
                if compressed: result += '.gz'
                shutil.copy(self.fastq, result)
            self.logFile.info("Successfully created new instance of FASTQ file: %s." % result)
        except RuntimeError as e:
            self.logFile.error("Unable to copy FASTQ file %s to working directory %s. Error received %s" %
                               (self.fastq, working_dir, e))
            raise RuntimeError
        if not compress and result.endswith(".gz"):
            try:
                proc = subprocess.Popen(["gunzip", result], stdout=subprocess.PIPE)
                output = proc.stdout.read()
                if output:
                    self.logFile.error("gunzip reported back the following message:\n%s" % output)
                    raise RuntimeWarning
                else:
                    result = result[:-3]
                    self.logFile.info("Sucessfully uncompressed new FASTQ instance %s!" % result)
            except RuntimeError as e:
                self.logFile.error("Unable to uncompress FASTQ file despite it ending with *.gz %s. Error received %s" % (result, e))
                raise RuntimeError
        elif compress and result.endswith(".fastq"):
            try:
                proc = subprocess.Popen(['gzip', result], stdout=subprocess.PIPE)
                output = proc.stdout.read()
                if output:
                    self.logFile.error("gzip reported back the following message:\n%s" % output)
                    raise RuntimeWarning
                else:
                    result = result + '.gz'
                    self.logFile.info("Sucessfully compressed new FASTQ instance %s!" % result)
            except RuntimeError as e:
                self.logFile.error("Unable to compress FASDTQ file despite it ending with .fastq %s. Error received %s" % (result, e))
                raise RuntimeError
        if change_reference:
            self.fastq = result
        return result

    def downsample(self, working_dir, bases=300000000, compress=False, change_reference=True):
        """ Subsample bases from FASTQ file. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            fastqfilter_version_info = open(external_wrappers_dir + 'fastqfilter.txt').readlines()[0]
            self.logFile.info('-' * 70)
            self.logFile.info("Subsampling %d bases from the sample %s using fastqfilter version %s." % (bases, self.fastq, fastqfilter_version_info))
            result = working_dir + self.sample_name + '.subsampled.fastq'
            stdout_file = working_dir + self.sample_name + '_fastqfilter.log'; stderr_file = working_dir + self.sample_name + "_fastqfilter.err"
            fastqfilter_cmd = ['/bin/bash', external_wrappers_dir + 'fastqfilter.sh', '-b', str(bases), '-o', result, self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(fastqfilter_cmd))
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(fastqfilter_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            assert (os.path.isfile(result))
            self.logFile.info("fastqfilter ran successfully! Resulting FASTQ can be found at: %s" % result)
            if compress:
                try:
                    os.system('gzip ' + result)
                    result += ".gz"
                    assert(os.path.isfile(result))
                    self.logFile.info("Successfully compressed resulting FASTQ file.")
                except:
                    self.logFile.error("Unable to compress output FASTQ file from fastqfilter! Raising exception ...")
                    raise RuntimeError
            if change_reference:
                self.fastq = result

            return result
        except RuntimeError as e:
            self.logFile.error("fastqfilter was not able to be run! Error message %s. Raising exception ..." % e)

    def subsample(self, working_dir, reads=100000, compress=False, change_reference=True):
        """ Subsample reads from FASTQ file. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            seqtk_version_info = open(external_wrappers_dir + 'seqtk_sample.txt').readlines()[0]
            self.logFile.info('-' * 70)
            self.logFile.info("Subsampling %d reads from the sample %s using seqtk version %s." % (reads, self.fastq, seqtk_version_info))
            working_dir = os.path.abspath(working_dir) + '/'
            seqtk_cmd = ['/bin/bash', external_wrappers_dir + 'seqtk_sample.sh', '-s100', self.fastq, str(reads)]
            self.logFile.info("Executing the command: %s" % ' '.join(seqtk_cmd))
            result = working_dir + self.sample_name + '.subsampled.fastq'; stderr_file = working_dir + self.sample_name + "_seqtk.err"
            stdout_handle = open(result, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(seqtk_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            try:
                assert (os.path.isfile(result) and os.path.getsize(result) > 1000)
                self.logFile.info("seqtk ran successfully! Resulting FASTQ can be found at: %s" % result)
            except:
                self.logFile.error("seqtk failed, likely ran out of memory!")
                raise RuntimeError
            if compress:
                try:
                    os.system('gzip ' + result)
                    result += ".gz"
                    assert(os.path.isfile(result))
                    self.logFile.info("Successfully compressed resulting FASTQ file.")
                except:
                    self.logFile.error("Unable to compress output FASTQ file from seqtk! Raising exception ...")
                    raise RuntimeError
            if change_reference:
                self.fastq = result

            return result
        except RuntimeError as e:
            self.logFile.error("seqtk was not able to be run! Error message %s. Raising exception ..." % e)

    def validate(self):
        """ Validate file is indeed a FASTQ using fqtools. """
        try:
            result = True
            assert (self.fastq.endswith(".fastq.gz") or self.fastq.endswith(".fastq"))
            fqtools_version_info = open(external_wrappers_dir + 'fqtools.txt').readlines()[0]
            self.logFile.info("-"*70)
            self.logFile.info("Validating FASTQ file %s using fqtools installation %s" % (self.fastq, fqtools_version_info))
            validate_cmd = ["/bin/bash", external_wrappers_dir + 'fqtools.sh', 'validate', self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(validate_cmd))
            proc = subprocess.Popen(validate_cmd,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            if out.decode('utf-8').strip() != "OK" or not out.decode('utf-8').strip():
                result = False
                self.logFile.warning("Unable to validate file %s as a FASTQ" % self.fastq)
            else:
                self.logFile.info("Was able to successfully validate file %s as a FASTQ" % self.fastq)
            return result
        except RuntimeError as e:
            self.logFile.error("fqtools was not able to be run! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_qc(self, working_dir, cores=1, options="--quiet"):
        """ Run FastQC on file. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            fastqc_version_info = open(external_wrappers_dir + "fastqc.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running FastQC on FASTQ file %s using installation %s" % (self.fastq, fastqc_version_info))
            fastqc_cmd = ['/bin/bash', external_wrappers_dir + 'fastqc.sh', options, '-o', working_dir, '-t', str(cores),
                 self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(fastqc_cmd))
            stdout_file = working_dir + self.sample_name + "_fastqc.log"; stderr_file = working_dir + self.sample_name + "_fastqc.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(fastqc_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            original_result_zip = working_dir + self.fastq.split('/')[-1].replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '').replace('.fq', '') + '_fastqc.zip'
            original_result_html = working_dir + self.fastq.split('/')[-1].replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '').replace('.fq', '') + '_fastqc.html'
            result_zip = working_dir + self.sample_name + '_fastqc.zip'
            result_html = working_dir + self.sample_name + '_fastqc.html'
            assert(os.path.isfile(original_result_zip) and os.path.isfile(original_result_html))
            os.system('mv %s %s' % (original_result_zip, result_zip))
            os.system('mv %s %s' % (original_result_html, result_html))
            self.logFile.info("FastQC ran successfully! The results can be found at: %s" % result_zip)
            return result_zip
        except RuntimeError as e:
            self.logFile.error("FastQC was not able to be run! Error message %s. Raising exception ..." % e)
            raise RuntimeError

    def trim_galore_adapter_trim(self, working_dir, options="", change_reference=False, compress=True):
        """ Run trim_galore adapter trimming. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            trim_galore_version_info = open(external_wrappers_dir + "trim_galore.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running trim_galore on FASTQ file %s using installation %s" % (self.fastq, trim_galore_version_info))
            result = working_dir + self.sample_name + '.adapter-trim.fastq.gz'
            trim_galore_cmd = ['/bin/bash', external_wrappers_dir + 'trim_galore.sh', options]
            if compress:
                trim_galore_cmd += ['--gzip']
                if not self.fastq.endswith('.gz'):  result += '.gz'
            trim_galore_cmd +=  ['-o', working_dir, self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(trim_galore_cmd))
            stdout_file = working_dir + self.sample_name + "_trimgalore.log"; stderr_file = working_dir + self.sample_name + "_trimgalore.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(trim_galore_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            os.system('mv %s/*.fq* %s' % (working_dir, result))
            assert(os.path.isfile(result))
            self.logFile.info("trim galore finished successfully!")
            if change_reference:
                self.logFile.info("Fastq Object changed reference from %s to %s" % (self.fastq, result))
                self.fastq = result
            return result
        except RuntimeError as e:
            self.logFile.error("trim_galore was not able to be run properly! Error message %s. Raising exception ..." % e)
            raise RuntimeError

    def cutadapt_adapter_trim(self, working_dir, options="-u 25", change_reference=False):
        """ Run cutadapt adapter trimming. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            cutadapt_version_info = open(external_wrappers_dir + "cutadapt.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running cutadapt on FASTQ file %s using installation %s" % (self.fastq, cutadapt_version_info))
            result = working_dir + self.sample_name + '.adapter-trim.fastq.gz'
            cutadapt_cmd = ['/bin/bash', external_wrappers_dir + 'cutadapt.sh', options, '-o', result, self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(cutadapt_cmd))
            stdout_file = working_dir + self.sample_name + "_cutadapt.log"; stderr_file = working_dir + self.sample_name + "_cutadapt.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(cutadapt_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stdout = ' '.join(open(stdout_file).readlines())
            assert ('Finished in ' in stdout)
            assert (os.path.isfile(result))
            self.logFile.info("Cutadapt finished successfully!")
            if change_reference:
                self.logFile.info("Fastq Object changed reference from %s to %s" % (self.fastq, result))
                self.fastq = result
            return result
        except RuntimeError as e:
            self.logFile.error("cutadapt was not able to be run properly! Error message %s. Raising exception ..." % e)
            raise RuntimeError

    def quality_trim(self, working_dir, options="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25",
                        change_reference=False, cores=1):
        """ Run Trimmomatic quality trimming. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            trimmomatic_version_info = open(external_wrappers_dir + "trimmomatic.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running Trimmomatic on FASTQ file %s using installation %s" % (self.fastq, trimmomatic_version_info))
            result = working_dir + self.sample_name + '.quality-trim.fastq.gz'
            trimmomatic_cmd = ['/bin/bash', external_wrappers_dir + 'trimmomatic.sh', 'SE', '-threads', str(cores), self.fastq, result] + options.split()
            self.logFile.info("Executing the command: %s" % ' '.join(trimmomatic_cmd))
            stdout_file = working_dir + self.sample_name + ".trimmomatic.log"; stderr_file = working_dir + self.sample_name + ".trimmomatic.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(trimmomatic_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())
            assert('Exception' not in stderr and 'Unknown' not in stderr)
            assert(os.path.isfile(result))
            self.logFile.info("Trimmomatic finished successfully!")
            if change_reference:
                self.logFile.info("Fastq Object changed reference from %s to %s" % (self.fastq, result))
                self.fastq = result
            return result
        except RuntimeError as e:
            self.logFile.error("Trimmomatic was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def kmerize(self, working_dir, options="-k 23"):
        """ Run StrainGE tool kmerseq to get k-mer fingerprint. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            strainge_version_info = open(external_wrappers_dir + "straingst_kmerize.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running straingst kmerize on FASTQ file %s using installation %s" % (self.fastq, strainge_version_info))
            result = working_dir + self.sample_name + '.hdf5'


            kmerseq_cmd = ['/bin/bash', external_wrappers_dir + 'straingst_kmerize.sh', options, '-o', result, self.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(kmerseq_cmd))
            stdout_file = working_dir + self.sample_name + ".kmerize.log"; stderr_file = working_dir + self.sample_name + ".kmerize.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(kmerseq_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stdout = ' '.join(open(stdout_file).readlines())

            assert('Writing output to ' in stdout)
            assert(os.path.isfile(result))
            self.logFile.info("straingst kmerize finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error("straingst kmerize was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def filter_ribo_rna(self, working_dir, sortmerna_database_dir, change_reference=False, compress=True, paired=False, cores=1):
        """ Run SortMeRNA rRNA filtering. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            sortmerna_version_info = open(external_wrappers_dir + "sortmerna.txt").readlines()[0]
            self.logFile.info(
                "Running SortMeRNA on FASTQ file %s using installation %s" % (self.fastq, sortmerna_version_info))

            self.logFile.info("SortMeRNA ribosomal databases/references could be found at %s" % sortmerna_database_dir)

            references = [sortmerna_database_dir + 'silva-bac-16s-id90.fasta', sortmerna_database_dir + 'silva-bac-23s-id98.fasta',
            sortmerna_database_dir + 'silva-arc-16s-id95.fasta', sortmerna_database_dir + 'silva-arc-23s-id98.fasta',
            sortmerna_database_dir + 'silva-euk-18s-id95.fasta', sortmerna_database_dir + 'silva-euk-28s-id98.fasta',
            sortmerna_database_dir + 'rfam-5s-database-id98.fasta', sortmerna_database_dir + 'rfam-5.8s-database-id98.fasta']

            result = working_dir + self.sample_name + '.rna-removed'
            sortmerna_cmd = ['/bin/bash', external_wrappers_dir + 'sortmerna.sh', '--workdir', working_dir,  '--reads',
                             self.fastq] + ['--ref ' + ref for ref in references] + ['--num_alignments' , '1',
                            '--threads', str(cores), '--fastx', '--other', result, '-v']
            self.logFile.info("Executing the command: %s" % ' '.join(sortmerna_cmd))

            stdout_file = working_dir + self.sample_name + '_sortmerna.log'; stderr_file = working_dir + self.sample_name + '_sortmerna.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(sortmerna_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            result = result + '.fastq'
            assert(os.path.isfile(result))
            self.logFile.info("Successfully ran SortMeRNA!")
            if compress:
                try:
                    os.system('gzip ' + result)
                    result += ".gz"
                    assert(os.path.isfile(result))
                    self.logFile.info("Successfully compressed resulting FASTQ file.")
                except :
                    self.logFile.error("Unable to compress output FASTQ file from SortMeRNA! Raising exception ...")
                    raise RuntimeError

            if change_reference:
                self.fastq = result

            return result
        except RuntimeError as e:
            self.logFile.error("ERROR: SortMeRNA was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def align_to_reference(self, working_dir, reference, options="", cores=1):
        """ Run BWA to align reads to reference genome. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            bwa_version_info = open(external_wrappers_dir + "bwa.txt").readlines()[0]
            self.logFile.info("Running BWA on FASTQ file %s using installation %s" % (self.fastq, bwa_version_info))
            assert(uF.fasta_validation(reference, check_bwa_index=True))
            result = working_dir + self.sample_name + '.sam'

            bwa_cmd = ['/bin/bash', external_wrappers_dir + 'bwa.sh', '-t', str(cores), options, reference, self.fastq, result]
            self.logFile.info("Executing the command: %s" % ' '.join(bwa_cmd))
            stdout_file = result; stderr_file = working_dir + self.sample_name + ".bwa.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(bwa_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())

            assert('[main] Real time:' in stderr)
            assert(os.path.isfile(result))
            self.logFile.info("bwa finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error("bwa was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def error_correction(self, working_dir, change_reference=False, compress=True):
        """ Run BayesHammer Error Correction. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-' * 70)
            spades_version_info = open(external_wrappers_dir + "spades.txt").readlines()[0]

            self.logFile.info("Running FreeBayes on FASTQ file %s using installation %s" % (self.fastq, spades_version_info))

            bayeshammer_cmd = ['/bin/bash', external_wrappers_dir + 'spades.sh', '-t', '1', '-m', '40', '--only-error-correction', '-s', self.fastq, '-o', working_dir]
            self.logFile.info("Executing the command: %s" % ' '.join(bayeshammer_cmd))
            stdout_file = working_dir + self.sample_name + '.bayes_hammer.log';
            stderr_file = working_dir + self.sample_name + '.bayes_hammer.err'
            stdout_handle = open(stdout_file, 'w');
            stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(bayeshammer_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()

            result = working_dir + self.sample_name + '.bayeshammer.fastq.gz'
            corrected_reads = [working_dir + 'corrected/' + f for f in os.listdir(working_dir + 'corrected/') if f.endswith('.fastq.gz')][0]
            os.system('mv %s %s' % (corrected_reads, result))

            if not compress and result.endswith('.gz'):
                try:
                    os.system('gunzip ' + result)
                    result = result[:-3]
                except:
                    raise RuntimeError
            if change_reference:
                self.fastq = result

            self.logFile.info("BayesHammer for Error Corrections finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error(
                "bayeshammer was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_metaphlan(self, working_dir, cores=1):
        """ Run Metaphlan for Taxonomic Profiling """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            metaphlan_version_info = open(external_wrappers_dir + "metaphlan2.txt").readlines()[0]
            self.logFile.info("Running Metaphlan2 on FASTQ file %s using installation %s" % (self.fastq, metaphlan_version_info))
            metaphlan_results = working_dir + self.sample_name + ".profiled_metagenome.txt"

            stderr_file = working_dir + self.sample_name + "_metaphlan2.err"
            stdout_handle = open(metaphlan_results, 'w'); stderr_handle = open(stderr_file, 'w')
            metaphlan_cmd = ['/bin/bash', external_wrappers_dir + 'metaphlan2.sh', self.fastq, '--input_type', 'fastq', '--nproc', str(cores)]
            self.logFile.info("Executing the command: %s" % ' '.join(metaphlan_cmd))
            metaphlan_process = subprocess.Popen(metaphlan_cmd, stdout=stdout_handle, stderr=stderr_handle)
            metaphlan_process.wait()
            stdout_handle.close(); stderr_handle.close()
            assert(os.path.isfile(metaphlan_results))
            self.logFile.info("Successfully ran Metaphlan2!")
            return metaphlan_results
        except RuntimeError as e:
            self.logFile.error("ERROR: Metaphlan was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def bin_taxonomically(self, working_dir, centrifuge_index, cores=1):
        """ Run Centrifuge for Taxonomic Binning """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            centrifuge_version_info = open(external_wrappers_dir + "centrifuge.txt").readlines()[0]
            self.logFile.info("Running Centrifuge on FASTQ file %s using installation %s" % (self.fastq, centrifuge_version_info))
            old_reference = None
            if self.fastq.endswith('.gz'):
                self.logFile.warning("FASTQ files are not compressed, time to create local instances.")
                old_reference = self.fastq
                self.create_new_instance(working_dir, compress=False, change_reference=True)
            else:
                self.logFile.info("Awesome! FASTQ files are already compressed!")

            centrifuge_results = working_dir + self.sample_name + "_centrifuge_results.txt"
            centrifuge_report = working_dir + self.sample_name + "_centrifuge_report.tsv"
            centrifuge_kreport = working_dir + self.sample_name + "_centrifuge_kraken_report.txt"

            stdout_file = working_dir + self.sample_name + "_centrifuge.log"; stderr_file = working_dir + self.sample_name + "_centrifuge.err"
            stdout_handle = open(stdout_file, 'wb'); stderr_handle = open(stderr_file, 'wb')
            centrifuge_cmd = ['/bin/bash', external_wrappers_dir + 'centrifuge.sh', "-p", str(cores), "-q", "-x", centrifuge_index, "-U", self.fastq, "-S", centrifuge_results, "--report-file", centrifuge_report]
            self.logFile.info("Executing the command: %s" % ' '.join(centrifuge_cmd))
            centrifuge_process = subprocess.Popen(centrifuge_cmd, stdout=stdout_handle, stderr=stderr_handle)
            centrifuge_process.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())
            assert('Calculating abundance' in stderr)
            stderr_file = working_dir + self.sample_name + "_centrifuge_kreport.err"
            stdout_handle = open(centrifuge_kreport, 'wb'); stderr_handle = open(stderr_file, 'wb')
            centrifuge_k_cmd = ['/bin/bash', external_wrappers_dir + "centrifuge_kreport.sh", "-x", centrifuge_index, centrifuge_results]
            self.logFile.info("Executing the command: %s" % ' '.join(centrifuge_k_cmd))
            centrifuge_k_process = subprocess.Popen(centrifuge_k_cmd, stdout=stdout_handle, stderr=stderr_handle)
            centrifuge_k_process.wait()
            stdout_handle.close(); stderr_handle.close()
            assert(os.path.isfile(centrifuge_results) and os.path.isfile(centrifuge_report) and os.path.isfile(centrifuge_kreport))
            self.logFile.info("Successfully ran centrifuge and generated kraken-like report file!")
            if not old_reference:
                os.system("rm -f %s" % self.fastq)
                self.fastq = old_reference
            return [centrifuge_results, centrifuge_report, centrifuge_kreport]
        except RuntimeError as e:
            self.logFile.error("ERROR: Centrifuge was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

class FastqPaired():
    def __init__(self, fastq_frw, fastq_rev, sample_name, logObject):
        """ The FastqPaired object consists of just three elements, a current path to forward and reverse FASTQ files and a logging object.
        A name can optionally be provided """
        self.fastqFrw = Fastq(fastq_frw, sample_name + '_R1', logObject)
        self.fastqRev = Fastq(fastq_rev, sample_name + '_R2', logObject)
        self.sample_name = sample_name
        self.logFile = logObject

    def validate(self):
        """ Validate file is indeed a FASTQS. """
        self.logFile.info('*'*70)
        self.logFile.info("Validating forward FASTQ file %s and reverse FASTQ file %s." % (self.fastqFrw.fastq, self.fastqRev.fastq))
        result = True
        frw_validation = self.fastqFrw.validate()
        rev_validation = self.fastqRev.validate()
        if not frw_validation or not rev_validation:
            self.logFile.warning("Unable to validate either/both forward or/and reverse FASTQ files as being in valid format.")
            result = False
        else:
            self.logFile.info("Successfully validated both forward and reverse files as being FASTQs!")
        return result

    def create_symlink(self, working_dir, change_reference=True):
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info("Creating symlinks of FASTQ files %s and %s in working directory %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, working_dir))
        frwFastqResult = self.fastqFrw.create_symlink(working_dir, change_reference=change_reference)
        revFastqResult = self.fastqRev.create_symlink(working_dir, change_reference=change_reference)
        self.logFile.info("Successfully created new instances of FASTQ files")
        return [frwFastqResult, revFastqResult]

    def create_new_instance(self, working_dir, compress=True, change_reference=True):
        """ Create new FASTQ instances """
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info("Creating new instance of FASTQ files %s and %s in working directory %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, working_dir))
        frwFastqResult = self.fastqFrw.create_new_instance(working_dir, compress=compress, change_reference=change_reference)
        revFastqResult = self.fastqRev.create_new_instance(working_dir, compress=compress, change_reference=change_reference)
        self.logFile.info("Successfully created new instances of FASTQ files")
        return [frwFastqResult, revFastqResult]

    def subsample(self, working_dir, reads=100000, compress=True, change_reference=True):
        """ Subsample FASTQ instances """
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info('*'*70)
        self.logFile.info("Creating subsampled FASTQs for %s and %s." % (self.fastqFrw.fastq, self.fastqRev.fastq))
        frwFastqResult = self.fastqFrw.subsample(working_dir, reads=reads, compress=compress, change_reference=change_reference)
        revFastqResult = self.fastqRev.subsample(working_dir, reads=reads, compress=compress, change_reference=change_reference)
        self.logFile.info("Successfully subsampled FASTQ files")
        return [frwFastqResult, revFastqResult]

    def downsample(self, working_dir, bases=300000000, compress=True, change_reference=True):
        """ Subsample FASTQ instances """
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info('*'*70)
        self.logFile.info("Creating subsampled FASTQs for %s and %s." % (self.fastqFrw.fastq, self.fastqRev.fastq))
        frwFastqResult = self.fastqFrw.downsample(working_dir, bases=bases, compress=compress, change_reference=change_reference)
        revFastqResult = self.fastqRev.downsample(working_dir, bases=bases, compress=compress, change_reference=change_reference)
        self.logFile.info("Successfully subsampled FASTQ files")
        return [frwFastqResult, revFastqResult]

    def run_qc(self, working_dir, cores=1, options="--quiet"):
        """ Run FastQC on file. """
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info('*'*70)
        self.logFile.info("Running FastQC on FASTQ files %s and %s" % (self.fastqFrw.fastq, self.fastqRev.fastq))
        results = [None, None]
        results[0] = self.fastqFrw.run_qc(working_dir, cores=cores, options=options)
        results[1] = self.fastqRev.run_qc(working_dir, cores=cores, options=options)
        self.logFile.info("Succesfully ran FastQC on both forward and reverse read files.\nResults for the forward readset can be found at %s.\nResults for the reverse readset can be found %s.")
        return results

    def cutadapt_adapter_trim(self, working_dir, options="-u 25", change_reference=False):
        """ Run cutadapt adapter trimming. """
        working_dir = os.path.abspath(working_dir) + '/'
        self.logFile.info('*'*70)
        self.logFile.info("Running cutadapt on FASTQ files %s and %s" % (self.fastqFrw.fastq, self.fastqRev.fastq))
        frwFastqResult = self.fastqFrw.adapter_trim(working_dir, options, change_reference=change_reference)
        revFastqResult = self.fastqRev.adapter_trim(working_dir, options, change_reference=change_reference)
        self.logFile.info("Successfully ran Cutadapt on both forward and reverse read files.\nAdapter trimmed results can be found at %s and %s." % (frwFastqResult, revFastqResult))
        return [frwFastqResult, revFastqResult]

    def trim_galore_adapter_trim(self, working_dir, options="", change_reference=False, compress=True):
        """ Run trim_galore adapter trimming. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            trim_galore_version_info = open(external_wrappers_dir + "trim_galore.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running trim_galore on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, trim_galore_version_info))
            result_1 = working_dir + self.sample_name + '_R1.adapter-trim.fastq'
            result_2 = working_dir + self.sample_name + '_R2.adapter-trim.fastq'
            trim_galore_cmd = ['/bin/bash', external_wrappers_dir + 'trim_galore.sh', '--paired', options]
            if compress:
                trim_galore_cmd += ['--gzip']
                result_1 += '.gz'; result_2 += '.gz'
            trim_galore_cmd += ['-o', working_dir, self.fastqFrw.fastq, self.fastqRev.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(trim_galore_cmd))
            stdout_file = working_dir + self.sample_name + "_trimgalore.log"; stderr_file = working_dir + self.sample_name + "_trimgalore.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(trim_galore_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            os.system('mv %s/*_R1*_val_1.fq* %s' % (working_dir, result_1))
            os.system('mv %s/*_R2*_val_2.fq* %s' % (working_dir, result_2))
            assert(os.path.isfile(result_1) and os.path.isfile(result_2))
            self.logFile.info("trim galore finished successfully!")
            if change_reference:
                self.logFile.info("Fastq Objects changed reference")
                self.fastqFrw.fastq = result_1
                self.fastqRev.fastq = result_2
            return [result_1, result_2]
        except RuntimeError as e:
            self.logFile.error("trim_galore was not able to be run properly! Error message %s. Raising exception ..." % e)
            raise RuntimeError

    def quality_trim(self, working_dir, options="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25", change_reference=False, cores=1, compress=True):
        """ Run Trimmomatic quality trimming. Note, only paired reads are retained if changing reference. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('*'*70)
            trimmomatic_version_info = open(external_wrappers_dir + "trimmomatic.txt").readlines()[0]
            self.logFile.info("Running Trimmomatic on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, trimmomatic_version_info))

            outfiles = [working_dir + self.fastqFrw.fastq.split('/')[-1].split('.fastq.gz')[0].split('.fastq')[
                            0] + '.trimmomatic_P.fastq.gz',
                        working_dir + self.fastqFrw.fastq.split('/')[-1].split('.fastq.gz')[0].split('.fastq')[
                            0] + '.trimmomatic_U.fastq.gz',
                        working_dir + self.fastqRev.fastq.split('/')[-1].split('.fastq.gz')[0].split('.fastq')[
                            0] + '.trimmomatic_P.fastq.gz',
                        working_dir + self.fastqRev.fastq.split('/')[-1].split('.fastq.gz')[0].split('.fastq')[
                            0] + '.trimmomatic_U.fastq.gz']

            stdout_file = working_dir + self.sample_name + '_R1_trimmomatic.log'; stderr_file = working_dir + self.sample_name + '_R1_trimmomatic.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            trimmomatic_cmd = ['/bin/bash', external_wrappers_dir + 'trimmomatic.sh', 'PE', '-threads', str(cores), self.fastqFrw.fastq, self.fastqRev.fastq] + outfiles + options.split()
            trimmomatic_proc = subprocess.Popen(trimmomatic_cmd, stdout=stdout_handle, stderr=stderr_handle)
            trimmomatic_proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())
            assert (not 'Exception' in stderr and not 'Unknown' in stderr)
            stdout_file_rev = working_dir + self.sample_name + '_R2_trimmomatic.log'; stderr_file_rev = working_dir + self.sample_name + '_R2_trimmomatic.err'
            os.system('cp %s %s' % (stdout_file, stdout_file_rev))
            os.system('cp %s %s' % (stderr_file, stderr_file_rev))
            self.logFile.info("Trimmomatic successfully ran on paired-end FASTQ set, now attempting to merge unpaired/paired file sets ...")
            outfile_tmp = [] + outfiles
            merged_frw = working_dir + self.sample_name + '_R1.quality-trim.fastq'
            merged_rev = working_dir + self.sample_name + '_R2.quality-trim.fastq'

            if not compress:
                os.system('zcat %s > %s' % (outfiles[0], merged_frw))
                os.system('zcat %s > %s' % (outfiles[2], merged_rev))
            else:
                merged_frw += '.gz'; merged_rev += '.gz'
                os.system('mv %s %s' % (outfiles[0], merged_frw))
                os.system('mv %s %s' % (outfiles[2], merged_rev))

            outfiles = [merged_frw, merged_rev]
            assert (os.path.getsize(merged_frw) > 0 and os.path.getsize(merged_rev) > 0)
            if change_reference:
                self.fastqFrw.fastq = outfiles[0]
                self.fastqRev.fastq = outfiles[1]
            os.system('rm -f ' + ' '.join(outfile_tmp))
            self.logFile.info("Successfully ran Trimmomatic and merged paired/un-paired result files.")
            return outfiles
        except:
            self.logFile.error("ERROR: Trimmomatic was not able to be run properly! Raising exception ...\n")
            raise RuntimeError

    def align_to_reference(self, working_dir, reference, options="", cores=1):
        """ Run BWA to align reads to reference genome. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            new_reference_flag = False
            self.logFile.info('-'*70)
            bwa_version_info = open(external_wrappers_dir + "bwa.txt").readlines()[0]
            self.logFile.info("Running BWA on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, bwa_version_info))
            try:
                assert(uF.fasta_validation(reference, check_bwa_index=True))
            except:
                new_reference = working_dir + reference.split('/')[-1]
                os.system('cp -r %s %s' % (reference, new_reference))
                reference = new_reference
                new_reference_flag = True
                bwa_index_cmd = ['/bin/bash', external_wrappers_dir + 'bwa.sh', 'index', reference]
                self.logFile.info("Executing the bwa index command: %s" % ' '.join(bwa_index_cmd))
                stdout_file = working_dir + self.sample_name + '.bwa-index.out'; stderr_file = working_dir + self.sample_name + ".bwa-index.err"
                stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
                proc = subprocess.Popen(bwa_index_cmd, stdout=stdout_handle, stderr=stderr_handle)
                proc.wait()
                stdout_handle.close(); stderr_handle.close()
                try:
                    assert(uF.fasta_validation(reference, check_bwa_index=True))
                except:
                    raise RuntimeError

            result = working_dir + self.sample_name + '.sam'

            bwa_cmd = ['/bin/bash', external_wrappers_dir + 'bwa.sh', 'mem', options, '-t', str(cores), reference, self.fastqFrw.fastq, self.fastqRev.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(bwa_cmd))
            stdout_file = result; stderr_file = working_dir + self.sample_name + ".bwa.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(bwa_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())

            assert('[main] Real time:' in stderr)
            assert(os.path.isfile(result))
            if new_reference_flag:
                os.system('rm -rf %s' % reference)
            self.logFile.info("bwa finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error("bwa was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def error_correction(self, working_dir, change_reference=True, compress=True):
        """ Run BayesHammer Error Correction. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('-'*70)
            spades_version_info = open(external_wrappers_dir + "spades.txt").readlines()[0]

            self.logFile.info("Running FreeBayes on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, spades_version_info))

            bayeshammer_cmd = ['/bin/bash', external_wrappers_dir + 'spades.sh', '--only-error-correction', '-t', '1', '-m', '40', '-1', self.fastqFrw.fastq, '-2', self.fastqRev.fastq, '-o', working_dir]
            self.logFile.info("Executing the command: %s" % ' '.join(bayeshammer_cmd))
            stdout_file = working_dir + self.sample_name + '.bayes_hammer.log'; stderr_file = working_dir + self.sample_name + '.bayes_hammer.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(bayeshammer_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            result1 = working_dir + self.sample_name + '_R1.bayeshammer.fastq.gz'
            result2 = working_dir + self.sample_name + '_R2.bayeshammer.fastq.gz'
            corrected_reads1 = [working_dir + 'corrected/' + f for f in os.listdir(working_dir + 'corrected/') if f.endswith('.fastq.gz') and '_R1.' in f][0]
            corrected_reads2 = [working_dir + 'corrected/' + f for f in os.listdir(working_dir + 'corrected/') if f.endswith('.fastq.gz') and '_R2.' in f][0]
            os.system('mv %s %s' % (corrected_reads1, result1))
            os.system('mv %s %s' % (corrected_reads2, result2))

            if not compress and corrected_reads1.endswith('.gz'):
                os.system('gunzip %s' % corrected_reads1)
                os.system('gunzip %s' % corrected_reads2)
                corrected_reads1 = corrected_reads1[:-3]
                corrected_reads2 = corrected_reads2[:-3]
            if change_reference:
                self.fastqFrw.fastq = corrected_reads1
                self.fastqRev.fastq = corrected_reads2
            self.logFile.info("BayesHammer for Error Corrections finished successfully!")
            return True
        except RuntimeError as e:
            self.logFile.error("bayeshammer was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_spades(self, working_dir, read_length=150, cores=4):
        """ Run SPAdes. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('-'*70)
            spades_version_info = open(external_wrappers_dir + "spades.txt").readlines()[0]

            self.logFile.info("Running SPAdes on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, spades_version_info))

            spades_options = '-k 21,33,55,77 --careful --memory 16 --threads %d' % cores
            if read_length == 250:
                spades_options = '-k 21,33,55,77,99,127 --careful --memory 16 --threads %d' % cores
            spades_cmd = ['/bin/bash', external_wrappers_dir + 'spades.sh', spades_options, '--pe1-1', self.fastqFrw.fastq, '--pe1-2', self.fastqRev.fastq, '-o', working_dir]
            self.logFile.info("Executing the command: %s" % ' '.join(spades_cmd))
            stdout_file = working_dir + self.sample_name + '.spades.log'; stderr_file = working_dir + self.sample_name + '.spades.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(spades_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            assert(os.path.isfile(working_dir + 'scaffolds.fasta'))
            self.logFile.info("SPAdes asssembly has finished successfully!")
            return True
        except RuntimeError as e:
            self.logFile.error("SPAdes was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_unicycler(self, working_dir, cores=4):
        """ Run Unicycler wrapper of SPAdes for added perks and pilon polishing. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('-'*70)
            unicycler_version_info = open(external_wrappers_dir + "unicycler.txt").readlines()[0]

            self.logFile.info("Running short read only Unicycler on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, unicycler_version_info))

            unicycler_cmd = ['/bin/bash', external_wrappers_dir + 'unicycler.sh', '-1', self.fastqFrw.fastq, '-2', self.fastqRev.fastq, '-o', working_dir, '-t', str(cores)]
            self.logFile.info("Executing the command: %s" % ' '.join(unicycler_cmd))
            stdout_file = working_dir + self.sample_name + '.unicycler.log'; stderr_file = working_dir + self.sample_name + '.unicycler.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(unicycler_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            assert(os.path.isfile(working_dir + 'assembly.fasta'))
            self.logFile.info("Unicycler asssembly has finished successfully!")
            return True
        except RuntimeError as e:
            self.logFile.error("Unicycler was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError


    def ariba(self, working_dir, name, ariba_db):
        """ Run ARIBA for finding genetic markers from short read sequencing reads directly. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('-'*70)
            ariba_version_info = open(external_wrappers_dir + "ariba.txt").readlines()[0]

            self.logFile.info("Running ARIBA on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, ariba_version_info))

            result_dir = working_dir + name
            if not result_dir.endswith('/'): result_dir += '/'

            ariba_cmd = ['/bin/bash', external_wrappers_dir + 'ariba.sh',  'run', ariba_db, self.fastqFrw.fastq, self.fastqRev.fastq, result_dir]

            self.logFile.info("Executing the command: %s" % ' '.join(ariba_cmd))
            stdout_file = working_dir + self.sample_name + '.ariba.log'; stderr_file = working_dir + self.sample_name + '.ariba.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(ariba_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            self.logFile.info("ARIBA finished successfully!")
            return True
        except RuntimeError as e:
            self.logFile.error("ARIBA was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def shortbred_amrp(self, working_dir, shortbred_markers, cores=1):
        """ Run shortBRED for ARMP. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'

            self.logFile.info('-'*70)
            shortBRED_version_info = open(external_wrappers_dir + "shortbred_quantify.txt").readlines()[0]

            self.logFile.info("Running shortBRED on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, shortBRED_version_info))

            result_file = working_dir + self.sample_name + '_shortBRED_results.txt'
            shortbred_tmp_dir = working_dir + 'shortbred_tmpdir/'
            os.system('mkdir %s' % shortbred_tmp_dir)

            fastq_files = [self.fastqFrw.fastq, self.fastqRev.fastq]
            wgs_fasta_input = working_dir + 'wgs.fna'
            validity = uF.fastqTofasta(fastq_files, wgs_fasta_input)
            assert(validity)

            shortbred_cmd = ['/bin/bash', external_wrappers_dir + 'shortbred_quantify.sh', '--threads', str(cores), '--markers', shortbred_markers, '--wgs', wgs_fasta_input, '--result', result_file, '--tmp', shortbred_tmp_dir]

            self.logFile.info("Executing the command: %s" % ' '.join(shortbred_cmd))
            stdout_file = working_dir + self.sample_name + '.shortbred_amrp.log'; stderr_file = working_dir + self.sample_name + '.shortbred_amrp.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(shortbred_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            self.logFile.info("shortBRED for AMRP finished successfully!")
            return True
        except RuntimeError as e:
            self.logFile.error("shortbred was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def filter_ribo_rna(self, working_dir, sortmerna_database_dir, change_reference=False, compress=True, cores=1):
        """ Run SortMeRNA rRNA filtering. Note, currently paired-in is set as the default setting, thus if only one of the paired reads aligns, both are regarded as ribosomal RNA. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            sortmerna_version_info = open(external_wrappers_dir + "sortmerna.txt").readlines()[0]
            self.logFile.info(
                "Running SortMeRNA on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, sortmerna_version_info))

            self.logFile.info("SortMeRNA ribosomal databases/references could be found at %s" % sortmerna_database_dir)

            references = [sortmerna_database_dir + 'silva-bac-16s-id90.fasta', sortmerna_database_dir + 'silva-bac-23s-id98.fasta',
            sortmerna_database_dir + 'silva-arc-16s-id95.fasta', sortmerna_database_dir + 'silva-arc-23s-id98.fasta',
            sortmerna_database_dir + 'silva-euk-18s-id95.fasta', sortmerna_database_dir + 'silva-euk-28s-id98.fasta',
            sortmerna_database_dir + 'rfam-5s-database-id98.fasta', sortmerna_database_dir + 'rfam-5.8s-database-id98.fasta']

            result = working_dir + self.sample_name + '.rna-removed.'
            sortmerna_cmd = ['/bin/bash', external_wrappers_dir + 'sortmerna.sh', '--workdir', working_dir,  '--reads',
                             self.fastqFrw.fastq, '--reads', self.fastqRev.fastq] + ['--ref ' + ref for ref in references] + \
                            ['--num_alignments' , '1', '--threads', str(cores), '--fastx',
                             '--other', result, '--out2', '--paired_in', '-v']
            self.logFile.info("Executing the command: %s" % ' '.join(sortmerna_cmd))

            stdout_file = working_dir + self.sample_name + '_sortmerna.log'; stderr_file = working_dir + self.sample_name + '_sortmerna.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(sortmerna_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            result_frw = result + '_fwd.fastq'
            result_rev = result + '_rev.fastq'
            assert(os.path.isfile(result_frw) and os.path.isfile(result_rev))
            self.logFile.info("Successfully ran SortMeRNA!")

            outfile_frw = working_dir + self.sample_name + '_R1.rna-removed.fastq'
            outfile_rev = working_dir + self.sample_name + '_R2.rna-removed.fastq'

            os.system('mv ' + result_frw + ' ' + outfile_frw)
            os.system('mv ' + result_rev + ' ' + outfile_rev)

            if compress:
                try:
                    os.system('gzip ' + outfile_frw + " " + outfile_rev)
                    outfile_frw += ".gz"
                    outfile_rev += ".gz"
                    assert (os.path.isfile(outfile_frw) and os.path.isfile(outfile_rev))
                    self.logFile.info("Successfully compressed resulting FASTQ files.")
                except:
                    self.logFile.error("Unable to compress output FASTQ files from SortMeRNA! Raising exception ...")
                    raise RuntimeError

            assert(os.path.isfile(outfile_frw) and os.path.isfile(outfile_rev))
            if change_reference:
                self.fastqFrw.fastq = outfile_frw
                self.fastqRev.fastq = outfile_rev
            self.logFile.info("Successfully ran SortMeRNA and paired end read set.")

            return [outfile_frw, outfile_rev]
        except:
            self.logFile.error("ERROR: SortMeRNA was not able to be run properly! Raising exception ...\n")
            raise RuntimeError

    def kmerize(self, working_dir, options="-k 23"):
        """ Run StrainGE kmerseq to get k-mer fingerprint. """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            strainge_version_info = open(external_wrappers_dir + "straingst_kmerize.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running straingst kmerize on paired end FASTQ files %s and %s using installation %s" % (self.fastqFrw, self.fastqRev, strainge_version_info))
            result = working_dir + self.sample_name + '.hdf5'

            kmerseq_cmd = ['/bin/bash', external_wrappers_dir + 'straingst_kmerize.sh', options, '-o', result, self.fastqFrw.fastq, self.fastqRev.fastq]
            self.logFile.info("Executing the command: %s" % ' '.join(kmerseq_cmd))
            stdout_file = working_dir + self.sample_name + ".kmerize.log"; stderr_file = working_dir + self.sample_name + ".kmerize.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(kmerseq_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stdout = ' '.join(open(stdout_file).readlines())

            assert(os.path.isfile(result) and os.path.getsize(result)>1000)
            self.logFile.info("straingst kmerize finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error("straingst kmerize was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def bin_taxonomically(self, working_dir, centrifuge_index, cores=1):
        """ Run Centrifuge for Taxonomic Binning """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('*'*70)
            centrifuge_version_info = open(external_wrappers_dir + "centrifuge.txt").readlines()[0]
            self.logFile.info("Running Centrifuge on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, centrifuge_version_info))

            old_frw_reference = None; old_rev_reference = None
            if self.fastqFrw.fastq.endswith('.gz') and self.fastqRev.fastq.endswith('.gz'):
                self.logFile.warning("FASTQ files are not compressed, uh time to create local instances.")
                old_frw_reference = "%s" % self.fastqFrw.fastq
                old_rev_reference = "%s" % self.fastqRev.fastq
                self.create_new_instance(working_dir, change_reference=True, compress=False)
            else:
                self.logFile.info("Awesome! FASTQ files are already compressed!")

            centrifuge_results = working_dir + self.sample_name + "_centrifuge_results.txt"
            centrifuge_report = working_dir + self.sample_name + "_centrifuge_report.tsv"
            centrifuge_kreport = working_dir + self.sample_name + "_centrifuge_kraken_report.txt"

            stdout_file = working_dir + self.sample_name + "_centrifuge.log"; stderr_file = working_dir + self.sample_name + "_centrifuge.err"
            stdout_handle = open(stdout_file, 'wb'); stderr_handle = open(stderr_file, 'wb')
            centrifuge_cmd = ['/bin/bash', external_wrappers_dir + 'centrifuge.sh', "-p", str(cores), "-q", "-x", centrifuge_index, "-1", self.fastqFrw.fastq, "-2", self.fastqRev.fastq, "-S", centrifuge_results, "--report-file", centrifuge_report]
            self.logFile.info("Executing the command: " + ' '.join(centrifuge_cmd))
            centrifuge_process = subprocess.Popen(centrifuge_cmd, stdout=stdout_handle, stderr=stderr_handle)
            centrifuge_process.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())
            assert('Calculating abundance' in stderr)
            stderr_file = working_dir + self.sample_name + "_centrifuge_kreport.err"
            stdout_handle = open(centrifuge_kreport, 'wb'); stderr_handle = open(stderr_file, 'wb')
            centrifuge_k_cmd = ['/bin/bash', external_wrappers_dir + "centrifuge_kreport.sh", "-x", centrifuge_index, centrifuge_results]
            self.logFile.info("Executing the command: " + ' '.join(centrifuge_k_cmd))
            centrifuge_k_process = subprocess.Popen(centrifuge_k_cmd, stdout=stdout_handle, stderr=stderr_handle)
            centrifuge_k_process.wait()
            stdout_handle.close(); stderr_handle.close()
            assert(os.path.isfile(centrifuge_results) and os.path.isfile(centrifuge_report) and os.path.isfile(centrifuge_kreport) and os.path.getsize(centrifuge_kreport) > 100 and os.path.getsize(centrifuge_report) > 100)
            self.logFile.info("Successfully ran centrifuge and generated kraken-like report file!")
            if old_rev_reference and old_frw_reference:
                os.system('rm -f %s %s' % (self.fastqFrw.fastq, self.fastqRev.fastq))
                self.fastqFrw.fastq = old_frw_reference
                self.fastqRev.fastq = old_rev_reference
            return [centrifuge_results, centrifuge_report, centrifuge_kreport]
        except RuntimeError as e:
            self.logFile.error("ERROR: Centrifuge was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_kneaddata(self, working_dir, options="", change_reference=False, compress=True, cores=1):
            """ Run Kneaddadta for metagenomic quality control. """
            try:
                working_dir = os.path.abspath(working_dir) + '/'
                self.logFile.info('-'*70)
                kneaddata_version_info = open(external_wrappers_dir + "kneaddata.txt").readlines()[0]
                self.logFile.info("Running KneadData on FASTQ files %s and %s using installation %s" % (self.fastqFrw.fastq, self.fastqRev.fastq, kneaddata_version_info))
                kneaddata_cmd = ['/bin/bash', external_wrappers_dir + 'kneaddata.sh', '-t', str(cores), options,
                           '--input', self.fastqFrw.fastq, '--input', self.fastqRev.fastq, '--output', working_dir, '--output-prefix', self.sample_name + "_kneaddata"]
                self.logFile.info("Executing the command: %s" % ' '.join(kneaddata_cmd))
                stdout_file = working_dir + self.sample_name + '.kneaddata.log'; stderr_file = working_dir + self.sample_name + ".kneaddata.err"
                stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
                proc = subprocess.Popen(kneaddata_cmd, stdout=stdout_handle, stderr=stderr_handle)
                proc.wait()
                stdout_handle.close(); stderr_handle.close()

                result_frw = None
                result_rev = None
                other_fastq = []
                try:
                    result_frw = [working_dir + f for f in os.listdir(working_dir) if f.endswith('paired_1.fastq')][0]
                    result_rev = [working_dir + f for f in os.listdir(working_dir) if f.endswith('paired_2.fastq')][0]
                    other_fastq += [working_dir + f for f in os.listdir(working_dir) if (not (f.endswith('paired_1.fastq') or (f.endswith('paired_2.fastq')))) and (f.endswith('.fastq') or f.endswith('.fastq.gz'))]
                except:
                    raise RuntimeError
                rename_result_frw = working_dir + self.sample_name + '_R1.kneaddata.fastq'
                rename_result_rev = working_dir + self.sample_name + '_R2.kneaddata.fastq'
                os.system('mv %s %s' % (result_frw, rename_result_frw))
                os.system('mv %s %s' % (result_rev, rename_result_rev))

                # clean up rest of fastqs in direcory
                os.system('rm -f ' + ' '.join(other_fastq))

                assert(os.path.isfile(rename_result_frw) and os.path.isfile(rename_result_rev))
                if compress:
                    try:
                        os.system('gzip ' + rename_result_frw + " " + rename_result_rev)
                        rename_result_frw += ".gz"
                        rename_result_rev += ".gz"
                        assert (os.path.isfile(rename_result_frw) and os.path.isfile(rename_result_rev))
                    except:
                        raise RuntimeError

                if change_reference:
                    self.fastqFrw.fastq = rename_result_frw
                    self.fastqRev.fastq = rename_result_rev
                self.logFile.info("kneaddata finished successfully!")
                return [rename_result_frw, rename_result_rev]
            except RuntimeError as e:
                self.logFile.error("kneaddata was not able to be run properly! Error message: %s. Raising exception ..." % e)
                raise RuntimeError

    def run_metaphlan(self, working_dir, cores=1):
        """ Run Metaphlan for Taxonomic Profiling """
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            self.logFile.info('-'*70)
            metaphlan2_version_info = open(external_wrappers_dir + "metaphlan2.txt").readlines()[0]
            self.logFile.info("Running Metaphlan2 on FASTQ file using installation %s" % (metaphlan2_version_info))
            metaphlan_results = working_dir + self.sample_name + ".profiled_metagenome.txt"

            stderr_file = working_dir + self.sample_name + "_metaphlan2.err"
            stdout_handle = open(metaphlan_results, 'w'); stderr_handle = open(stderr_file, 'w')
            metaphlan_cmd = ['/bin/bash', external_wrappers_dir + 'metaphlan2.sh', self.fastqFrw.fastq + ',' +self.fastqRev.fastq, '--bowtie2out', working_dir + self.sample_name + '.bowtie2.bz2', '--input_type', 'fastq', '--nproc', str(cores)]
            self.logFile.info("Executing the command: %s" % ' '.join(metaphlan_cmd))
            metaphlan_process = subprocess.Popen(metaphlan_cmd, stdout=stdout_handle, stderr=stderr_handle)
            metaphlan_process.wait()
            stdout_handle.close(); stderr_handle.close()
            assert(os.path.isfile(metaphlan_results))
            self.logFile.info("Successfully ran Metaphlan2!")
            return metaphlan_results
        except RuntimeError as e:
            self.logFile.error("ERROR: Metaphlan was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError
