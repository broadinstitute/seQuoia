import os
import sys
import subprocess
import shutil
import h5py

external_wrapper_dir = '/'.join(
    os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

class Alignment():
    def __init__(self, alignment_file, sample_name, logObject):
        self.alignment_file = alignment_file
        self.sample_name = sample_name
        self.logFile = logObject
        self.validate()

    """
    Because these files are extremly large, I incorporate features to allow for instant clean up rather than saving to 
    clear intermediate files with the Cleanup_task module.
    """

    def create_new_instance(self, working_dir, change_reference=True):
        """ Create new SAM/BAM instance """
        self.logFile.info("Creating new instance of Alignment file %s in working directory %s" % (self.alignment_file, working_dir))
        result = None
        suffix = self.parse_suffix()
        try:
            working_dir = os.path.abspath(working_dir) + '/'
            if not os.path.isfile(working_dir + self.alignment_file.split('/')[-1]):
                shutil.copy(self.alignment_file, working_dir)
                result = working_dir + os.path.basename(self.alignment_file)
            else:
                result = working_dir + self.sample_name + '.copy.' + suffix
                shutil.copy(self.alignment_file, result)
            self.logFile.info("Successfully created new instance of Alignment file: %s." % result)
        except RuntimeError as e:
            self.logFile.error("Unable to copy Alignment file %s to working directory %s. Error received %s" % (self.alignment_file, working_dir, e))
            raise RuntimeError
        if change_reference:
            self.alignment_file = result
        self.logFile.info("Successfully created new instance of alignment file!")
        return result

    def parse_suffix(self):
        suffix = None
        try:
            if '.sorted.bam' in self.alignment_file: suffix = '.sorted.bam'
            elif '.bam' in self.alignment_file: suffix = '.bam'
            elif '.sam' in self.alignment_file: suffix = '.sam'
            return suffix
        except:
            self.logFile.error("Unable to parse out suffix of file.")
            raise RuntimeError

    def validate(self):
        """ Simply validate alignment file exists and that it endswith an appropriate suffix. """
        self.logFile.info("Attempting to validate whether alignment file %s is in SAM/BAM/CRAM format." % (self.alignment_file))
        try:
            # currently just checking for reasonable suffices, samtools quickcheck and picards validatesamfile are not really too great
            assert(os.path.isfile(self.alignment_file))
            assert(self.alignment_file.endswith(".sam") or self.alignment_file.endswith(".bam") or self.alignment_file.endswith(".cram"))
            self.logFile.info("Successfully validated alignment file format!")
        except RuntimeError as e:
            self.logFile.error("Unable to validate file as being in SAM/BAM/CRAM format! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def compress_sam(self, working_dir, clean=False, change_reference=True):
        """ Convert SAM file to BAM file. """
        try:
            samtools_version = open(external_wrapper_dir + 'samtools.txt').readlines()[0].strip()
            assert(self.alignment_file.endswith('.sam'))
            self.logFile.info("Attempting to compress SAM file %s to BAM format using %s." % (self.alignment_file, samtools_version))
            result = working_dir + self.sample_name + '.bam'
            stdout_file = working_dir + self.sample_name + ".samtools_view.log"; stderr_file = working_dir + self.sample_name + ".samtools_view.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle= open(stderr_file, 'w')
            samtools_view_cmd = ["/bin/sh", external_wrapper_dir + "samtools.sh", 'view', '-S', '-b', '-o', result, self.alignment_file]
            self.logFile.info("Executing the command: %s" % ' '.join(samtools_view_cmd))
            proc = subprocess.Popen(samtools_view_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines()).strip()
            if clean:
                os.system('rm -f %s' % self.alignment_file)
                self.alignment_file = result
            elif change_reference:
                self.alignment_file = result
            self.logFile.info("Successfully compressed sam to bam!")
            return result
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to compress sam to bam! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def sort_bam(self, working_dir, clean=False, change_reference=True, cores=1):
        """ Sort BAM file. """
        try:
            samtools_version = open(external_wrapper_dir + 'samtools.txt').readlines()[0].strip()
            assert(self.alignment_file.endswith('.bam'))
            self.logFile.info("Attempting to sort BAM file %s using %s." % (self.alignment_file, samtools_version))
            result = working_dir + self.sample_name + '.sorted.bam'
            stdout_file = working_dir + self.sample_name + ".samtools_sort.log"; stderr_file = working_dir + self.sample_name + ".samtools_sort.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle= open(stderr_file, 'w')
            samtools_sort_cmd = ["/bin/sh", external_wrapper_dir + "samtools.sh", 'sort', '-@', str(cores), '-o', result, self.alignment_file]
            self.logFile.info("Executing the command: %s" % ' '.join(samtools_sort_cmd))
            proc = subprocess.Popen(samtools_sort_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            assert(os.path.getsize(result) > 0)
            if clean:
                os.system('rm -f %s' % self.alignment_file)
                self.alignment_file = result
            elif change_reference:
                self.alignment_file = result
            self.logFile.info("Successfully sorted bam!")
            return result
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to sort bam! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def index_bam(self, working_dir, make_copy=False):
        """ Index BAM file. """
        try:
            samtools_version = open(external_wrapper_dir + 'samtools.txt').readlines()[0].strip()
            assert(self.alignment_file.endswith('.bam'))
            if make_copy:
                self.create_new_instance(working_dir)
            self.logFile.info("Attempting to index BAM file %s using %s." % (self.alignment_file, samtools_version))
            stdout_file = working_dir + self.sample_name + ".samtools_index.log"; stderr_file = working_dir + self.sample_name + ".samtools_index.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            samtools_index_cmd = ["/bin/sh", external_wrapper_dir + "samtools.sh", 'index', self.alignment_file]
            self.logFile.info("Executing the command: %s" % ' '.join(samtools_index_cmd))
            proc = subprocess.Popen(samtools_index_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines()).strip()
            self.logFile.info("Successfully indexed bam!")
            return True
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to index bam. Error message: %s. Raising exception ..." % e)

    def mark_dups(self, working_dir, change_reference=True, clean=False):
        """ Mark duplicates using Picard"""
        try:
            picard_version = open(external_wrapper_dir + 'picard.txt').readlines()[0].strip()
            assert (self.alignment_file.endswith('.bam'))
            self.logFile.info("Attempting to mark duplicates in BAM file %s, using %s." % (self.alignment_file, picard_version))
            result = working_dir + self.sample_name + '.dedup.bam'
            metrics_file = working_dir + self.sample_name + '.dedup.metrics.txt'
            stdout_file = working_dir + self.sample_name + ".picard_markdups.log"; stderr_file = working_dir + self.sample_name + ".picard_markdups.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            picard_markdups_cmd = ["/bin/sh", external_wrapper_dir + "picard.sh", 'MarkDuplicates', 'INPUT=' + self.alignment_file, 'OUTPUT=' + result, 'METRICS_FILE=' + metrics_file]
            self.logFile.info("Executing the command: %s" % ' '.join(picard_markdups_cmd))
            proc = subprocess.Popen(picard_markdups_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines()).strip()
            assert (not 'ERROR' in stderr)
            if clean:
                os.system('rm -f %s' % self.alignment_file)
                self.alignment_file = result
            elif change_reference:
                self.alignment_file = result
            self.logFile.info("Successfully marked duplicates!")
            return True
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to mark duplicates. Error message: %s. Raising exception ..." % e)

    def extract_fastqs(self, working_dir, include_non_pf_reads=False, compress=True):
        """ Extract FASTQ(s) from BAM file using Picard SamToFastq """
        try:
            samtools_version = open(external_wrapper_dir + 'samtools.txt').readlines()[0].strip()
            assert(self.alignment_file.endswith(".bam"))
            self.logFile.info("Attempting to extract FASTQ(s) from BAM file %s, using %s." % (self.alignment_file, samtools_version))
            result_r1 = working_dir + self.sample_name + "_R1.fastq"
            result_r2 = working_dir + self.sample_name + "_R2.fastq"
            stdout_file = working_dir + self.sample_name + '.extract_fastqs.log'; stderr_file = working_dir  + self.sample_name + '.extract_fastqs.err'
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            extract_cmd = ["/bin/sh", external_wrapper_dir + "picard.sh", "SamToFastq", "INTERLEAVE=true", "INCLUDE_NON_PF_READS=" + str(include_non_pf_reads).lower(), "I=" + self.alignment_file, "FASTQ=" + result_r1]
            self.logFile.info("Executing the command: %s" % ' '.join(extract_cmd))
            proc = subprocess.Popen(extract_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            assert(os.path.isfile(result_r1) and os.path.getsize(result_r1) > 0)

            # attempt to extract paired end read 2 files:
            split_cmd = ["grep", "'^@.*/2$'", "-A", "3", "--no-group-separator", result_r1, '>', result_r2]
            self.logFile.info("Executing the command: %s" % ' '.join(split_cmd))
            os.system(' '.join(split_cmd))

            if not os.path.isfile(result_r2) or os.path.getsize(result_r2) == 0:
                os.system('rm -f %s' % result_r2)
            else:
                tmp_fastq = working_dir + 'tmp.fq'
                split_cmd = ["grep", "'^@.*/1$'", "-A", "3", "--no-group-separator", result_r1, '>', tmp_fastq]
                self.logFile.info("Executing the command: %s" % ' '.join(split_cmd))
                os.system(' '.join(split_cmd))
                assert(os.path.isfile(tmp_fastq) and os.path.getsize(tmp_fastq) > 0)
                os.system('mv %s %s' % (tmp_fastq, result_r1))

            fastqs = [result_r1, result_r2]
            if compress:
                os.system('gzip %s' % result_r1)
                if os.path.isfile(result_r2): os.system('gzip %s' % result_r2)
                fastqs[0] = fastqs[0] + '.gz'; fastqs[1] = fastqs[1] + '.gz'

            fastqs = [fq for fq in fastqs if os.path.isfile(fq)]
            return fastqs
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to extract fastq files from Bam using samtools. Error message: %s. Raising exception ..." % e)

    def run_pilon(self, working_dir, reference, options='--variant'):
        """ Run Pilon variant calling. Assumes alignment is Illumina based. """
        try:
            pilon_version = open(external_wrapper_dir + 'pilon.txt').readlines()[0].strip()
            assert (self.alignment_file.endswith('.bam'))
            self.logFile.info(
                "Attempting to call variants on the BAM file %s, using %s." % (self.alignment_file, pilon_version))
            output_directory = os.path.abspath(working_dir) + '/results/'
            os.system('mkdir %s' % output_directory)
            stdout_file = working_dir + self.sample_name  + ".pilon.log"; stderr_file = working_dir + self.sample_name  + ".pilon.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            assert(os.path.isfile(reference))
            pilon_cmd = ["/bin/sh", external_wrapper_dir + "pilon.sh", options, '--genome', reference,
                                   '--frags', self.alignment_file, '--outdir', output_directory]
            self.logFile.info("Executing the command: %s" % ' '.join(pilon_cmd))
            proc = subprocess.Popen(pilon_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close();
            stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines()).strip()
            assert (not 'ERROR' in stderr)
            self.logFile.info("Successfully called variants!")
            return True
        except RuntimeError as e:
            self.logFile.error("ERROR: Unable to call variants. Error message: %s. Raising exception ..." % e)