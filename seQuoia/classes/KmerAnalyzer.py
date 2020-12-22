import os
import sys
import subprocess
import shutil
import h5py

external_wrappers_dir = '/'.join(
    os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

class Kmer():
    # The Fastq object consists of just two elements, a current path to some FASTQ file and a logging object.
    def __init__(self, kmer_file, sample_name, logObject):
        self.kmer_file = kmer_file
        self.sample_name = sample_name
        self.logFile = logObject
        self.assert_hdf5()

    def assert_hdf5(self):
        self.logFile.info("Attempting to validate whether kmer file %s is in valid format." % self.kmer_file)
        try:
            hf = h5py.File(self.kmer_file, 'r')
            hf.close()
            self.logFile.info("Successfully validated file as having HDF5 format!")
        except RuntimeError as e:
            self.logFile.error("Unable to validate file as having HDF5 format! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    def run_straingst(self, working_dir, db, options=""):
        """ Run StrainGST. """
        try:
            strainge_version_info = open(external_wrappers_dir + "straingst_run.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Running StrainGST on kmer HDF5 formatted file %s using installation %s" % (self.kmer_file, strainge_version_info))
            result = working_dir + self.sample_name + '.straingst_result.tsv'

            straingst_cmd = ['/bin/bash', external_wrappers_dir + 'straingst_run.sh', options, '-o', result, db, self.kmer_file]
            self.logFile.info("Executing the command: %s" % ' '.join(straingst_cmd))
            stdout_file = working_dir + self.sample_name + ".straingst.log"; stderr_file = working_dir + self.sample_name + ".straingst.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(straingst_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()
            stderr = ' '.join(open(stderr_file).readlines())

            assert('INFO:root:Done.' in stderr)
            assert(os.path.isfile(result))
            self.logFile.info("StrainGST finished successfully!")
            return result
        except RuntimeError as e:
            self.logFile.error("StrainGST was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError

    """
    IN PROGRESS
    def create_histogram(self, working_dir, plot=True):
        # Create a histogram of the k-mer spectrum.
        try:
            strainge_version_info = open(external_wrappers_dir + "strainge.txt").readlines()[0]
            self.logFile.info('-'*70)
            self.logFile.info("Extracting k-mer spectrum histogram of kmer HDF5 formatted file %s using installation %s" % (self.kmer_file, strainge_version_info))
            result_txt = working_dir + self.sample_name + '.kmer_spectrum.histogram.txt'
            strainge_stats_cmd = ['/bin/bash', external_wrappers_dir + 'kmerseq.sh', '--histogram', result_txt]
            if plot:
                result_png = working_dir + self.sample_name + '.kmer_spectrum.histogram.png'
                strainge_stats_cmd += ['-s', result_png]
            strainge_stats_cmd += ['--kmerset', self.kmer_file]
            self.logFile.info("Executing the command: %s" % ' '.join(strainge_stats_cmd))
            stdout_file = working_dir + self.sample_name + ".strainge_stats.log"; stderr_file = working_dir + self.sample_name + ".strainge_stats.err"
            stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
            proc = subprocess.Popen(strainge_stats_cmd, stdout=stdout_handle, stderr=stderr_handle)
            proc.wait()
            stdout_handle.close(); stderr_handle.close()

            self.logFile.info("StrainGE Stats & Plot for K-mer Spectrum finished successfully!")
        except RuntimeError as e:
            self.logFile.error("StrainGE was not able to be run properly! Error message: %s. Raising exception ..." % e)
            raise RuntimeError
    """