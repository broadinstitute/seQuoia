import os
import subprocess
from seQuoia.other import usefulFunctions as uF

external_wrappers_dir = '/'.join(
os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['external_wrappers']) + '/'

"""
To keep lightweight these functions assume the input has already been vetted.
"""

def recursive_search(root_dir, barcode, concat_fastq_file, logFile):
    if '/workspace/fail/' in root_dir or 'fast5' in root_dir: return
    if 'barcode' in root_dir and not barcode in root_dir: return
    if barcode in root_dir and not '/workspace/pass/' in root_dir: return
    root_dir = os.path.abspath(root_dir) + '/'
    for sub in os.listdir(root_dir):
        if os.path.isdir(root_dir + sub):
            sub_dir = root_dir + sub + '/'
            recursive_search(sub_dir, barcode, concat_fastq_file, logFile)
        elif os.path.isfile(root_dir + sub):
            sub_file = root_dir + sub
            if sub.endswith('.fastq') or sub.endswith('.fq'):
                logFile.info('Including file %s' % sub_file)
                os.system('cat %s >> %s' % (sub_file, concat_fastq_file))
            elif sub.endswith('.fastq.gz') or sub.endswith('.fq.gz'):
                logFile.info('Including file %s' % sub_file)
                os.system('zcat %s >> %s' % (sub_file, concat_fastq_file))
    return

def concat_fastqs(fastq_dir, working_dir, sample_name, logFile, compress=True, barcode=None):
    """ Concatenate FASTQs within directory supplied in a directory. Will search recursively through
        subdirectories when possible. If barcode is provided, "/barcode_id/" will need to be in file paths for
        fastq files to be considered as belonging to sample. Additionally '/workspace/fail/' cannot be in the path.
    """
    try:
        assert(os.path.isdir(fastq_dir))
        fastq_dir = os.path.abspath(fastq_dir) + '/'
        logFile.info("Creating concatenated FASTQ files from those located in directory %s.\nCurrent working directory:%s" % (fastq_dir, working_dir))

        concat_fastq_file = working_dir + sample_name + '.fastq'
        recursive_search(fastq_dir, barcode, concat_fastq_file, logFile)

        assert(validate(concat_fastq_file, logFile))

        if compress:
            try:
                os.system('gzip ' + concat_fastq_file)
                concat_fastq_file += ".gz"
                assert(os.path.isfile(concat_fastq_file))
                logFile.info("Successfully compressed resulting FASTQ file.")
            except:
                logFile.error("Unable to compress output FASTQ file! Raising exception ...")
                raise RuntimeError
        return concat_fastq_file
    except RuntimeError as e:
        logFile.error("Unable to concatenate FASTQ, from files located in directory at: %s\nError received:%s" % (fastq_dir, e))
        raise RuntimeError

def validate(fastq, logFile):
    """ Validate file is indeed a FASTQ using fqtools. """
    try:
        result = True
        assert (fastq.endswith(".fastq.gz") or fastq.endswith(".fastq"))
        fqtools_version_info = open(external_wrappers_dir + 'fqtools.txt').readlines()[0]
        logFile.info("-" * 70)
        logFile.info("Validating FASTQ file %s using fqtools installation %s" % (fastq, fqtools_version_info))
        validate_cmd = ["/bin/bash", external_wrappers_dir + 'fqtools.sh', 'validate', fastq]
        logFile.info("Executing the command: %s" % ' '.join(validate_cmd))
        proc = subprocess.Popen(validate_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = proc.communicate()
        if out.decode('utf-8').strip() != "OK" or not out.decode('utf-8').strip():
            result = False
            logFile.warning("Unable to validate file %s as a FASTQ" % fastq)
        else:
            logFile.info("Was able to successfully validate file %s as a FASTQ" % fastq)
        return result
    except RuntimeError as e:
        logFile.error("fqtools was not able to be run! Error message: %s. Raising exception ..." % e)
        raise RuntimeError

def filter_sequence_summary(sequencing_summary, barcode, workspace, logFile):
    try:
        assert(os.path.isfile(sequencing_summary))
        logFile.info("Starting attempt to filter sequencing summary file %s for rows pertaining to barcode %s" % (sequencing_summary, barcode))
        outf = open(workspace + 'sample_sequencing_summary.txt', 'w')
        with open(sequencing_summary) as oss:
            for i, line in enumerate(oss):
                line = line.strip()
                if i == 0:
                    outf.write(line + '\n')
                elif ('\t' + barcode + '\t') in line:
                    outf.write(line + '\n')
        outf.close()
        logFile.info("Successfully filtered sequencing summary file!")
        return workspace + 'sample_sequencing_summary.txt'
    except:
        logFile.error("Was unable to filter sequencing summary text file for specified barcode.")
        raise RuntimeError

def run_nanoplot_qc(input, workspace, logFile, input_type="fastq_rich", cores=1):
    """ Run NanoPlot to generate QC summary report. """
    try:
        resdir = workspace + 'nanoplot/'
        nanoplot_version_info = open(external_wrappers_dir + 'nanoplot.txt').readlines()[0]
        logFile.info("-" * 70)
        logFile.info("Generating a QC report with Nanoplot version %s" % (nanoplot_version_info))
        logFile.info("Results will be written to %s" % (resdir))
        nanoplot_cmd = ["/bin/bash", external_wrappers_dir + 'nanoplot.sh', '--plots hex dot kde pauvre', '--loglength', '-t', str(cores), '-o', resdir, '--' + input_type, input]

        logFile.info("Executing the command: %s" % ' '.join(nanoplot_cmd))

        stdout_file = workspace + "nanoplot.log"
        stderr_file = workspace + "nanoplot.err"
        stdout_handle = open(stdout_file, 'w')
        stderr_handle = open(stderr_file, 'w')
        proc = subprocess.Popen(nanoplot_cmd, stdout=stdout_handle, stderr=stderr_handle)
        proc.wait()
        stdout_handle.close();
        stderr_handle.close()
    except RuntimeError as e:
        logFile.error("Unable to run NanoPlot for QC. Raising exception ..." % e)
        raise RuntimeError

def run_minion_qc(sequencing_summary, workspace, logFile):
    """ Run minion_qc to generate QC summary report. """
    try:
        assert(os.path.isfile(sequencing_summary))
        #resdir = workspace
        resdir = workspace + 'minion_qc/'
        minionqc_version_info = open(external_wrappers_dir + 'minion_qc.txt').readlines()[0]
        logFile.info("-" * 70)
        logFile.info("Generating a QC report with minion_qc version %s" % (minionqc_version_info))
        logFile.info("Results will be written to %s" % (resdir))
        minionqc_cmd = ["/bin/bash", external_wrappers_dir + 'minion_qc.sh',  '-i', sequencing_summary, '-o', resdir]

        logFile.info("Executing the command: %s" % ' '.join(minionqc_cmd))

        stdout_file = workspace + "minion_qc.log"
        stderr_file = workspace + "minion_qc.err"
        stdout_handle = open(stdout_file, 'w')
        stderr_handle = open(stderr_file, 'w')
        proc = subprocess.Popen(minionqc_cmd, stdout=stdout_handle, stderr=stderr_handle)
        proc.wait()
        stdout_handle.close(); stderr_handle.close()
        stdout = ' '.join(open(stdout_file).readlines())
        #os.system('mv %s %s' % (resdir + '/NanoQC', final_resdir))
    except RuntimeError as e:
        logFile.error("Unable to run minion_qc for QC. Raising exception ..." % e)
        raise RuntimeError

def run_pilon(assembly, sample_name, workspace, logFile, nanopore_bam=None, illumina_bam=None, options='', cores=1):
    try:
        assert (os.path.isfile(illumina_bam) or os.path.isfile(nanopore_bam))
        workspace = os.path.abspath(workspace) + '/'
        assert(os.path.isdir(workspace))
        pilon_version = open(external_wrappers_dir + 'pilon.txt').readlines()[0]
        logFile.info("-" * 70)
        logFile.info("Running pilon for illumina and/or nanopore to refine assembly, using version %s" % pilon_version)

        result_assembly = workspace + sample_name + '.fasta'
        result_changes = workspace + sample_name + '.changes'
        pilon_cmd = ["/bin/bash", external_wrappers_dir + "pilon.sh", options, '--genome', assembly, '--changes', '--threads', str(cores), '--output', sample_name, '--outdir', workspace]
        if nanopore_bam:
            pilon_cmd += ['--nanopore', nanopore_bam]
        if illumina_bam:
            pilon_cmd += ['--frags', illumina_bam]

        logFile.info("Executing the command: %s" % ' '.join(pilon_cmd))
        stdout_file = workspace + 'pilon.log'
        stderr_file = workspace + "pilon.err"
        stdout_handle = open(stdout_file, 'w')
        stderr_handle = open(stderr_file, 'w')
        proc = subprocess.Popen(pilon_cmd, stdout=stdout_handle, stderr=stderr_handle)
        proc.wait()
        stdout_handle.close(); stderr_handle.close()

        assert(os.path.isfile(result_assembly))
        logFile.info("Results are available at %s" % result_assembly)
        return [result_assembly, result_changes]
    except RuntimeError as e:
        logFile.error("Unable to run minimap alignment. Raising exception ..." % e)
        raise RuntimeError