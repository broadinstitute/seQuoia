import os
import sys
import subprocess
from collections import defaultdict
import logging
from seQuoia.classes.hUGE import hUGE
import numpy as np
from Bio import SeqIO
import glob
import time
import subprocess

true_set = set(['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh'])


def check_completion(checkpoint_file):
    try:
        assert (os.path.isfile(checkpoint_file))
    except:
        time.sleep(60)
        if not os.path.isfile(checkpoint_file):
            raise RuntimeWarning();
            sys.exit(1)

def compare_assemblies(assembly_1, assembly_2):
    assemblies_same = True

    seqs_1 = set([])
    seqs_2 = set([])

    try:
        with open(assembly_1) as oa1:
            for rec in SeqIO.parse(oa1, 'fasta'):
                seqs_1.add(str(rec.seq))

        with open(assembly_2) as oa2:
            for rec in SeqIO.parse(oa2, 'fasta'):
                seqs_2.add(str(rec.seq))
    except:
        sys.stderr.write('error with reading in assemblies as fastas. exiting now.')
        raise RuntimeError()

    diff = seqs_1.symmetric_difference(seqs_2)
    if diff:
        assemblies_same = False

    return assemblies_same

def checkIfFileExists(check_file):
    try:
        with open(check_file) as file:
            pass
        return True
    except IOError as e:
        return False

def depth(d):
    level = 1
    try:
        for i, sub_item in enumerate(d):
            if isinstance(d[sub_item], list) and len(d[sub_item]) > 0 and isinstance(d[sub_item][0], dict): level +=1
            break
        return level
    except:
        sys.stderr.write('Error with getting depth of config dictionary. Exiting now...')
        raise RuntimeError()

def createLoggerObject(log_file):
    logger = logging.getLogger('task_logger')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    #logger.handlers[0].stream = sys.stderr
    return logger

def closeLoggerObject(logObject):
    handlers = logObject.handlers[:]
    for handler in handlers:
        handler.close()
        logObject.removeHandler(handler)

def wrapAppropriately(string):
    if string:
        return "'\"" + string + "\"'"
    else:
        return "'\"\"'"

def setupDirectory(parent_dir, subdir_name, panic_if_exists=True):
    try:
        assert (os.path.isdir(parent_dir))
    except:
        sys.stderr.write("Parent repository/directory %s does not exist! Exiting now\n" % parent_dir)
        raise RuntimeError
    resdir = os.path.abspath(parent_dir) + '/' + subdir_name.rstrip('/') + '/'
    try:
        os.makedirs(resdir)
    except:
        if panic_if_exists:
            sys.stderr.write("Directory %s already exists. Exiting now.\n" % resdir)
            raise RuntimeError
    return resdir

def runSingleJob(cmd, name, sample, workspace, threads=1, mem=5, timelimit="02:00:00", cluster="Local", project=None):
    if type(cmd) is list:
        cmd = ' '.join(cmd)

    log_file = workspace + 'j' + sample + '_' + name + '.log'
    pbs_script = workspace + 'j' + sample + '_' + name + '.pbs'

    controller = hUGE(cmd, pbs_script, log_file, memory=mem, cpus=threads, throttle=None, platform=cluster, runtime=timelimit, project=project)
    controller.run_cmd()

def runSingleJobBackground(cmd, name, sample, workspace, threads=1, mem=5, timelimit="02:00:00", cluster="Local", project=None):
    if type(cmd) is list:
        cmd = ' '.join(cmd)

    log_file = workspace + 'j' + sample + '_' + name + '.log'
    pbs_script = workspace + 'j' + sample + '_' + name + '.pbs'

    controller = hUGE(cmd, pbs_script, log_file, memory=mem, cpus=threads, throttle=None, platform=cluster, runtime=timelimit, project=project)
    job_id = controller.run_cmd_no_wait()
    return(job_id)

def extractResultingAlignment(module, sample_workspace):
    resultingFile = None
    if (module == 'reference_aligner'):
        module_dir = sample_workspace + 'ReferenceAlignment/'
        list_of_files = glob.glob(module_dir + '*')
        for f in sorted(list_of_files, key=os.path.getctime):
            if f.endswith('.bam') or f.endswith('.sam'):
                resultingFile = f
    else:
        raise RuntimeWarning("module is not appropriate.")
    return resultingFile

def fastqTofasta(fastqs, fasta):
    finished = False
    try:
        for fq in fastqs:
            assert(os.path.isfile(fq))
            os.system("zcat %s | sed -n '1~4s/^@/>/p;2~4p' >> %s" % (fq, fasta))
        finished = True
    except:
        raise RuntimeWarning("unable to properly merge and convert FASTQs to FASTA")
    return finished

def extractResultingFastqs(module, sample_workspace, singleFastq, forwardFastq, reverseFastq):
    resultingFiles = ["", "", ""]
    sample_name = sample_workspace.split('/')[-2]
    if (module.strip().lower() == "gp_process"):
        singleFastq = sample_workspace + "ProcessGPDirectory/" + sample_name + '_R1.fastq'
        forwardFastq = sample_workspace + "ProcessGPDirectory/" + sample_name + '_R1.fastq'
        reverseFastq = sample_workspace + "ProcessGPDirectory/" + sample_name + "_R2.fastq"
    elif (module.strip().lower() == 'symlink'):
        sample_name = sample_workspace.split('/')[-2]
        singleFastq = sample_workspace + "Symlink_Input/" + sample_name + '_R1.fastq'
        forwardFastq = sample_workspace + "Symlink_Input/" + sample_name + '_R1.fastq'
        reverseFastq = sample_workspace + "Symlink_Input/" + sample_name + "_R2.fastq"
    elif (module.strip().lower() == "adapter_trim"):
        singleFastq = sample_workspace + "AdapterTrim/" +sample_name + '_R1.adapter-trim.fastq'
        forwardFastq = sample_workspace + "AdapterTrim/" + sample_name + '_R1.adapter-trim.fastq'
        reverseFastq = sample_workspace + "AdapterTrim/" + sample_name + '_R2.adapter-trim.fastq'
    elif (module.strip().lower() == "quality_trim"):
        singleFastq = sample_workspace + "QualityTrim/" + sample_name + '_R1.quality-trim.fastq'
        forwardFastq = sample_workspace + "QualityTrim/" + sample_name + '_R1.quality-trim.fastq'
        reverseFastq = sample_workspace + "QualityTrim/" + sample_name + '_R2.quality-trim.fastq'
    elif (module.strip().lower() == "sortmerna"):
        singleFastq = sample_workspace + "SortMeRNA/" + sample_name + "_R1.rna-removed.fastq"
        forwardFastq = sample_workspace + "SortMeRNA/" + sample_name + "_R1.rna-removed.fastq"
        reverseFastq = sample_workspace + "SortMeRNA/" + sample_name + "_R2.rna-removed.fastq"
    elif (module.strip().lower().startswith("subsample")):
        singleFastq = sample_workspace + "Subsample/" + sample_name + '_R1.subsampled.fastq'
        forwardFastq = sample_workspace + "Subsample/" + sample_name + '_R1.subsampled.fastq'
        reverseFastq = sample_workspace + "Subsample/" + sample_name + '_R2.subsampled.fastq'
    elif (module.strip().lower() == "bayeshammer"):
        singleFastq = sample_workspace + "BayesHammer/" + sample_name + '_R1.bayeshammer.fastq'
        forwardFastq = sample_workspace + "BayesHammer/" + sample_name + '_R1.bayeshammer.fastq'
        reverseFastq = sample_workspace + "BayesHammer/" + sample_name + '_R2.bayeshammer.fastq'
    elif (module.strip().lower() == 'kneaddata'):
        singleFastq = sample_workspace + "KneadData/" + sample_name + '_R1.kneaddata.fastq'
        forwardFastq = sample_workspace + "KneadData/" + sample_name + '_R1.kneaddata.fastq'
        reverseFastq = sample_workspace + "KneadData/" + sample_name + '_R2.kneaddata.fastq'
    else:
        raise RuntimeWarning("module is not appropriate.")
        return [singleFastq, forwardFastq, reverseFastq]
    if os.path.isfile(singleFastq): resultingFiles[0] = singleFastq
    if os.path.isfile(forwardFastq): resultingFiles[1] = forwardFastq
    if os.path.isfile(reverseFastq): resultingFiles[2] = reverseFastq
    singleFastq += ".gz"; forwardFastq += ".gz"; reverseFastq += ".gz"
    if os.path.isfile(singleFastq): resultingFiles[0] = singleFastq
    if os.path.isfile(forwardFastq): resultingFiles[1] = forwardFastq
    if os.path.isfile(reverseFastq): resultingFiles[2] = reverseFastq
    if not (os.path.isfile(resultingFiles[2]) and os.path.isfile(resultingFiles[1])):
        resultingFiles[1] = ""
        resultingFiles[2] = ""
    else:
        resultingFiles[0] = ""
    return resultingFiles

def cleanUp(directory):
    try:

        if os.path.isdir(directory):
            os.system('rm -rf ' + directory)
    except:
        sys.stderr.write("Unable to remove directory")
        raise RuntimeError

def is_number(x):
    try: float(x); return True
    except: return False

def kraken_report_parser(kraken_report_file):
    useful_taxonomic_levels = {'D': 'Domain', 'P': 'Phylum'}
    unclassified_reads = 0
    classified_reads = 0
    human_reads = 0
    ecoli_reads = 0
    taxonomic_counts = {}
    with open(kraken_report_file) as okrf:
        for line in okrf:
            line = line.rstrip('\n').strip()
            ls = line.split()
            total_read_count = int(ls[1])
            taxonomy_level = ls[3]
            taxonomy = ' '.join([x for x in ls[5:] if x])
            if taxonomy_level in useful_taxonomic_levels:
                taxonomy_level = useful_taxonomic_levels[taxonomy_level]
                if not taxonomy_level in taxonomic_counts:
                    taxonomic_counts[taxonomy_level] = defaultdict(int)
                taxonomic_counts[taxonomy_level][taxonomy] = total_read_count
            if taxonomy == 'unclassified': unclassified_reads = total_read_count
            elif taxonomy == 'root': classified_reads = total_read_count
            elif taxonomy == 'Homo sapiens': human_reads = total_read_count
            elif taxonomy == 'Escherichia coli': ecoli_reads = total_read_count
    return [unclassified_reads, classified_reads, human_reads, ecoli_reads, taxonomic_counts]

""" 
Following functions are taken from stackoverflow forum:
https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
"""

def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    import warnings
    warnings.filterwarnings('ignore')
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh

def percentile_based_outlier(data, threshold=95):
    diff = (100 - threshold) / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])
    return (data < minval) | (data > maxval)

def parse_best_straingst_reference(straingst_result_file):
    best_reference_fasta = None
    try:
        with open(straingst_result_file) as osrf:
            for i, line in enumerate(osrf):
                if i == 3:
                    line = line.rstrip('\n')
                    ls = line.split('\t')
                    assert(ls[0] == '0')
                    best_reference_fasta = ls[1]
    except:
        sys.stderr.write('unable to parse StrainGST results.')
    return best_reference_fasta

def fasta_validation(file_path, check_bwa_index=False):
    try:
        with open(file_path, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
        if check_bwa_index:
            bwa_index_suffices = ['.sa', '.amb', '.ann', '.pac', '.bwt']
            for bis in bwa_index_suffices:
                assert(os.path.isfile(file_path + bis))
        return True
    except:
        sys.stderr.write('unable to validate FASTA file exists/in correct format or confirm that it contains the needed BWWA index.')
        return False

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh