#!/usr/local/bin/python3

"""
sheppard
Infectious Diseases and Microbiome Platform
Bacterial Genomics
01-24-2017
Rauf Salamzade
"""

import os
import sys
import argparse
import multiprocessing
from collections import defaultdict
import importlib.machinery
import subprocess
import time
sys.tracebacklimit = None
from seQuoia.other import usefulFunctions as uF

# global variables set by multiprocessing_handler
output_directory, grid_dir, seqc_dir = [None]*3
final_cluster, final_project, workflow_module = [None]*3

def workflow_process(sample_data):
    sample, ildat, npdat, metdat, config_params = sample_data
    singleFastq = ""; forwardFastq = ""; reverseFastq = ""; bamLocation = ""
    npFastq = ""; npFast5 = ""; npSeqSum = ""; npBarcode = ""
    if ildat:
        singleFastq = ildat['single']
        forwardFastq = ildat['forward']
        reverseFastq = ildat['reverse']
        bamLocation = ildat['bam-location']
    if npdat:
        npFastq = npdat['fastq']
        npFast5 = npdat['fast5']
        npSeqSum = npdat['seqsum']
        npBarcode = npdat['barcode']
    metadata = metdat
    try:
        WorkflowInstance = workflow_module.Workflow(sample, output_directory, metadata, grid_dir, seqc_dir, singleFastq=singleFastq,
                                                    forwardFastq=forwardFastq, reverseFastq=reverseFastq, bamLocation=bamLocation,
                                                    nanoporeFastq=npFastq, nanoporeFast5=npFast5, nanoporeSeqSummary=npSeqSum, nanoporeBarcode=npBarcode,
                                                    cluster=final_cluster, project=final_project, config_parameters=config_params)
        WorkflowInstance.runWorkflow()
    except RuntimeWarning as e:
        sys.stderr.write("Warning %s with running sample %s" % (e, sample) + '\n')


def read_and_store_metadata(inputfile, logObject):
    """read input data file and store meta information for each strain in
    an ugly modified json-like string for easy transferability.
    """
    metadata_warped_json = {}
    metadata_parsed = {}
    column_name = {}
    with open(inputfile) as of:
        for i, line in enumerate(of):
            line = line.rstrip('\n')
            ls = line.split('\t')
            if i == 0:
                for j, val in enumerate(ls[1:]):
                    column_name[j] = val
            else:
                sample = ls[0]
                metadata_warped_json[sample] = []
                for j, val in enumerate(ls[1:]):
                    colname = column_name[j]
                    metadata_warped_json[sample].append("QUOTES" + colname + "QUOTES:QUOTES" + val + "QUOTES")

    for s in metadata_warped_json:
        metadata_parsed[s] = metadata_parsed[s] = "{" + "-COMMA-".join(metadata_warped_json[s]) + "}"

    logObject.info("Successfully read in meta data containing input!")
    sys.stdout.write("Successfully read in meta data containing input!\n")

    return metadata_parsed

def prepare_for_submission(cluster, project, outdir, ildata_parsed, npdata_parsed, metadata_parsed, workflow, config_params_dict):
    global final_cluster; global final_project; global workflow_module; global output_directory; global grid_dir; global seqc_dir

    final_cluster = cluster
    final_project = project
    workflow_module = importlib.machinery.SourceFileLoader('workflow_module', workflow).load_module()
    output_directory = outdir
    grid_dir = outdir + 'grid_submission_files/'
    if not os.path.isdir(grid_dir):
        os.makedirs(grid_dir)
    seqc_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/seQuoia/'
    sample_data = []
    runs_to_sample = defaultdict(set)
    for sample in metadata_parsed:
        if config_params_dict and sample in config_params_dict:
            for i, param_combo in enumerate(config_params_dict[sample]):
                data = [sample]
                if ildata_parsed and sample in ildata_parsed:
                    data.append(ildata_parsed[sample])
                else:
                    data.append(None)
                if npdata_parsed and sample in npdata_parsed:
                    data.append(npdata_parsed[sample])
                else:
                    data.append(None)
                data.append(metadata_parsed[sample])
                modified_sample_name = sample
                if len(config_params_dict[sample]) > 1:
                    modified_sample_name += '_RUN-' + str(i+1)
                param_combo_naming = ([''.join(x[0].strip().replace('_', ' ').split()) + '_' + '-'.join(x[1].strip().replace('_', ' ').split()).strip('-') for x in param_combo.items() if not ('memory' in x[0].lower() or 'threads' in x[0].lower() or 'cpu' in x[0].lower())])
                data.append(param_combo)
                data[0] = modified_sample_name
                runs_to_sample[sample].add(tuple([modified_sample_name, '_'.join([x.strip() for x in param_combo_naming]).replace('---', '-').replace('--', '-')]))
                sample_data.append(tuple(data))
        else:
            data = [sample]
            if ildata_parsed and sample in ildata_parsed:
                data.append(ildata_parsed[sample])
            else:
                data.append(None)
            if npdata_parsed and sample in npdata_parsed:
                data.append(npdata_parsed[sample])
            else:
                data.append(None)
            data.append(metadata_parsed[sample])
            if not config_params_dict:
                data.append(None)
            elif uF.depth(config_params_dict) == 1:
                data.append(config_params_dict)
            else:
                data.append(None)
            sample_data.append(tuple(data))

    return sample_data, runs_to_sample

def process_nanopore_data(np_data, metadata_parsed, logObject):
    samples = set(metadata_parsed.keys())

    sample_fastq_data = {}; sample_fast5_data = {}; sample_summary_data = {}; sample_barcode_data = {}

    try:
        with open(np_data) as ond:
            for line in ond:
                line = line.strip()
                sample, albacore_dir, barcode = line.split('\t')
                assert(os.path.isdir(albacore_dir))
                albacore_dir = os.path.abspath(albacore_dir) + '/'
                sequence_summary = albacore_dir + 'sequencing_summary.txt'
                assert(os.path.isfile(sequence_summary))
                sample_summary_data[sample] = sequence_summary
                sample_barcode_data[sample] = barcode
                sample_fast5_data[sample] = albacore_dir
                sample_fastq_data[sample] = albacore_dir
    except:
        logObject.error('Issues when parsing nanopore listing file! Exiting now ...\n')
        sys.exit(1)

    nanoporedata_parsed = {}
    for s in samples:
        nanoporedata_parsed[s] = {'fastq': "", 'fast5': "", 'seqsum': "", 'barcode': ""}
        nanoporedata_parsed[s]['fastq'] = sample_fastq_data[s]
        nanoporedata_parsed[s]['fast5'] = sample_fast5_data[s]
        nanoporedata_parsed[s]['seqsum'] = sample_summary_data[s]
        nanoporedata_parsed[s]['barcode'] = sample_barcode_data[s]

    logObject.info("Successfully assigned nanopore data for each sample!")
    sys.stdout.write("Successfully assigned nanopore data for each sample!\n")

    return nanoporedata_parsed

def process_illumina_data(ildata, datatype, metadata_parsed, nanopore_present, logObject):
    samples = set(metadata_parsed.keys())
    ildata_parsed = {}

    tmp_data = defaultdict(list)
    if os.path.isdir(ildata):
        for sf in os.listdir(ildata):
            seqdat = ildata + sf
            for s in samples:
                flag = False
                if s == sf.split('.bam'): flag = True
                elif s == sf.split('.frags')[0]: frags = True
                elif s == sf.split('.1.fq')[0] or s == sf.split('.2.fq')[0]: flag = True
                elif s == sf.split('.1.fastq')[0] or s == sf.split('.2.fastq')[0]: flag = True
                elif s == sf.split('_R1.')[0] or s == sf.split('_R2.')[0]: flag = True
                elif os.path.isdir(seqdat) and sf == s: flag = True
                if flag:
                    tmp_data[s].append(seqdat)
    else:
        with open(ildata) as osd:
            for i, line in enumerate(osd):
                line = line.rstrip()
                ls = line.split('\t')
                sample = ls[0]
                data = [os.path.abspath(x) for x in ls[1:]]
                tmp_data[sample] = data

    missing_samples = False
    # structure illumina sequencing data for each sample
    for s in samples:
        ildata_parsed[s] = defaultdict(str)
        logObject.info("Sample %s has %s data: %s\n" % (s, datatype, ', '.join(tmp_data[s])))
        try:
            if not tmp_data[s]:
                logObject.warning("Unable to find illumina sequencing data for sample %s." % s)
                sys.stderr.write("Unable to find illumina sequencing data for sample %s.\n" % s)
                missing_samples = True
                continue
            if datatype == 'bam':
                assert (len(tmp_data[s]) == 1)
                ildata_parsed[s]['bam-location'] = tmp_data[s][0]
            elif datatype == 'gp-directory':
                assert (len(tmp_data[s]) == 1)
                ildata_parsed[s]['bam-location'] = tmp_data[s][0]
            elif datatype == "illumina-single":
                assert (len(tmp_data[s]) == 1)
                ildata_parsed[s]["single"] = tmp_data[s][0]
            elif datatype == "illumina-paired":
                if '_R1' in tmp_data[s][0] and '_R2' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][0]
                    ildata_parsed[s]["reverse"] = tmp_data[s][1]
                elif '_R2' in tmp_data[s][0] and '_R1' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][1]
                    ildata_parsed[s]["reverse"] = tmp_data[s][0]
                elif '.1.f' in tmp_data[s][0] and '.2.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][0]
                    ildata_parsed[s]["reverse"] = tmp_data[s][1]
                elif '.2.f' in tmp_data[s][0] and '.1.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][1]
                    ildata_parsed[s]["reverse"] = tmp_data[s][0]
                elif '.r1.f' in tmp_data[s][0] and '.r2.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][0]
                    ildata_parsed[s]["reverse"] = tmp_data[s][1]
                elif '.r2.f' in tmp_data[s][0] and '.r1.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][1]
                    ildata_parsed[s]["reverse"] = tmp_data[s][0]
                elif '_1.f' in tmp_data[s][0] and '_2.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][0]
                    ildata_parsed[s]["reverse"] = tmp_data[s][1]
                elif '_2.f' in tmp_data[s][0] and '_1.f' in tmp_data[s][1]:
                    ildata_parsed[s]["forward"] = tmp_data[s][1]
                    ildata_parsed[s]["reverse"] = tmp_data[s][0]
                else:
                    logObject.warning("Unable to find illumina sequencing data for sample %s." % s)
                    sys.stderr.write("Unable to find illumina sequencing data for sample %s.\n" % s)
        except:
            if nanopore_present:
                logObject.warning("Unable to find illumina sequencing data for some samples or trouble parsing input. Since nanopore data present, continuing on with workflow submissions")
                sys.stderr.write("Unable to find illumina sequencing data for some samples or trouble parsing input. Since nanopore data present, continuing on with workflow submissions\n")

            else:
                logObject.error("Unable to find illumina sequencing data for some samples in meta info file. Since nanopore data absent and this is an illumina only analysis, exiting now!")
                sys.exit(1)

    if missing_samples and not nanopore_present:
        logObject.error("Unable to find illumina sequencing data for some samples in meta info file. Since nanopore data absent and this is an illumina only analysis, exiting now!")
        sys.exit(1)

    logObject.info("Successfully assigned illumina data for samples where possible!")
    sys.stdout.write("Successfully assigned illumina data for samples where possible!\n")

    return ildata_parsed

def read_config_file(config, workflow, logObject):
    logObject.info("Reading configuration parameter settings and checking if variables exist in workflow!")

    config_parameters = {}
    all_workflow_variables = set([])

    try:
        with open(workflow) as ow:
            for line in ow:
                line = line.strip()
                if not line.endswith(':') and '=' in line:
                    possible_variable_name = line.split('=')[0].strip()
                    all_workflow_variables.add(possible_variable_name)

        with open(config) as oc:
            for line in oc:
                line = line.strip()
                if not line.startswith("#") and not line.startswith('//'):
                    var_name = line.split('=')[0].strip()
                    assert(var_name in all_workflow_variables)
                    var_value = '='.join(line.split('=')[1:]).strip().strip('"')
                    config_parameters[var_name] = var_value
    except:
        logObject.error("Most likely a variable was specified in the config file which does not exist in the workflow. Please remove or adapt the config. Exiting now ...")
        exit(1)

    return config_parameters

def read_sample_config_file(sample_config, workflow, logObject):
    """

    """
    logObject.info("Reading sample configuration parameter settings and checking if variables exist in workflow!")

    config_parameters = {}
    all_workflow_variables = set([])

    try:
        with open(workflow) as ow:
            for line in ow:
                line = line.strip()
                if not line.endswith(':') and '=' in line:
                    possible_variable_name = line.split('=')[0].strip()
                    all_workflow_variables.add(possible_variable_name)

        var_index_to_name = {}
        with open(sample_config) as oc:
            for i, line in enumerate(oc):
                line = line.strip()
                ls = line.split('\t')
                if i == 0:
                    for j, val in enumerate(ls[1:]):
                        assert(val in all_workflow_variables)
                        var_index_to_name[j] = val
                else:
                    sample = ls[0]
                    if not sample in config_parameters:
                        config_parameters[sample] = []
                    params = {}
                    for j, val in enumerate(ls[1:]):
                        var_name = var_index_to_name[j]
                        params[var_name] = val.strip().strip('"')
                    config_parameters[sample].append(params)
    except:
        logObject.error("Most likely a variable was specified in the config file which does not exist in the workflow. Please remove or adapt the config. Exiting now ...")
        exit(1)

    return config_parameters

def capture_external_versioning(outdir, logObject):
    """ Log info for provenance of external software versions. """
    logObject.info("Comprehensive capturing of the version and environments used for external programs references.")

    vande_dir = outdir + 'software_version_information/'
    if not os.path.isdir(vande_dir):
        os.system('mkdir %s' % vande_dir)
        seqc_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/'
        os.system('cp -r %s %s' % (seqc_dir + 'environment_provenance/current/', vande_dir + 'environment_info/'))
        os.system('cp -r %s %s' % (seqc_dir + 'prog_to_environment.txt', vande_dir))
        os.system('cp -r %s %s' % (seqc_dir + 'conda_environments.txt', vande_dir))

def clearAnyZombieJobs(outdir, cluster, logObject, wait=30):
    logObject.info("Killing/qdeling any zombie jobs related to the seQc Repo!")

    lsof_cmd = ['lsof', '+D', outdir]
    proc = subprocess.Popen(lsof_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    for i, line in enumerate(out.decode('utf-8').split('\n')):
        if i == 0: continue
        line = line.strip()
        ls = line.split()
        if len(ls) < 2: continue
        if not ls[1].strip() == str(os.getpid()):
            sys.stderr.write("Warning: running kill -9 %s\n" % ls[1])
            os.system('kill -9 %s' % ls[1])

    grid_dir = outdir + 'grid_submission_files/'
    if os.path.isdir(grid_dir):
        possible_pbs_jobs = set([])
        for f in os.listdir(grid_dir):
            if f.endswith('.job-id'):
                pbs_log_file = grid_dir + f
                for line in open(pbs_log_file):
                    line = line.strip()
                    if len(line.split()) == 1 and uF.is_number(line):
                        possible_pbs_jobs.add(line)

        if cluster == 'UGE':
            user_jobs = set([])
            qstat_cmd = ['qstat']
            proc = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            for i, line in enumerate(out.decode('utf-8').split('\n')):
                line = line.strip()
                ls = line.split()
                if i > 1 and len(ls) > 4:
                    if not 'E' in ls[4]:
                        user_jobs.add(ls[0])


            qdeling_jobs = False
            for j in user_jobs.intersection(possible_pbs_jobs):
                sys.stderr.write('Warning: running qdel %s\n' % j)
                os.system('qdel %s' % j)
                qdeling_jobs = True
            if qdeling_jobs:
                time.sleep(wait)
        elif cluster == 'SLURM':
            user_jobs = set([])
            squeue_cmd = ['squeue']
            proc = subprocess.Popen(squeue_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            for i, line in enumerate(out.decode('utf-8').split('\n')):
                line = line.strip()
                ls = line.split()
                if i > 0 and len(ls) > 4:
                    if not 'E' in ls[4]:
                        user_jobs.add(ls[0])

            scanceling_jobs = False
            for j in user_jobs.intersection(possible_pbs_jobs):
                sys.stderr.write('Warning: running scancel %s\n' % j)
                os.system('scancel %s' % j)
                scanceling_jobs = True
            if scanceling_jobs:
                time.sleep(wait)

def seQuoiaSheppard(meta, illumina_data, illumina_data_format, nanopore, workflow, outdir, poolsize, cluster, project, config, sample_config):
    # assert that input file exists and there is a directory with sequencing data
    try:
        assert(os.path.isfile(meta))
    except:
        sys.stderr.write("ERROR: Unable to locate input meta-specifying data file. Now exiting ...\n")
        sys.exit(1)

    meta_file = os.path.abspath(meta)
    outdir = os.path.abspath(outdir) + '/'

    # set up batch analysis directory structure
    if not os.path.isdir(outdir):
        os.system('mkdir %s' % outdir)
    else:
        sys.stderr.write('-*-'*20 + '\n')
        sys.stderr.write('Warning: Results directory already exists! Will delete failed steps for samples, and pick up what needs to still be run!\n')
        sys.stderr.write('-*-'*20 + '\n')

    # create logging object
    log_file = outdir + 'seQc.log'
    logObject = uF.createLoggerObject(log_file)

    # clear any currently running zombie jobs associated with the seQuoia output repository
    # and capture version information for external software
    clearAnyZombieJobs(outdir, cluster, logObject, wait=30)
    capture_external_versioning(outdir, logObject)

    logObject.info("Starting batch seQuoia analysis for meta-input file at: %s" % meta_file)
    sys.stdout.write("Starting batch seQuoia analysis for meta-input file at:\n%s\n" % meta_file)

    # Parse configuration file
    config_params_dict = {}
    try:
        assert(os.path.isfile(workflow))
        os.system('cp %s %s' % (workflow, outdir))
        if config and sample_config:
            logObject.error("Please provide only a config or sample_config file. The simultaneous use of both is not supported. Exiting now ...")
            sys.exit(1)
        elif config:
            assert(os.path.isfile(config))
            logObject.info("The following workflow is being run with config parameters: %s, making copy for future reference." % workflow)
            logObject.info("Config parameters provided in %s. Also copied to seQc repo directory." % config)
            os.system('cp %s %s' % (config, outdir))
            config_params_dict = read_config_file(config, workflow, logObject)
        elif sample_config:
            assert (os.path.isfile(sample_config))
            logObject.info("The following workflow is being run with sample specfic config parameters: %s, making copy for future reference." % workflow)
            logObject.info("Config parameters provided in %s. Also copied to seQc repo directory." % sample_config)
            os.system('cp %s %s' % (sample_config, outdir))
            config_params_dict = read_sample_config_file(sample_config, workflow, logObject)
        else:
            logObject.info("The following workflow is being run with default parameters: %s, making copy for future reference." % workflow)
    except:
        logObject.error("Problem locating workflow or reading configurations file. Exiting now ...")
        sys.exit(1)


    # Parse Illumina input sequencing data
    illumina_present = False
    nanopore_present = False
    il_data = None; np_data = None
    if illumina_data or illumina_data_format:
        # assert that the datatype is a valid specification and that the
        try:
            valid_illumina_datatypes = set(['illumina-paired', 'illumina-single', 'gp-directory', 'bam'])
            assert(illumina_data_format in valid_illumina_datatypes)
            assert(os.path.isfile(illumina_data) or os.path.isdir(illumina_data))
            il_data = os.path.abspath(illumina_data)
            if os.path.isdir(il_data):
                il_data += '/'
            illumina_present = True
            logObject.info("Short read Illumina data specified in %s format.\nLocation of this data is provided/listed at: %s" % (illumina_data_format, il_data))
            sys.stdout.write("Short read Illumina data specified in %s format.\nLocation of this data is provided/listed at:\n%s\n" % (illumina_data_format, il_data))
        except:
            logObject.error("The datatype of the short read sequencing data is not a valid entry or illumina data path provided did not correspond to a file / directory. Please fix and retry. Exiting now ...\n")
            sys.exit(1)

    # Parse Oxford Nanopore input sequencing data.
    if nanopore:
        try:
            assert(os.path.isfile(nanopore))
            np_data = os.path.abspath(nanopore)
            nanopore_present = True
            logObject.info("Nanopore data provided and can be found listed at: %s" % np_data)
            sys.stdout.write("Nanopore data provided and can be found listed at:\n%s\n" % np_data)
        except:
            logObject.error("Nanopore data is not provided in expected format. Exiting now ...")
            sys.exit(1)

    # read input data file and store meta information for each strain.
    metadata_parsed = read_and_store_metadata(meta_file, logObject)

    # match samples specified in input file to illumina sequencing data files.
    ildata_parsed = None
    if illumina_present:
        ildata_parsed = process_illumina_data(il_data, illumina_data_format, metadata_parsed, nanopore_present, logObject)

    # match samples specified in input file to nanopore sequencing data files.
    npdata_parsed = None
    if nanopore_present:
        npdata_parsed = process_nanopore_data(np_data, metadata_parsed, logObject)

    # prepare for submission of workflows with multiprocessing pool
    sample_data, sample_runs = prepare_for_submission(cluster, project, outdir, ildata_parsed, npdata_parsed, metadata_parsed, workflow, config_params_dict)

    # catalog multiple runs for same sample with different parameters, if needed.
    sample_to_runs_file = output_directory + 'sample_to_runs.txt'
    sample_to_runs_handle = open(sample_to_runs_file, 'w')
    for s in sample_runs.items():
        if len(s[1]) > 1:
            for i in s[1]:
                sample_to_runs_handle.write(s[0] + '\t' + i[0] + '\t' + i[1] + '\n')
    sample_to_runs_handle.close()
    if os.path.getsize(sample_to_runs_file) < 10:
        os.system('rm -f %s' % sample_to_runs_file)

    # create pool and submit data to pool worker
    logObject.info("Starting pool submission!!!")
    sys.stdout.write("Starting pool submission!!!\n")

    if len(sample_data) < poolsize:
        poolsize = len(sample_data)

    p = multiprocessing.Pool(poolsize)

    p.map(workflow_process, sample_data)

    try:
        p.map(workflow_process, sample_data)
    except:
        clearAnyZombieJobs(outdir, cluster, logObject, wait=0)
        logObject.error("User prompted KeyboardInterrupt! Running jobs successfully killed, Exiting!")
        p.close()
        os.exit(1)
    else:
        p.close()

    # seQcSheppard exiting!
    logObject.info("Done!")
    sys.stdout.write("\nDone!\n")
    sys.exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    This program serves as a wrapper to setup a general sequencing 
    batch repository and acts as a sheppard for a pool of workflows running simultaneously.
    
    If started in a cluster node, please make sure to "use" (Dotkit lingo) that cluster for the individual processes as well. 
    
    There are currently four ways to input short read Illumina data:
        - illumina-paired: this is when you have forward and reverse reads
        - illumina-single: this is when you have single end reads
        - gp-directory : this is when you have a directory from the Broad Genomics Platform.
        - bam : this is when you just have the path to a bam (aligned/unlaigned).
                
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-m', '--meta', help="Input sample meta-data table. Tab-delimited. Headers necessary and first column must be called \"sample_id\" baring the names of samples up until _R1.fastq (and _R2.fastq if paired-end data).",
                        required=True)
    parser.add_argument('-i', '--illumina', help="Input directory containing sequencing data. Names should match sample_id column in input. Alternatively, user can provide a tab-delimited file listing the input files. No column headers, but first column must match sampfrom meta",
                        required=False)
    parser.add_argument('-f', '--illumina_data_format', help="Specify the format for the illumina data.", choices=['illumina-paired', 'illumina-single', 'bam', 'gp-directory'], required=False)
    parser.add_argument('-n', '--nanopore', help="A tab-delimited file with three columns: (1) sample_id, (2) Albacore results directory for run, and (3) sample barcode ID for run.",
                        required=False)
    parser.add_argument('-w', '--workflow', help="Provide path to established workflow. Alternatively, provide your own workflow file!",
                        required=True)
    parser.add_argument('-o', '--outdir', help="Path to the output directory.",
                        required=True)
    parser.add_argument('-p', '--poolsize', type=int, help="Pool size. Number of samples to process simultaneously. Default is 1.", default=1,
                        required=False)
    parser.add_argument('-a', '--config', help="Specify parameter configurations to change in workflow. Must match IDs used in workflow.", required=False, default=None)
    parser.add_argument('-s', '--sample_config', help="Specify parameter configurations to change in workflow, per sample. Config parameters must match IDs used in workflow and sample IDs must match those provided in meta input.", required=False, default=None)
    parser.add_argument('-c', '--cluster', help='specify cluster to use. Default is "UGE". Other options include genomics "SLURM" and "Local"', default="UGE",
                        required=False)
    args = parser.parse_args()

    seQuoiaSheppard(args.meta, args.illumina, args.illumina_data_format, args.nanopore, args.workflow, args.outdir, args.poolsize, args.cluster, None, args.config, args.sample_config)