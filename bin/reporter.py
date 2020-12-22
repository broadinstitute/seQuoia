# reporter
# Infectious Diseases and Microbiome Platform
# Bacterial Genomics
# 02-05-2017
# Rauf Salamzade

import os
import sys
import argparse
import json
import subprocess
import numpy as np
import xlsxwriter
import pandas as pd
from collections import defaultdict

from seQuoia.other import usefulFunctions as uF
multiqc_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/seQuoia/multiqc/'
ext_wrapper_dir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/seQuoia/external_wrappers/'
curr_dir = os.getcwd()

def identify_searching_spaces(input, logObject):
    """ Get the seQuoia Sheppard QC directory locations or parent locations provided. """
    repo_dir = [os.path.abspath(input) + '/']
    if os.path.isfile(input):
        repo_dir = []
        with open(input) as io:
            for line in io:
                line = line.rstrip('\n')
                if os.path.isdir(line): repo_dir.append(os.path.abspath(line) + '/')
                else: sys.stderr.write("Path %s is not a directory, skipping ...")
    logObject.info("Will be using the following repo directories and subdirectories thereof:\n%s" % '\n'.join(repo_dir))
    logObject.info('-'*80)
    return repo_dir

def run_MultiQC(repo_dir, title, outdir, logObject):
    """ Run and extract some info from MultiQC. """

    title = '_'.join(title.split())

    inmportant_dirs = ['FastQC', 'ProcessGPDirectory', 'AdapterTrim', 'QualityTrim']
    outf = open(outdir + 'file_list_for_multiqc.txt', 'w')
    paired_end_flag = False
    for repo in repo_dir:
        for subdir, dirs, files in os.walk(repo, followlinks=True):
            for file in files:
                full_file_path = os.path.join(subdir, file)
                in_important_dir = False
                for imp in inmportant_dirs:
                    if '/' + imp + '/' in full_file_path:
                        in_important_dir = True
                if '/Symlink_Input/' in full_file_path and (full_file_path.endswith('_R2.fastq') or full_file_path.endswith('_R2.fastq.gz')): paired_end_flag = True
                if in_important_dir and not file.endswith(".hdf5") and not file.endswith('.fastq') and not file.endswith('.fastq.gz') and not file.endswith('.bam') and not file.endswith('.sam') and not file.endswith('.vcf'):
                    outf.write(full_file_path + '\n')
    outf.close()

    # Run MultiQC
    stdout_file = outdir + 'multiqc.log'; stderr_file = outdir + 'multiqc.err'
    stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
    os.chdir(outdir)
    multiqc_cmd = ['bash', ext_wrapper_dir + 'multiqc.sh', '--force', '--title', title, '--config', multiqc_dir + 'seQuoia_multiqc.yaml']
    if not paired_end_flag: multiqc_cmd += ['--ignore-samples', '*_R2']
    modules = ['fastqc', 'picard', 'trimmomatic', 'cutadapt'] # metadata should only be for SPARC
    for mod in modules:
        multiqc_cmd += ['-m', mod]
    multiqc_cmd += ['--file-list', outdir + 'file_list_for_multiqc.txt']
    logObject.info('Running the following cmd for MultiQC: %s' %  ' '.join(multiqc_cmd))
    proc = subprocess.Popen(multiqc_cmd, stdout=stdout_handle, stderr=stderr_handle)
    proc.wait()
    stdout_handle.close(); stderr_handle.close()
    os.chdir(curr_dir)

    try: assert(len([x for x in os.listdir(outdir) if x.endswith('.html')]) > 0); logObject.info("MultiQC successfully ran!")
    except: logObject.error("Something went wrong with MultiQC. Please check the log/err files associated."); sys.exit(1)

    # Parse MultiQC general stats
    mqc_results_dir = outdir + title + '_multiqc_report_data/'
    mqc_general_stats_file = None
    if os.path.isfile(mqc_results_dir + "multiqc_general_stats.txt"): mqc_general_stats_file = mqc_results_dir + "multiqc_general_stats.txt"
    elif os.path.isfile(mqc_results_dir + "multiqc_fastqc.txt"): mqc_general_stats_file =  mqc_results_dir + "multiqc_fastqc.txt"
    mqc_data = defaultdict(dict)
    col_to_name = {}
    mqc_columns = []
    if not mqc_general_stats_file: raise RuntimeError()
    for i, line in enumerate(open(mqc_general_stats_file)):
        line = line.rstrip('\n')
        ls = line.split('\t')
        if i == 0:
            for j, val in enumerate(ls[1:]):
                col_to_name[j] = val
                mqc_columns.append(val)
        else:
            samp = ls[0]
            for j, val in enumerate(ls[1:]):
                mqc_data[samp][col_to_name[j]] = val

    logObject.info("MultiQC data successfully extracted!\n")
    logObject.info('-' * 70)
    return [mqc_data, mqc_columns]

def extract_metadata_and_centrifuge(repo_dir, outdir, logObject):
    """ Recursively look through input directories """
    meta_columns = []
    metadata = {}
    centrifugedata= {}
    cent_file = outdir + "centrifuge_summary.txt"
    cent_handle = open(cent_file, 'w')
    for repo in repo_dir:
        for subdir, dirs, files in os.walk(repo, followlinks=True):
            for file in files:
                if file == 'meta_information.json':
                    sample = subdir.split('/')[-1]
                    with open(os.path.join(subdir, file)) as of:
                        data = json.load(of)
                        metadata[sample] = data
                        meta_columns = data.keys()
                elif file.endswith("_kraken_report.txt"):
                    sample = subdir.split('/')[-2]
                    centrifuge_results = uF.kraken_report_parser(os.path.join(subdir, file))
                    centrifugedata[sample] = centrifuge_results[:-1]
                    for tc in centrifuge_results[-1]:
                        for t in centrifuge_results[-1][tc]:
                            cent_handle.write('\t'.join([sample, tc, t, str(float(centrifuge_results[-1][tc][t]) / centrifuge_results[1])]) + '\n')
    cent_handle.close()
    logObject.info("Meta-data and Centrifuge results (if run) have been parsed!")
    return [meta_columns, metadata, centrifugedata]

def create_high_level_stats_ont(metadata, meta_columns, sample_assembly_metrics, sample_GAEMR_htmls, data_headers, outdir, logObject):
    final_file = outdir + "high_level_stats_table.txt"
    final_handle = open(final_file, 'w')
    headers = ['sample_id', 'assembly'] + data_headers + ['GAEMR_index_locations']
    final_handle.write('\t'.join(headers) + '\n')
    for sample in sample_assembly_metrics:
        for assembly in sample_assembly_metrics[sample]:
            printlist = [sample, assembly]
            for h in data_headers:
                printlist.append(sample_assembly_metrics[sample][assembly][h])
            printlist.append(sample_GAEMR_htmls[sample][assembly])
            final_handle.write('\t'.join(printlist) + '\n')
    final_handle.close()
    logObject.info("High level statistics file created and ready for visualization with ShinyGAEMR application!")

def create_high_level_stats(metadata, mqc_data, mqc_columns, meta_columns, centrifugedata, outdir, logObject):
    """ Create high level statistic file with all the meta data and just the important sequencing stats. """
    final_file = outdir + "high_level_stats_table.txt"
    final_handle = open(final_file, 'w')
    headers =['sample_id', 'sample'] + sorted(mqc_columns) + sorted(meta_columns)
    if len(centrifugedata) > 0: headers += ['final_total_reads', 'unclassified', 'classified', 'H._sapiens_mapping', 'E._coli_mapping']
    final_handle.write('\t'.join(headers) + '\n')
    data = [headers]
    for readset in mqc_data:
        sample = readset.split('_R1')[0].split('_R2')[0]
        printlist = [readset, sample]
        for c in sorted(mqc_columns):
            printlist.append(mqc_data[readset][c])
        for c in sorted(meta_columns):
            value = ''
            if c in metadata[sample]:
                value = metadata[sample][c]
            printlist.append(value)
        if len(centrifugedata) > 0:
            unclassified_reads, classified_reads, human_reads, ecoli_reads = ["NA"]*4
            if sample in centrifugedata.keys():
                unclassified_reads, classified_reads, human_reads, ecoli_reads = centrifugedata[sample]
            total_reads = unclassified_reads + classified_reads
            printlist += [str(x) for x in [total_reads, unclassified_reads, classified_reads, human_reads, ecoli_reads]]
        data.append(printlist)
        printlist_for_r = [x.replace("'", "").replace("#", '') for x in printlist]
        final_handle.write('\t'.join(printlist_for_r) + '\n')
    final_handle.close()
    logObject.info("High level statistics file created and ready for visualization with seeQc Shiny app!")
    return [data, headers]

def determine_outliers(data, logObject):
    """ Determine outliers for numeric statistics using MAD approach """
    col_to_samp = {}
    outlier_samps = {}
    for i, vals in enumerate(zip(*data)):
        #print (vals)
        if i == 0 or i == 1:
            if i == 0:
                for j, val in enumerate(vals[1:]):
                    col_to_samp[j] = val
        elif i > 1:
            num_numeric = 0
            num_nonnumeric = 0
            for val in vals[1:]:
                if val != 'NA' and uF.is_number(val):
                    num_numeric += 1
                elif val != 'NA': num_nonnumeric += 1
            if num_numeric > 0 and num_nonnumeric == 0:
                new_vals_existant = []
                new_vals_existant_samps = []
                for j, x in enumerate(vals[1:]):
                    if x != 'NA': new_vals_existant.append(float(x)); new_vals_existant_samps.append(j)
                new_vals_existant = np.array(new_vals_existant)
                try:
                    outliers = uF.mad_based_outlier(new_vals_existant)
                    outlier_samps[vals[0]] = set([])
                    for j, outlier_call in enumerate(outliers):
                        if outlier_call: outlier_samps[vals[0]].add(col_to_samp[new_vals_existant_samps[j]])
                except RuntimeWarning:
                    # looks like no outliers were found
                    pass
    logObject.info("Outlier samples have been detected using the MAD approach!")
    return outlier_samps

def generate_excel_spreadsheet(outlier_samps, headers, data, outdir, logObject):
    """ Create an excel spreadsheet with outliers marked in red and headers nicely formatted!"""
    data_frame = pd.DataFrame(data, columns=headers)
    writer = pd.ExcelWriter(outdir + 'all_general_stats.xlsx', engine='xlsxwriter')
    data_frame.to_excel(writer, sheet_name="All_General_Stats", startrow=0, header=False, index=False,startcol=1)

    workbook = writer.book
    worksheet = writer.sheets['All_General_Stats']

    header_format = workbook.add_format({
        'bold': True,
        'font_size': 14,
        'text_wrap': True,
        'valign': 'top',
        'align': 'center',
        'fg_color': '#cce6ff',
        'border': 1})

    outlier_format = workbook.add_format({
        'bold': True,
        'font_size': 12,
        'text_wrap': True,
        'fg_color': '#ff6666',
        'border': 1})

    for col in range(100):
        worksheet.set_column(col, col, 20)

    for col_num, value in enumerate(data_frame.columns.values):
        worksheet.write(0, col_num+1, value, header_format)

    for (row_num, row) in data_frame.iterrows():
        if row_num > 0:
            sample_file = row.loc['sample_id']
            for col_num, col in enumerate(data_frame.columns.values):
                if col in outlier_samps and sample_file in outlier_samps[col]:
                   worksheet.write(row_num, col_num+1, row.loc[col], outlier_format)
    logObject.info("Excel spreadsheet created with high level statistics and numeric categories having any outliers marked!")

def parse_nano_assembly_structure(repo_dir, logObject):
    """ Parse stats from seQuoia sheppard repositories for Bacterial Nanopore Assembly after reorganization. """

    unorganized_naming = {'full-np': 'Unicycler_All-ONT', 'sub-np': 'Unicycler_Subsampled-ONT', 'canu': 'Canu_Pure-ONT'}

    sample_metrics = {}
    gaemr_paths = {}
    for repo in repo_dir:
        logObject.info("Beginning parsing of seQuoia repository %s" % repo)
        for s in os.listdir(repo):
            sample_path = repo + s + '/'
            sample_checkpoint = sample_path + 'SETUP.txt'
            reorg_checkpoint = sample_path + 'REORGANIZATION.txt'
            gaemr_dirs = {}
            if os.path.isfile(sample_checkpoint) and os.path.isfile(reorg_checkpoint):
                assembly_dir = sample_path + 'Assembly_Results/'
                for a in os.listdir(assembly_dir):
                    gaemr_QC_dir = assembly_dir + a + '/QC/'
                    gaemr_dirs[a] = gaemr_QC_dir
            elif os.path.isfile(sample_checkpoint):
                for a in [x for x in os.listdir(sample_path) if x.startswith('GAEMR')]:
                    suffix = a.split('GAEMR_')[1]
                    gaemr_checkpoint = sample_path + 'GAEMR_' + suffix + '.txt'
                    if os.path.isfile(gaemr_checkpoint):
                        aname = unorganized_naming[suffix]
                        gaemr_QC_dir = sample_path + '/' + a + '/QC/'
                        gaemr_dirs[aname] = gaemr_QC_dir
            else:
                continue
            sample_metrics[s] = {}
            gaemr_paths[s] = {}
            for a, apath in gaemr_dirs.items():
                if not os.path.isdir(apath): continue
                gaemr_qc_path = apath + '/gaemr/'
                gaemr_html_path = apath + '/gaemr/html/index_ont.html'
                sample_metrics[s][a] = defaultdict(lambda: "NA")

                seq_cov_f = gaemr_qc_path + 'table/assembly.scaffolds.bam_seq_cvg_stats.table.txt'
                assembly_stats_f = gaemr_qc_path + 'table/assembly.basic_assembly_stats.table.txt'
                rrna_f = gaemr_qc_path + 'table/assembly.rna_analysis_summary.table.txt'
                frag_map_f = gaemr_qc_path + 'table/assembly.scaffolds.simple_bam_stats.table.txt'

                try:
                    assert(os.path.isfile(seq_cov_f) and os.path.isfile(assembly_stats_f) and os.path.isfile(rrna_f) and os.path.isfile(frag_map_f) and os.path.isfile(gaemr_html_path))
                except:
                    logObject.warning("Skipping sample %s, because one of the expected GAEMR files is missing" % s)
                    continue
                else:
                    logObject.info("Harnessing GAEMR QC stats for sample %s" % s)

                gaemr_paths[s][a] = gaemr_html_path

                # get seq_cvg
                index_1 = None
                index_2 = None
                with open(seq_cov_f) as seq_cov_of:
                    for i, line in enumerate(seq_cov_of):
                        line = line.strip()
                        ls = line.split(' | ')
                        if line.startswith('#Stat'):
                            if 'Fragments' in ls[1]: index_1 = 'frag'
                            elif 'Long' in ls[1]: index_1 = 'long'
                            if len(ls) > 2 and 'Fragments' in ls[2]: index_2 = 'frag'
                            elif len(ls) > 2 and  'Long' in ls[2]: index_2 = 'long'
                        elif line.startswith('Median Cvg'):
                            if index_1: median_cvg_1 = line.split(' | ')[1]; sample_metrics[s][a]['seq_cvg_' + index_1] = median_cvg_1
                            if index_2: median_cvg_2 = line.split(' | ')[2]; sample_metrics[s][a]['seq_cvg_' + index_2] = median_cvg_2

                # get contigs, contig_max, contig_ln, contig_n50
                with open(assembly_stats_f) as assembly_stats_of:
                    for i, line in enumerate(assembly_stats_of):
                        line = line.strip()
                        ls = line.split(' | ')
                        if line.startswith("Contigs"):
                            sample_metrics[s][a]['contigs'] = ls[1]
                        elif line.startswith("Max Contig"):
                            sample_metrics[s][a]['contig_max'] = ls[1]
                        elif line.startswith("Total Contig Length"):
                            sample_metrics[s][a]['contig_ln'] = ls[1]
                        elif line.startswith("Contig N50"):
                            sample_metrics[s][a]['contig_n50'] = ls[1]

                # get 16s & 16s_copies
                sample_metrics[s][a]['16s'] = []
                sample_metrics[s][a]['16s_copies'] = []
                with open(rrna_f) as rrna_of:
                    for i, line in enumerate(rrna_of):
                        line = line.strip()
                        ls = line.split(' | ')
                        if not line.startswith('#'):
                            if ls[2] == 'genus':
                                taxonomy = ls[-1]
                                copy_number = ls[1]
                                sample_metrics[s][a]['16s'].append(taxonomy)
                                sample_metrics[s][a]['16s_copies'].append(copy_number)
                sample_metrics[s][a]['16s'] = '; '.join(sample_metrics[s][a]['16s'])
                sample_metrics[s][a]['16s_copies'] = '; '.join(sample_metrics[s][a]['16s_copies'])

                # get pct_map
                with open(frag_map_f) as frag_map_of:
                    for i, line in enumerate(frag_map_of):
                        line = line.strip()
                        ls = line.split(' | ')
                        if ls[0] == 'Mapped':
                            if index_1: sample_metrics[s][a]['pct_' + index_1 + '_map'] = ls[1].split()[1].replace('(', '').replace(')', '')
                            if index_2: sample_metrics[s][a]['pct_' + index_2 + '_map'] = ls[2].split()[1].replace('(', '').replace(')', '')
    data_headers = ['16s', '16s_copies', 'contigs', 'contig_max', 'contig_ln', 'contig_n50', 'pct_frag_map', 'pct_long_map', 'seq_cvg_frag', 'seq_cvg_long']
    return [sample_metrics, gaemr_paths, data_headers]

def Reporter(input, outdir, analysis_type, title):
    try: assert(os.path.isfile(input) or os.path.isdir(input))
    except: sys.stderr.write("Input was neither a file nor directory. Exiting now ...")

    # create results directory
    outdir = os.path.abspath(outdir) + '/'
    if not os.path.isdir(outdir): os.system('mkdir ' + outdir)
    else: sys.stderr.write('Warning: Results directory already exists! Press control-C multiple times repeatedly and panic!\n... or just let the data be overwritten.\n')

    # create logging object
    log_file = outdir + 'Reporter.log'
    logObject = uF.createLoggerObject(log_file)

    # get the sheppard QC directory locations or parent locations provided
    repo_dir = identify_searching_spaces(input, logObject)

    # read meta data and centrifuge results
    meta_columns, metadata, centrifugedata = extract_metadata_and_centrifuge(repo_dir, outdir, logObject)

    if analysis_type == 'basic':
        # run MultiQC and extract some data from the summary text files
        mqc_data, mqc_columns = run_MultiQC(repo_dir, title, outdir, logObject)

        # create high level statistics file for visualization in R Markdown or seeQc Shiny Application
        general_stats_data, general_stats_categories = create_high_level_stats(metadata, mqc_data, mqc_columns, meta_columns, centrifugedata, outdir, logObject)

        # determine outliers using MAD approach
        outlier_samps = determine_outliers(general_stats_data, logObject)

        # create Pandas data frame and write to excel spreadsheet.
        generate_excel_spreadsheet(outlier_samps, general_stats_categories, general_stats_data, outdir, logObject)

    elif analysis_type == 'nano-assembly':
        # parse seQuoia repos for GAEMR metrics
        sample_assembly_metrics, sample_GAEMR_htmls, data_headers = parse_nano_assembly_structure(repo_dir, logObject)

        # create high level statistics file for visualization in Shiny Application
        create_high_level_stats_ont(metadata, meta_columns, sample_assembly_metrics, sample_GAEMR_htmls, data_headers, outdir, logObject)


if __name__ == '__main__':
    # Pull out the arguments.

    parser = argparse.ArgumentParser(description="""
    This is the primary seQuoia program for downstream processing analysis of seQuoia repo data. It primarily performs three
    tasks:
    
    1.) Generates a MultiQC report for easy sharing of basic sequencing QC metrics
    
    2.) Outlier detection and reporting in comprehensive Excel spreadsheet 
    
    3.) Provides a tab-delimited file [ reults/high_level_stats_table.txt ] for downstream visualizing with R Markdown or R Shiny
    
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', help="Either the output directory of a batch sheppard.py run or a list of individual seQuoia repositories as generated by seQuoia sheppard - one per line.", required=True)
    parser.add_argument('-o', '--outdir', help="Output directory where to generate data.", required=True)
    parser.add_argument('-t', '--analysis_type', help="Specify the analysis type.", choices=['basic', 'nano-assembly'], required=True)
    parser.add_argument('-n', '--title', help="MultiQC report title! Default is \"MultiQC Report\".", default="MultiQC Report", required=False)

    args = parser.parse_args()

    Reporter(args.input, args.outdir, args.analysis_type, args.title)