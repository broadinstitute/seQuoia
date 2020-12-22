import os
import sys
import argparse
import zipfile

from collections import defaultdict
from seQuoia.other import usefulFunctions as uF
from fadapa import Fadapa

def reorganize(sample_dir):
    try: assert(os.path.isdir(sample_dir))
    except: sys.stderr.write("ERROR: Sample directory doesn't seem to exist! Exiting now ...\n"); raise RuntimeError

    sample_dir = os.path.abspath(sample_dir) + '/'

    # set up directory structure
    workspace_name = "LSARP_Results/"
    workspace = sample_dir + workspace_name
    if not os.path.isdir(workspace):
        workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'LSARP_Table_Creation.log'
    logObject = uF.createLoggerObject(log_file)

    sample = sample_dir.split('/')[-2]
    logObject.info("Creating easy upload formats for sample %s", sample)
    logObject.info("-"*80)

    # FASTQC Tables

    logObject.info('Creating FastQC Data Tables.')
    logObject.info('-'*80)

    FastQC_results = 'FastQC/'
    FastQC_results_workspace = workspace + FastQC_results
    fastqc_modules = ['Per base sequence quality', 'Per tile sequence quality', 'Per sequence quality scores',
                      'Per base sequence content', 'Per sequence GC content', 'Per base N content',
                      'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences',
                      'Adapter Content']

    try:
        fastqc_zipped_data_dirs = [sample_dir + 'FastQC/' + zd for zd in os.listdir(sample_dir + 'FastQC/') if zd.endswith('.zip')]
        assert(len(fastqc_zipped_data_dirs) > 0)
        for zd in fastqc_zipped_data_dirs: assert(os.path.isfile(zd))
        if not os.path.isdir(FastQC_results_workspace):
            FastQC_results_workspace = uF.setupDirectory(workspace, FastQC_results)
    except: logObject.error('No FastQC results available or path is unable to be determined!')
    else:
        for zd in fastqc_zipped_data_dirs:
            with zipfile.ZipFile(zd) as z:
                for filename in z.namelist():
                    if filename.split('/')[-1] == 'fastqc_data.txt':
                        with z.open(filename) as fh:
                            FastQC_tmp_out = open(FastQC_results_workspace + 'tmp.txt', 'wb')
                            for line in fh:
                                FastQC_tmp_out.write(line)
                            FastQC_tmp_out.close()
                            fadapa = Fadapa(FastQC_results_workspace + 'tmp.txt')
                            for module in fastqc_modules:
                                try:
                                    table_file = '_'.join(module.split())
                                    cleaned_module_data = fadapa.clean_data(module)
                                    if cleaned_module_data:
                                        table_handle = open(FastQC_results_workspace + table_file + '.table.txt', 'w')
                                        for i, split_line in enumerate(cleaned_module_data):
                                            if i == 0:
                                                split_line = ['sample', 'read'] + split_line
                                            else:
                                                split_line = [sample_dir.split('/')[-2], zd.split('/')[-1].split(sample_dir.split('/')[-2] + '_')[1].split('_fastqc.zip')[0].split('.')[0]] + split_line
                                            table_handle.write('\t'.join(split_line) + '\n')
                                        table_handle.close()
                                except:
                                    pass
                            os.system('rm -f %s' % FastQC_results_workspace + 'tmp.txt')
    logObject.info('*'*80)

    # Centrifuge Tables

    logObject.info('Creating Centrifuge Data Tables.')
    logObject.info('-' * 80)

    Centrifuge_results = 'Centrifuge/'
    Centrifuge_results_workspace = workspace + Centrifuge_results

    centrifuge_report_file = sample_dir + 'Centrifuge/' + sample_dir.split('/')[-2] + '_centrifuge_report.tsv'
    kraken_report_file = sample_dir + 'Centrifuge/' + sample_dir.split('/')[-2] + '_centrifuge_kraken_report.txt'

    try:
        assert (os.path.isfile(centrifuge_report_file) and os.path.isfile(kraken_report_file))
        if not os.path.isdir(Centrifuge_results_workspace):
            Centrifuge_results_workspace = uF.setupDirectory(workspace, Centrifuge_results)

        centrifuge_report_table_file = Centrifuge_results_workspace + 'centrifuge_report.table.txt'
        centrifuge_report_table_handle = open(centrifuge_report_table_file, 'w')

        centrifuge_report_data = defaultdict(lambda: ['NA']*6)
        for i, line in enumerate(open(centrifuge_report_file)):
            if i > 0:
                line = line.rstrip('\n')
                name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance = line.split('\t')
                centrifuge_report_data[name] = [taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance]

        header = ['sample', 'taxonomy_name', 'taxonomy_level', 'taxonomy_rank', 'taxonomy_id', 'genome_size',
                  'centrifuge_abundance', 'percentage_of_fragments_recursively_covered',
                  'number_of_fragments_recursively_included', 'number_of_fragments_direct']
        centrifuge_report_table_handle.write('\t'.join(header) + '\n')
        for i, line in enumerate(open(kraken_report_file)):
            line = line.rstrip('\n')
            prop, frag_recurse, frag_direct, tax_level, tax_id = line.split()[:5]
            tax = ' '.join(line.split()[5:]).strip()
            taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance = centrifuge_report_data[tax]
            centrifuge_report_table_handle.write('\t'.join([sample_dir.split('/')[-2], tax, tax_level, taxRank, taxID, genomeSize, abundance, prop, frag_recurse, frag_direct]) + '\n')

        centrifuge_report_table_handle.close()
    except:
        logObject.error('No Centrifuge results available!')

    logObject.info('*' * 80)

    # AMRP Tables

    logObject.info('Moving Results from ARIBA and ShortBRED AMR Searches.')
    logObject.info('-' * 80)

    AMRP_results = 'AMRP_Searches/'
    AMRP_results_workspace = workspace + AMRP_results

    try:
        AMRP_dir = sample_dir + 'AMRP/'
        assert(os.path.isdir(AMRP_dir))
        if not os.path.isdir(workspace + AMRP_dir):
            AMRP_results_workspace = uF.setupDirectory(workspace, AMRP_results)
        for sd in os.listdir(AMRP_dir):
            ariba_dir = AMRP_dir + sd + '/'
            ariba_report = ariba_dir + 'report.tsv'
            if os.path.isfile(ariba_report):
                ariba_result = AMRP_results_workspace + sample_dir.split('/')[-2] + '_' + sd + '_ariba_results.txt'
                os.system('cp %s %s' % (ariba_report, ariba_result))
    except:
        logObject.error('Unable to create AMR prediction data tables.')# Raising exception now ...')

    logObject.info('*' * 80)

    # MLST Tables

    logObject.info('Creating MLST Data Tables.')
    logObject.info('-' * 80)

    MLST_results = 'MLST/'
    MLST_results_workspace = workspace + MLST_results

    try:
        MLST_dir = sample_dir + 'MLST/'
        MLST_result_file = MLST_dir + 'ariba_mlst/mlst_report.tsv'

        if not os.path.isdir(MLST_results_workspace):
            MLST_results_workspace = uF.setupDirectory(workspace, MLST_results)
        os.system('cp %s %s' % (MLST_result_file, MLST_results_workspace))

    except:
        logObject.error('Unable to create MLST call data tables.')# Raising exception now ...')
        #raise RuntimeError

    logObject.info('*' * 80)

    # De Novo Assembly Storage

    logObject.info('Moving de novo assembly to results directory.')
    logObject.info('-' * 80)

    Assembly_results = 'Assembly/'
    Assembly_results_workspace = workspace + Assembly_results

    try:
        Assembly_dir = sample_dir + 'Assembly/'
        Assembly_original_location = Assembly_dir + 'assembly.fasta'
        if not os.path.isfile(Assembly_original_location):
            Assembly_original_location = Assembly_dir + 'scaffolds.fasta'
        assert(os.path.isfile(Assembly_original_location))
        if not os.path.isdir(Assembly_results_workspace):
            Assembly_results_workspace = uF.setupDirectory(workspace, Assembly_results)
        Assembly_new_location = Assembly_results_workspace + sample_dir.split('/')[-2] + '.genome.fa'
        os.system('cp %s %s' % (Assembly_original_location, Assembly_new_location))
    except:
        logObject.error('Unable to move assembly to results directory.')

    logObject.info('*' * 80)

    # Assembly QC Storage

    logObject.info('Moving GAEMR assembly QC to results directory.')
    logObject.info('-' * 80)

    try:
        Assembly_QC_new_location = workspace + 'Assembly_QC/'
        Assembly_QC_original_dir = sample_dir + 'GAEMR/QC/'
        assert(os.path.isdir(Assembly_QC_original_dir))
        os.system('cp -r %s %s' % (Assembly_QC_original_dir, Assembly_QC_new_location))

    except:
        logObject.error('Unable to move GAEMR assembly QC to results directory.')

    logObject.info('*' * 80)

    # Pilon Results Storage

    logObject.info('Moving Pilon output to results directory.')
    logObject.info('-' * 80)

    try:
        Pilon_new_dir = workspace + 'Reference_Assembly_and_Variant_Calling/'
        Pilon_original_dir = sample_dir + 'Pilon/results/'
        assert(os.path.isdir(Pilon_original_dir))
        os.system('cp -r %s %s' % (Pilon_original_dir, Pilon_new_dir))
        os.system('gzip %s*' % Pilon_new_dir)

    except:
        logObject.error('Unable to move Pilon output to results directory.')

    logObject.info('*' * 80)

    # StrainGST Results Storage

    logObject.info('Moving StrainGST output to results directory.')
    logObject.info('-' * 80)

    try:
        Straingst_result_file = sample_dir + 'StrainGST/' + sample + '.straingst_result.tsv'
        assert(os.path.isfile(Straingst_result_file))
        Straingst_new_dir = 'StrainGST/'
        Straingst_results_workspace = workspace + Straingst_new_dir 
        if not os.path.isdir(workspace + Straingst_new_dir):
            Straingst_results_workspace = uF.setupDirectory(workspace, Straingst_new_dir)
        os.system('cp %s %s' % (Straingst_result_file, Straingst_results_workspace))

    except:
        logObject.error('Unable to move StrainGST output to results directory.')

    logObject.info('*' * 80)

    uF.closeLoggerObject(logObject)

    # create successful completion file if steps completed!
    conf_file = open(sample_dir + "LSARP.txt", 'w')
    conf_file.write("LSARP Table Creation: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    This module creates easy-to-parse tables for LSARP pipeline.
    """)

    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    args = parser.parse_args()

    reorganize(args.sample_dir)
