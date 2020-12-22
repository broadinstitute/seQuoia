import os
import sys
import argparse
import gzip
from seQuoia.other import usefulFunctions as uF
from seQuoia.classes import AssemblyAnalyzer

def runAssemblyQC(assembly, sample_name, sample_dir, format_options, qc_options, cores, illumina_frw, illumina_rev, picard_insert_file, nanopore_fastq, gaemr_ont, identifier):
    try: assert(os.path.isfile(assembly))
    except:
        sys.stderr.write("ERROR: Assembly does not exist. Raising exception\n")
        raise RuntimeError

    try: assert((not nanopore_fastq) or (nanopore_fastq and os.path.isfile(nanopore_fastq)))
    except:
        sys.stderr.write("ERROR: Nanopore FASTQ provided but the path does not exist. Raising exception\n")
        raise RuntimeError

    format_options = format_options.strip('"')
    qc_options = qc_options.strip('"')

    # set up directory structure
    workspace_name = "GAEMR"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(sample_dir, workspace_name)

    # create logging object
    log_file = workspace + 'GAEMR.log'
    logObject = uF.createLoggerObject(log_file)

    # Format Assembly for GAEMR QC Analysis
    workspace_a = uF.setupDirectory(workspace, "Formatting")
    AssemblyObj = AssemblyAnalyzer.Assembly(assembly, sample_name, logObject)
    AssemblyObj.run_gaemr_formatter(workspace_a, options=format_options, reference_change=True)

    # Generate read list file for QC
    read_list = workspace + 'read_list.txt'
    outf = open(read_list, 'w')
    outf.write("#name,lib_type,mean_read_length,dir,insert_size,files\n")
    frw_read = None
    rev_read = None
    illumina_avg_read_length = 250
    illumina_avg_insert_length = 400
    if illumina_frw and illumina_rev and os.path.isfile(illumina_frw) and os.path.isfile(illumina_rev):
        readlengths = []
        if illumina_frw.endswith(".gz"):
            with gzip.open(illumina_frw, 'rt') as ofr:
                for i, line in enumerate(ofr):
                    if i > 40000: continue
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))
        else:
            with open(illumina_frw) as ofr:
                for i, line in enumerate(ofr):
                    if i > 40000: continue
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))
        frw_read = illumina_frw
        rev_read = illumina_rev
        illumina_avg_read_length = sum(readlengths)/float(len(readlengths))

    elif os.path.isdir(sample_dir + 'Subsample'):
        frw_read = [sample_dir + 'Subsample/' + x for x in os.listdir(sample_dir + 'Subsample') if '_R1.' in x and (x.endswith('.fastq.gz') or x.endswith('.fastq'))][0]
        rev_read = [sample_dir + 'Subsample/' + x for x in os.listdir(sample_dir + 'Subsample') if '_R2.' in x and (x.endswith('.fastq.gz') or x.endswith('.fastq'))][0]

        readlengths = []
        if frw_read.endswith(".gz"):
            with gzip.open(frw_read, 'rt') as ofr:
                for i, line in enumerate(ofr):
                    if i > 40000: continue
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))
        else:
            with open(frw_read) as ofr:
                for i, line in enumerate(ofr):
                    if i > 40000: continue
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))

        illumina_avg_read_length = sum(readlengths)/float(len(readlengths))

    if picard_insert_file and os.path.isfile(picard_insert_file):
        with open(picard_insert_file) as oisf:
            for line in oisf:
                line = line.strip()
                ls = line.split('\t')
                if len(ls) > 0 and ls[0].startswith("MEDIAN_INSERT_SIZE"): flag_header_observed = True; continue
                if flag_header_observed:
                    illumina_avg_insert_length = int(float(ls[5]))
                    break

    elif os.path.isdir(sample_dir + 'ProcessGPDirectory'):
        insert_stats_file_query = [sample_dir + 'ProcessGPDirectory/' + x for x in os.listdir(sample_dir + 'ProcessGPDirectory') if x.endswith('.insert_size_metrics')]
        if len(insert_stats_file_query) == 1:
            insert_stats_file = insert_stats_file_query[0]
            flag_header_observed = False
            with open(insert_stats_file) as oisf:
                for line in oisf:
                    line = line.strip()
                    ls = line.split('\t')
                    if len(ls) > 0 and ls[0].startswith("MEDIAN_INSERT_SIZE"): flag_header_observed = True; continue
                    if flag_header_observed:
                        illumina_avg_insert_length = int(float(ls[5]))
                        break

    if frw_read and rev_read:
        outf.write('Fragments,fragment,%d,fr,%d,%s,%s\n' % (illumina_avg_read_length, illumina_avg_insert_length, frw_read, rev_read))

    if nanopore_fastq and gaemr_ont:
        readlengths = []
        if nanopore_fastq.endswith(".gz"):
            with gzip.open(nanopore_fastq, 'rt') as ofr:
                for i, line in enumerate(ofr):
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))
        else:
            with open(nanopore_fastq) as ofr:
                for i, line in enumerate(ofr):
                    if i > 0 and (i+1) % 2 == 0 and (i+1) % 4 != 0:
                        readlengths.append(len(line.strip()))

        nanopore_avg_read_length = 2000
        if len(readlengths) > 0:
            nanopore_avg_read_length = sum(readlengths)/len(readlengths)

        nanopore_avg_insert_size = nanopore_avg_read_length

        outf.write('Long,unpaired,%d,,%d,%s\n' % (nanopore_avg_read_length, nanopore_avg_insert_size, nanopore_fastq))

    outf.close()

    # Run GAEMR QC
    workspace_b = uF.setupDirectory(workspace, "QC")
    if nanopore_fastq and gaemr_ont:
        AssemblyObj.run_gaemr_qc_ont(read_list, workspace_b, options=qc_options, cores=cores)
    else:
        AssemblyObj.run_gaemr_qc(read_list, workspace_b, options=qc_options, cores=cores)

    # create successful completion file if steps completed!
    conf_file_name = sample_dir + "GAEMR"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("GAEMR: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for Assembly QC with GAEMR. Uses nanoporeFunctions.py.
    """)

    parser.add_argument('-a', '--assembly', help='location of (?hybrid) ONT assembly.', required=True)
    parser.add_argument('-s', '--sample', help="Sample name.", required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-gf', '--format_options', help='formatting options for GAEMR.', required=False, default="-g 1 -c 100 -r")
    parser.add_argument('-gq', '--qc_options', help='QC options for GAEMR.', required=False, default="--force --analyze_rna")
    parser.add_argument('-c', '--cores', type=int, help='number of cores to allocate.', required=False, default=1)
    parser.add_argument('-frw', '--fastq_frw', help='illumina frw reads', required=False, default=None)
    parser.add_argument('-rev', '--fastq_rev', help='illumina rev reads', required=False, default=None)
    parser.add_argument('-p', '--picard_insert_file', help='provide the picard_insert_file', required=False, default=None)
    parser.add_argument('-n', '--nanopore_fastq', help='nanopore fastq', required=False, default=None)
    parser.add_argument('-f', '--gaemr_ont', action="store_true", help='whether to run GAEMR_ont.py instead of GAEMR.py. Currently only available on Broad servers.', required=False, default=False)
    parser.add_argument('-i', '--identifier', help="identifier, in case multiple runs in workflow.", required=False, default=None)
    args = parser.parse_args()

    runAssemblyQC(args.assembly, args.sample, args.sample_dir, args.format_options, args.qc_options, args.cores, args.fastq_frw, args.fastq_rev, args.picard_insert_file, args.nanopore_fastq, args.gaemr_ont, args.identifier)
