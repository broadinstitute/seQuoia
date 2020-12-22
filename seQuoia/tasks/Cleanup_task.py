import os
import sys
import argparse
import json
from collections import defaultdict
from seQuoia.other import usefulFunctions as uF

DEFAULT_SUFFIX_PURGE_FILE = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/other/purge_filetypes.txt'

def is_json(myjson):
    try: json_object = json.loads(myjson)
    except ValueError as e: return False
    return True

def cleanup(suffix_file, sample_dir):
    try: assert(os.path.isdir(sample_dir))
    except: sys.stderr.write("ERROR: Sample directory doesn't seem to exist! Exiting now ...\n"); raise RuntimeError
    try: assert(not suffix_file or os.path.isfile(suffix_file))
    except: sys.stderr.write("ERROR: Input file for purging filetypes doesn't exist. Exiting now ...\n"); raise RuntimeError

    sample_dir = os.path.abspath(sample_dir) + '/'

    if not suffix_file:
        suffix_file = DEFAULT_SUFFIX_PURGE_FILE

    # create logging object
    log_file = sample_dir + 'CLEANUP.txt'
    logObject = uF.createLoggerObject(log_file)

    sample = sample_dir.split('/')[-3]
    logObject.info("Cleanup initiating for sample %s", sample)
    logObject.info("-"*80)
    purging_filetypes = defaultdict(set)
    dot_depth = 0
    with open(suffix_file) as osf:
        for line in osf:
            line = line.rstrip('\n').strip()
            ls = line.split()
            grouping = ls[0]
            for l in ls[1:]:
                purging_filetypes[grouping].add(l)
                type_dot_depth = len(l.split('.'))
                if type_dot_depth > dot_depth:
                    dot_depth = type_dot_depth

    for subdir, dirs, files in os.walk(sample_dir):
        for file in files:
            file_ends = set([])
            file_suffices_len = len(file.split('.'))-1
            neg_counter = -1
            for i in range(0, dot_depth):
                if i < file_suffices_len:
                    if i == 0:
                        file_ends.add('.' + file.split('.')[neg_counter])
                    else:
                        file_ends.add('.' + '.'.join(file.split('.')[neg_counter:]))
                    neg_counter -= 1
            for group in purging_filetypes.items():
                if len(file_ends.intersection(group[1])) > 0 and (group[0].lower().strip() == "all" or group[0].lower().strip() == os.path.abspath(subdir).split()[-1].lower().strip()):
                    full_file_path = os.path.join(subdir, file)
                    if not '/Input/' in full_file_path and not '/KneadData/' in full_file_path:
                        logObject.info('purging file %s' % full_file_path)
                        os.system("rm -f " + full_file_path)

    uF.closeLoggerObject(logObject)

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    This module cleans up and purges any heavy item files ending with specific suffices 
    the user can provide in a file.
    
    The purge descriptions file should be space/tab separated with the first column regarding the subdirectory, 
    e.g. FastQC or Trimmomatic, while the rest regard suffices. The first column is not case sensitive.
    Additionally, "all" can be used for the first column for instance if you want to remove all files.
    """)

    parser.add_argument('-s', '--suffix_file', help='path to file with suffices of file types to remove. Default is: /path/to/seQc/other/purge_filetypes.txt', required=False, default=None)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    args = parser.parse_args()

    cleanup(args.suffix_file, args.sample_dir)