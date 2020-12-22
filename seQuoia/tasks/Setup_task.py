import os
import sys
import argparse
import json
from seQuoia.other import usefulFunctions as uF

def is_json(myjson):
    try: json_object = json.loads(myjson)
    except ValueError as e: return False
    return True

def setup(sample_name, parent_dir, meta):
    meta = meta.strip('"')
    try: assert(os.path.isfile(meta) or is_json(meta.replace("QUOTES", "\"").replace("-COMMA-", ",")))
    except: sys.stderr.write("ERROR: Meta file does not exist! Exiting now ...\n"); raise RuntimeError

    meta = meta.replace("QUOTES", "\"").replace("-COMMA-", ",")

    # set up directory structure
    workspace = uF.setupDirectory(parent_dir, sample_name)
    meta_info_file = open(workspace + "meta_information.json", 'w')

    meta_information = {}
    if os.path.isfile(meta):
        with open(meta) as om:
                for line in om:
                    line = line.rstrip('\n')
                    key, value = line.split('\t')
                    meta_information[key] = value
        meta_info_file.write(json.dumps(meta_information))
    else:
        meta_info_file.write(meta)

    # create successful completion file if steps completed!
    conf_file = open(workspace + "SETUP.txt", 'w')
    conf_file.write("Setup: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    This module sets up the seQuoia repo for a sample. 
    """)

    parser.add_argument('-s', '--sample', help="Sample name.", required=True)
    parser.add_argument('-o', '--parent_dir', help='path to the output directory.', required=True)
    parser.add_argument('-m', '--meta', help='Meta information file/json-string.', required=True)
    args = parser.parse_args()

    setup(args.sample, args.parent_dir, args.meta)