import os
import sys
import argparse
from seQuoia.classes.AssemblyAnalyzer import Assembly
from seQuoia.other import usefulFunctions as uF

def runMLST(assembly, sample_name, parent_dir, identifier):
    try: assert( os.path.isfile(assembly) )
    except: sys.stderr.write("ERROR: Assembly does not seem to exist.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Assembly_MLST"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Assembly_MLST.log'
    logObject = uF.createLoggerObject(log_file)

    # Initialize Assembly Object
    AssemblyObj = Assembly(assembly, sample_name, logObject)

    # Run MLST
    AssemblyObj.run_mlst(workspace)

    # create successful completion file if steps completed!
    conf_file_name = parent_dir + "ASSEMBLY_MLST"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("Assembly MLST: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module for running multi-locus subtyping analysis. Uses FastqAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-a', '--assembly', help='location of the assembly in FASTA format.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-i', '--identifier', help='identifier, in case multiple runs in workflow.', required=False, default=None)
    args = parser.parse_args()

    runMLST(args.assembly, args.sample, args.sample_dir, args.identifier)