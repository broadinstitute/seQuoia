import os
import sys
import argparse
from seQuoia.classes.AssemblyAnalyzer import Assembly
from seQuoia.other import usefulFunctions as uF

def runAssemblyAdapterRemoval(assembly, sample_name, parent_dir, run_guinan, gaemr_options, guinan_options, size_filter, identifier):
    try: assert( os.path.isfile(assembly) )
    except: sys.stderr.write("ERROR: Assembly does not seem to exist.\n"); raise RuntimeError

    # set up directory structure
    workspace_name = "Assembly_Adapter_Removal"
    if identifier:
        workspace_name += '_' + identifier
    workspace = uF.setupDirectory(parent_dir, workspace_name)

    # create logging object
    log_file = workspace + 'Assembly_Adapter_Removal.log'
    logObject = uF.createLoggerObject(log_file)

    # Initialize Assembly Object
    AssemblyObj = Assembly(assembly, sample_name, logObject)

    # Run GAEMR formatting program to generate assembly graph
    workspace_a = uF.setupDirectory(workspace, "Assembly_Formatted/")
    AssemblyObj.run_gaemr_formatter(workspace_a, reference_change=True)

    if run_guinan:
        # Run GAEMR based adapters in assembly
        guinan_commands_file = AssemblyObj.detect_adapters(workspace, options=gaemr_options)

        # Run guinan suite to remove detected adapters from assembly
        AssemblyObj.remove_adapters(guinan_commands_file, workspace, options=guinan_options)

    # Run assembly filter by contig size
    AssemblyObj.filter_contigs_by_size(workspace, size_filter=size_filter)

    # create successful completion file if steps completed!
    conf_file_name = parent_dir + "ASSEMBLY_ADAPTER_REMOVAL"
    if identifier: conf_file_name += '_' + identifier
    conf_file_name += ".txt"
    conf_file = open(conf_file_name, 'w')
    conf_file.write("Assembly Adapter Removal: Module Completed Succesfully!")
    conf_file.close()

if __name__ == '__main__':
    #Pull out the arguments.

    parser = argparse.ArgumentParser(description="""
    seQc module for running adapter screening of assemblies and subsequent removal. Also filter assemblies on contig size. Uses AssemblyAnalyzer.py class OOP framework.
    """)

    parser.add_argument('-a', '--assembly', help='location of the assembly in FASTA format.', required=True)
    parser.add_argument('-s', '--sample', help='sample name.', required=True)
    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-r', '--run_guinan', action='store_true', help='run adapter trimming of assembly using guinan suite. Currently only available on Broad servers.', required=False, default=False)
    parser.add_argument('-g', '--gaemr_options', help='options for GAEMR run_contamination_screen.py.', default="-r -e -f Directive -v --no_protection --adaptorsdb", required=False)
    parser.add_argument('-u', '--guinan_options', help='path to the output directory.', required=False, default="--sort_by_length --rename")
    parser.add_argument('-f', '--size_filter', type=int, help='contig size filter.', required=False, default=0)
    parser.add_argument('-i', '--identifier', help='identifier, in case multiple runs in workflow.', required=False, default=None)
    args = parser.parse_args()

    runAssemblyAdapterRemoval(args.assembly, args.sample, args.sample_dir, args.run_guinan, args.gaemr_options, args.guinan_options, args.size_filter, args.identifier)