import os
import sys
import argparse
from seQuoia.other import usefulFunctions as uF

def reorganize(sample_dir, reorganize, sample_name):
    try:
        assert (os.path.isdir(sample_dir))
    except:
        sys.stderr.write("ERROR: Sample directory doesn't seem to exist! Exiting now ...\n"); raise RuntimeError

    sample_dir = os.path.abspath(sample_dir) + '/'

    if reorganize:
        # set up directory structure
        workspace_name = "Assembly_Results/"
        workspace = uF.setupDirectory(sample_dir, workspace_name, panic_if_exists=False)

        # create logging object
        log_file = workspace + 'Reorganization.log'
        logObject = uF.createLoggerObject(log_file)

        logObject.info("Creating easy upload formats for sample %s", sample_name)
        logObject.info("-" * 80)

        runs = ['full-np', 'sub-np', 'canu']
        names = ['Unicycler_All-ONT', 'Unicycler_Subsampled-ONT', 'Canu_Pure-ONT']
        for i, run in enumerate(runs):
            run_name = names[i]

            # De Novo Assembly + QC Storage

            logObject.info('Moving GAEMR folder for run %s to results directory.' % run)
            logObject.info('-' * 80)

            try:
                new_location = os.path.abspath(workspace + run_name)
                original_dir = os.path.abspath(sample_dir + 'GAEMR_' + run) + '/'
                assert(os.path.isdir(original_dir))
                os.system('mv %s %s' % (original_dir, new_location))
            except:
                logObject.warning('Unable to move GAEMR directory for run %s to results directory.' % run_name)

            # MLST Results Storage

            logObject.info('Moving MLST folder for run %s to results directory.' % run)
            logObject.info('-' * 80)

            try:
                new_location = os.path.abspath(workspace + run_name)
                original_dir = os.path.abspath(sample_dir + 'Assembly_MLST_' + run) + '/'
                assert (os.path.isdir(original_dir))
                os.system('mv %s %s' % (original_dir, new_location))
            except:
                logObject.warning('Unable to move Assembly_MLST directory for run %s to results directory.' % run_name)

        logObject.info('*' * 80)

        intermediate_workspace_name = 'Intermediate_Subdirectories/'
        intermediate_workspace = uF.setupDirectory(sample_dir, intermediate_workspace_name)

        logObject.info("Moving intermediate subdirectories of workflow to directory %s" % intermediate_workspace)

        try:
            for sub in os.listdir(sample_dir):
                sub_dir = os.path.abspath(sample_dir + sub) + '/'
                if os.path.isdir(sub_dir) and sub != 'Assembly_Results':
                    os.system('mv %s %s' % (sub_dir, intermediate_workspace))
        except:
            logObject.error("Something went wrong when moving intermediate directories.")
            raise RuntimeError()

        logObject.info('*' * 80)

        checkpoint_workspace_name = 'Checkpoint_Files/'
        checkpoint_workspace = uF.setupDirectory(sample_dir, checkpoint_workspace_name)

        logObject.info("Moving checkpoint files of workflow to directory %s" % checkpoint_workspace)

        try:
            for f in os.listdir(sample_dir):
                checkpoint_file = os.path.abspath(sample_dir + f)
                if os.path.isfile(checkpoint_file) and f.endswith('.txt'):
                    os.system('mv %s %s' % (checkpoint_file, checkpoint_workspace))
        except:
            logObject.error("Something went wrong when moving checkpoint files.")
            raise RuntimeError()

        uF.closeLoggerObject(logObject)

        # create successful reorganization file if steps completed!
        conf_file = open(sample_dir + "REORGANIZATION.txt", 'w')
        conf_file.write("Reorganization was Completed Successfully!")
        conf_file.close()

    # create successful completion file if steps completed!
    conf_file = open(sample_dir + "COMPLETION.txt", 'w')
    conf_file.write("Completion: Module Completed Successfully!")
    conf_file.close()


if __name__ == '__main__':
    #Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Module to wrap up the Bacterial Assembly Workflow.
    """)

    parser.add_argument('-o', '--sample_dir', help='path to the output directory.', required=True)
    parser.add_argument('-r', '--reorganize', action='store_true', help='reorganize directory structure to make more human readable?', required=False, default=False)
    parser.add_argument('-s', '--sample', help='the sample name', required=True)
    args = parser.parse_args()

    reorganize(args.sample_dir, args.reorganize, args.sample)