from setuptools import setup
import os
import sys

def generate_executable_wrappers(conda_env_file, prog_to_conda, activate_path, deactivate_path, external_wrappers_dir, other_dir):
    try:
        assert (os.path.isdir(external_wrappers_dir) and os.path.isdir(other_dir) and os.path.isfile(conda_env_file) and os.path.isfile(prog_to_conda))
    except:
        raise RuntimeError("Setup issues which shouldn't happen assuming proper download of git repo.")

    external_wrappers_dir = os.path.abspath(external_wrappers_dir) + '/'
    env_paths = {}
    with open(conda_env_file) as ocef:
        for line in ocef:
            line = line.strip()
            print(line.split('\t'))
            env_name, env_path = line.split('\t')[:2]
            env_paths[env_name] = env_path

    # create main environment loading file
    outf = open(other_dir + 'activate_main_env.sh', 'w')
    outf.write('source %s %s' % (activate_path, env_paths["main_env"]))
    outf.close()

    # create main environment unloading file
    outf = open(other_dir + 'deactivate_main_env.sh', 'w')
    outf.write('source %s %s' % (deactivate_path, env_paths["main_env"]))
    outf.close()

    with open(prog_to_conda) as optc:
        for line in optc:
            if line.startswith("#"): continue
            line = line.strip()
            print(line.split('\t'))
            prog, env_name, prog_version = line.split('\t')

            env_path = env_paths[env_name]

            prog_name = '_'.join(prog.split()).replace('-', '_').split('.py')[0].split('.R')[0].split('.sh')[0].lower()
            prog_shell_script = external_wrappers_dir + prog_name + '.sh'
            prog_version_text = external_wrappers_dir + prog_name + '.txt'

            pss_handle = open(prog_shell_script, 'w')
            pss_handle.write('#!/usr/bin/env bash\n\n')
            pss_handle.write('source %s %s\n' % (activate_path, env_path))
            if prog_name == 'kneaddata':
                pss_handle.write('%s --trimmomatic %s $@\n' % (prog, env_paths["main_env"] + 'share/trimmomatic-*/'))
            else:
                pss_handle.write('%s $@\n' % prog)
            pss_handle.write('source %s\n' % (deactivate_path))
            pss_handle.close()

            pvt_handle = open(prog_version_text, 'w')
            pvt_handle.write('%s (%s)' % (prog_version, env_name))
            pvt_handle.close()

conda_inst_file = os.path.abspath('./path_to_conda_installation.sh')
conda_envs_file = os.path.abspath('./conda_environments.txt')
prog_to_cenv_file = os.path.abspath('./prog_to_environment.txt')
external_wrappers_dir = os.path.abspath('./seQuoia/external_wrappers') + '/'
other_dir = os.path.abspath('./seQuoia/other') + '/'

if not os.path.isdir(external_wrappers_dir):
	os.system('mkdir %s' % external_wrappers_dir)

### Parse conda installation path
conda_dir = None; activate_path = None; deactivate_path = None
with open(conda_inst_file) as ocif:
    for line in ocif:
        conda_dir = os.path.abspath(line.strip()) + '/'
        activate_path = conda_dir + 'bin/activate'
        deactivate_path = conda_dir + 'bin/deactivate'

try:
    assert(os.path.isdir(conda_dir) and os.path.isfile(activate_path) and os.path.isfile(deactivate_path))
except:
    raise RuntimeError

# create executable shell scripts
generate_executable_wrappers(conda_envs_file, prog_to_cenv_file, activate_path, deactivate_path, external_wrappers_dir, other_dir)

try:
    setup(name='seQuoia', version='0.1', description='Workflow framework for bacterial genomics / metagenomics.', url='https://github.com/broadinstitute/seQuoia', author='Rauf Salamzade', author_email='rsalamza@broadinstitute.org', license='GNU 3', packages=['seQuoia'], zip_safe=False)
except:
    raise RuntimeError("Error with setup! Please email/raise issue with author.")
