import os
import sys
import datetime

# PLEASE DO NOT RUN IN THE CIL SHED INSTANCE OF SEQC AND ONLY IN LOCAL GIT REPOSITORIES
# FROM WHICH YOU PUSH TO THE SEQC GIT REPO. CIL SHED INSTANCE MERELY FOR PULLING!

main_directory = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/'
ep_directory = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')) + '/'
Rscript_program = ep_directory + 'Rlist.R'
miniconda_activate_path = open(main_directory + 'path_to_conda_installation.sh').readlines()[0].strip() + 'bin/activate'
conda_env_listings_file = main_directory + 'conda_environments.txt'

currentDT = datetime.datetime.now()
date_time_id = (currentDT.strftime("%Y-%m-%d_%H-%M-%S"))
curr_directory = ep_directory + 'archive/' + date_time_id + '/'
if not os.path.isdir(curr_directory):
	os.system('mkdir %s' % curr_directory)
print('Storing current versioining info in the directory :%s' % curr_directory) 

with open(conda_env_listings_file) as ocelf:
	for line in ocelf:
		if line.startswith('#'): continue
		line = line.strip()
		ls = line.split('\t')
		env_name, env_location = ls[:2]
		print('Gathering versioning info for %s' % env_name)
		conda_yml_outf = curr_directory + env_name + '.conda-yml.txt'
		conda_yml_agn_outf = curr_directory + env_name + '.conda-cross-platform-yml.txt'
		conda_explicit_outf = curr_directory + env_name + '.conda-explicit-list.txt'
		conda_list_outf = curr_directory + env_name + '.conda-list.txt'
		pip_freeze_outf = curr_directory + env_name + '.pip-freeze.txt'
		R_session_info_outf = curr_directory + env_name + '.R-libraries.txt'
		load_cmd = 'source %s %s' % (miniconda_activate_path, env_location)
		conda_yml_cmd = 'conda env export --prefix %s > %s' % (env_location, conda_yml_outf)
		conda_yml_agn_cmd = 'conda env export --prefix %s --from-history > %s' % (env_location, conda_yml_agn_outf)
		conda_exp_list_cmd = 'conda list --explicit > %s' % conda_explicit_outf
		conda_list_cmd = 'conda list > %s' % conda_list_outf
		pip_freeze_cmd = 'pip freeze --all > %s' % pip_freeze_outf
		r_sess_list_cmd = ''
		if len(ls) == 3 and ls[2].upper().strip() == 'R':
			r_sess_list_cmd = 'Rscript %s %s' % (Rscript_program, R_session_info_outf)
		os.system('; '.join([load_cmd, conda_yml_cmd, conda_yml_agn_cmd, conda_exp_list_cmd, conda_list_cmd, pip_freeze_cmd, r_sess_list_cmd]))

curr_ref_directory = ep_directory + 'current/'
os.system('rm -rf %s' % curr_ref_directory)
os.system('cp -r %s %s' % (curr_directory, curr_ref_directory))
datetime_info_file = open(curr_ref_directory + 'datetime_info.txt', 'w')
datetime_info_file.write(date_time_id + '\n')
datetime_info_file.close()
