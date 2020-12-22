import os
import sys
import subprocess
import time
import signal

try:
    assert (sys.version_info[0] == 3)
except:
    sys.stderr.write("Please use Python-3.4 to run this program. Exiting now ...\n");
    sys.exit(1)

other_dir = '/'.join(
    os.path.dirname(os.path.realpath(__file__)).split('/')[:-1] + ['other']) + '/'

def is_number(x):
    try: float(x); return True
    except: return False

class hUGE:
    def __init__(self, commands, pbs_script, log_file, platform='UGE',environment_variables=[], memory=1, cpus=1, throttle=False, runtime='48:00:00', project='gscid'):
        self.commands = commands
        self.pbs_script = pbs_script
        self.environment_variables = environment_variables
        self.platform = platform
        self.log_file = log_file
        self.throttle = throttle
        self.runtime = runtime
        self.memory = memory
        self.cpus = cpus
        self.project = project

    def run_cmd(self):

        if self.platform == 'UGE':

            pbs_handle = open(self.pbs_script, 'w')
            pbs_handle.write("#!/bin/bash\n")
            pbs_handle.write('#$ -cwd\n')
            pbs_handle.write('#$ -l h_rt=%s\n' % self.runtime)
            pbs_handle.write('#$ -j y\n')
            pbs_handle.write('#$ -o %s\n' % self.log_file)
            if self.cpus != 1:
                try:
                    assert(is_number(self.cpus))
                    pbs_handle.write('#$ -pe smp %s\n#$ -R y\n' % (str(self.cpus)))
                except:
                    sys.stderr.write('Cores specified was non-numeric, using just a single core.')
            if self.memory != 1:
                try:
                    assert(is_number(self.memory))
                    memory_per_core = self.memory
                    memory_full = memory_per_core * self.cpus
                    pbs_handle.write('#$ -l h_vmem=%sG\n' % str(int(memory_full)))
                except:
                    sys.stderr.write('Memory specified was non-numeric, using just a single gigabyte.')

            pbs_handle.write('hostname\n')
            pbs_handle.write('ENV_FILE="' + other_dir + 'activate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE)\n')
            pbs_handle.write('eval $cmd_env\n')
            pbs_handle.write(self.commands + '\n')
            pbs_handle.write('ENV_FILE_UNLOAD="' + other_dir + 'deactivate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE_UNLOAD)\n')
            pbs_handle.close()

            qsub_cmd = ['qsub', self.pbs_script]
            qsub_submission_success_flag = False
            qsub_submission_id = None
            while not qsub_submission_success_flag:
                proc = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                time.sleep(10)
                if out.decode('utf-8').strip() == "": continue
                qsub_submission_id = out.decode('utf-8').split()[2].split('.')[0]
                if is_number(qsub_submission_id): qsub_submission_success_flag = True

            os.system("echo '%s\n' >> %s" % (qsub_submission_id, self.log_file.split('.log')[0] + '.job-id'))
            time.sleep(10)

            flag = True
            qstat_cmd = ['qstat']

            while (flag):
                proc = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                loc_flag = False
                if out.decode('utf-8').strip() == "": continue
                for i, line in enumerate(out.decode('utf-8').split('\n')):
                    line = line.strip()
                    ls = line.split()
                    if i > 1 and len(ls) > 4:
                        if ls[0] == qsub_submission_id and not 'E' in ls[4]:
                            loc_flag = True
                if not loc_flag:
                    flag = False
                time.sleep(30)

        elif self.platform == 'SLURM':
            try:
                assert (is_number(self.memory) and is_number(self.cpus))
            except:
                sys.stderr.write('Memory and/or CPU specified was non-numeric.')
                raise RuntimeError

            pbs_handle = open(self.pbs_script, 'w')
            pbs_handle.write("#!/bin/bash\n")
            pbs_handle.write('#SBATCH --cpus-per-task=%s\n' % str(self.cpus))
            pbs_handle.write('#SBATCH --mem=%sGB\n' % str(self.memory))
            pbs_handle.write('#SBATCH -t %s\n' % self.runtime)
            pbs_handle.write('#SBATCH -o %s\n' % self.log_file)
            pbs_handle.write('#SBATCH -e %s\n' % self.log_file)
            pbs_handle.write('hostname\n')
            pbs_handle.write('ENV_FILE="' + other_dir + 'activate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE)\n')
            pbs_handle.write('eval $cmd_env\n')
            pbs_handle.write(self.commands + '\n')
            pbs_handle.write('ENV_FILE_UNLOAD="' + other_dir + 'deactivate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE_UNLOAD)\n')
            pbs_handle.close()

            slurm_cmd = ['sbatch', self.pbs_script]
            qsub_submission_id = None
            qsub_submission_success_flag = False
            while not qsub_submission_success_flag:
                proc = subprocess.Popen(slurm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                time.sleep(10)
                if out.decode('utf-8').strip() == "": continue
                qsub_submission_id = out.decode('utf-8').split()[-1].strip()
                if is_number(qsub_submission_id): qsub_submission_success_flag = True

            os.system("echo '%s\n' >> %s" % (qsub_submission_id, self.log_file.split('.log')[0] + '.job-id'))
            time.sleep(10)

            flag=True
            squeue_cmd = ['squeue']
            while (flag):
                proc = subprocess.Popen(squeue_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                loc_flag = False
                for i, line in enumerate(out.decode('utf-8').split('\n')):
                    line = line.strip()
                    ls = line.split()
                    if i > 0 and len(ls) > 4:
                        if ls[0] == qsub_submission_id and not 'E' in ls[4]:
                            loc_flag = True
                if not loc_flag:
                    flag = False
                time.sleep(30)
        elif self.platform == 'Local':
            activate_cmd = open(other_dir + 'activate_main_env.sh').readlines()[0].split()
            deactivate_cmd = open(other_dir + 'deactivate_main_env.sh').readlines()[0].split()
            all_cmds = activate_cmd + [';'] + self.commands.split() + [';'] + deactivate_cmd
            log_handle = open(self.log_file, 'w')
            proc = subprocess.Popen(' '.join(all_cmds), stdout=log_handle, stderr=log_handle, shell=True)
            proc.wait()
            log_handle.close()

    def run_cmd_no_wait(self):

        if self.platform == 'UGE':

            pbs_handle = open(self.pbs_script, 'w')
            pbs_handle.write("#!/bin/bash\n")
            pbs_handle.write('#$ -cwd\n')
            pbs_handle.write('#$ -l h_rt=%s\n' % self.runtime)
            pbs_handle.write('#$ -j y\n')
            pbs_handle.write('#$ -o %s\n' % self.log_file)
            if self.cpus != 1:
                try:
                    assert(is_number(self.cpus))
                    pbs_handle.write('#$ -pe smp %s\n#$ -R y\n' % (str(self.cpus)))
                except:
                    sys.stderr.write('Cores specified was non-numeric, using just a single core.')
            if self.memory != 1:
                try:
                    assert(is_number(self.memory))
                    memory_per_core = self.memory
                    memory_full = memory_per_core * self.cpus
                    pbs_handle.write('#$ -l h_vmem=%sG\n' % str(int(memory_full)))
                except:
                    sys.stderr.write('Memory specified was non-numeric, using just a single gigabyte.')

            pbs_handle.write('hostname\n')
            pbs_handle.write('ENV_FILE="' + other_dir + 'activate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE)\n')
            pbs_handle.write('eval $cmd_env\n')
            pbs_handle.write(self.commands + '\n')
            pbs_handle.write('ENV_FILE_UNLOAD="' + other_dir + 'deactivate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE_UNLOAD)\n')
            pbs_handle.close()

            qsub_cmd = ['qsub', self.pbs_script]
            qsub_submission_success_flag = False
            qsub_submission_id = None
            while not qsub_submission_success_flag:
                proc = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                time.sleep(10)
                if out.decode('utf-8').strip() == "": continue
                qsub_submission_id = out.decode('utf-8').split()[2].split('.')[0]
                if is_number(qsub_submission_id): qsub_submission_success_flag = True

            os.system("echo '%s\n' >> %s" % (qsub_submission_id, self.log_file.split('.log')[0] + '.job-id'))
            time.sleep(10)

            return qsub_submission_id
        elif self.platform == 'SLURM':
            try:
                assert (is_number(self.memory) and is_number(self.cpus))
            except:
                sys.stderr.write('Memory and/or CPU specified was non-numeric.')
                raise RuntimeError

            pbs_handle = open(self.pbs_script, 'w')
            pbs_handle.write("#!/bin/bash\n")
            pbs_handle.write('#SBATCH --cpus-per-task=%s\n' % str(self.cpus))
            pbs_handle.write('#SBATCH --mem=%sGB\n' % str(self.memory))
            pbs_handle.write('#SBATCH -t %s\n' % self.runtime)
            pbs_handle.write('#SBATCH -o %s\n' % self.log_file)
            pbs_handle.write('#SBATCH -e %s\n' % self.log_file)
            pbs_handle.write('hostname\n')
            pbs_handle.write('ENV_FILE="' + other_dir + 'activate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE)\n')
            pbs_handle.write('eval $cmd_env\n')
            pbs_handle.write(self.commands + '\n')
            pbs_handle.write('ENV_FILE_UNLOAD="' + other_dir + 'deactivate_main_env.sh"\n')
            pbs_handle.write('cmd_env=$(head -n 1 $ENV_FILE_UNLOAD)\n')
            pbs_handle.close()

            slurm_cmd = ['sbatch', self.pbs_script]
            qsub_submission_success_flag = False
            qsub_submission_id = None
            while not qsub_submission_success_flag:
                proc = subprocess.Popen(slurm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (out, err) = proc.communicate()
                time.sleep(10)
                if out.decode('utf-8').strip() == "": continue
                qsub_submission_id = out.decode('utf-8').split()[-1].strip()
                if is_number(qsub_submission_id): qsub_submission_success_flag = True

            os.system("echo '%s\n' >> %s" % (qsub_submission_id, self.log_file.split('.log')[0] + '.job-id'))
            time.sleep(10)
            return qsub_submission_id
        elif self.platform == 'Local':
            activate_cmd = open(other_dir + 'activate_main_env.sh').readlines()[0].split()
            deactivate_cmd = open(other_dir + 'deactivate_main_env.sh').readlines()[0].split()
            all_cmds = activate_cmd + [';'] + self.commands.split() + [';'] + deactivate_cmd
            log_handle = open(self.log_file, 'w')
            proc = subprocess.Popen(' '.join(all_cmds), stdout=log_handle, stderr=log_handle, shell=True)
            proc.wait()
            log_handle.close()