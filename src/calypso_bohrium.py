#/usr/bin/env python
import dpdispatcher
from dpdispatcher import Machine, Resources, Task, Submission
import os
import sys
import shutil
import glob
import oss2


# os.system("echo 'Split = T' >> /root/run_test/Test_bak/input.dat")

def get_value(key):
    temp = (os.popen(f'grep {key} input.dat').read().split('=')[-1])
    return temp 

MaxStep = int(get_value('MaxStep'))
PopSize = int(get_value('PopSize'))
N_INCAR = int(get_value('NumberOfLOcalOptim'))

machine = Machine.load_from_json('machine.json')
resources = Resources.load_from_json('resources.json')

command = '''source /opt/intel/oneapi/setvars.sh;
cp INCAR_1 INCAR;
mpirun -n 16 vasp_std 1>>log 2>>err;
'''

for step in range(1, MaxStep+1):

        os.system("/root/run_test/calypso5_test/calypso.x | tee caly.log")

            task_list = []
            for pop in range(1, PopSize + 1):
                task_dir = "step%03d.pop%03d"% (step,pop)
                os.mkdir(task_dir)
                shutil.copyfile("POSCAR_%d"%pop, os.path.join(task_dir, "POSCAR"))
                shutil.copyfile("POSCAR_%d"%pop, os.path.join(task_dir, "POSCAR.ori"))
                for n_incar in range(1, N_INCAR + 1):
                    shutil.copyfile("INCAR_%d"%n_incar , os.path.join(task_dir, "INCAR_%d"%n_incar))
                    shutil.copyfile("POTCAR", os.path.join(task_dir, "POTCAR"))

                    task = Task(command = command, task_work_path = task_dir, forward_files = ['POSCAR', 'INCAR_1', 'INCAR_2', 'INCAR_3', 'POTCAR'], backward_files = ['CONTCAR', 'OUTCAR', 'log', 'err'])
                    task_list.append(task)
                    submission = Submission(work_base = os.getcwd(), machine= machine, resources =resources, task_list = task_list)
                    submission.run_submission()

                    for pop in range(1, PopSize + 1):
                        task_dir = "step%03d.pop%03d"% (step,pop)
                        shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d"%pop )
                        shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d"%pop)
                        
                        os.system("/root/run_test/calypso5_test/calypso.x | tee caly.log")

os.system("/root/run_test/calypso5_test/calypso.x | tee caly.log")
