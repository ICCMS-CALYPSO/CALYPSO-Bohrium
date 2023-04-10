#/usr/bin/env python
import dpdispatcher
from dpdispatcher import Machine, Resources, Task, Submission
import os
import sys
import shutil
import glob
import oss2


def get_value(key):
    temp = (os.popen(f'grep {key} input.dat').read().split('=')[-1])
    return temp 

machine = Machine.load_from_json('machine.json')
resources = Resources.load_from_json('resources.json')

_ncpu = d['remote_profile']['input_data']['scass_type'].split('_')[0].strip('c')
ncpu = _ncpu if _ncpu != '' else 32

MaxStep = int(get_value('MaxStep'))
PopSize = int(get_value('PopSize'))
N_INCAR = int(get_value('NumberOfLOcalOptim'))

if get_value('PickUp') == 'T':
    PickStep = int(get_value('PickStep'))
    os.system(f'echo {PickStep} > step')
else:
    if os.path.exists('step'):
        PickStep = int(os.popen('cat step').read())
    else:
        PickStep = 1

os.system("echo 'Split = T' >> ./input.dat") if get_value('Split') == '' else pass

command_intel = 'source /opt/intel/oneapi/setvars.sh;'
# prepare a mapping of dft run command 
command_runvasp_list = [f'cp INCAR_{idx} INCAR; mpirun -n {ncpu} vasp_std 1 >> vasp.log 2 >> err' for idx in range(1, N_INCAR+1)]  # cpu number how to detect
command_runvasp = ';'.join(command_runvasp_list)

command = command_intel + command_runvasp


for step in range(PickStep, MaxStep + 1):

    os.system("./calypso.x | tee caly.log")

    task_list = []
    for pop in range(1, PopSize + 1):
        task_dir = "step%03d.pop%03d"% (step,pop)
        os.mkdir(task_dir)

        # depand on the type of dft software
        shutil.copyfile("POSCAR_%d"%pop, os.path.join(task_dir, "POSCAR"))
        shutil.copyfile("POSCAR_%d"%pop, os.path.join(task_dir, "POSCAR.ori"))
        for n_incar in range(1, N_INCAR + 1):
            shutil.copyfile("INCAR_%d"%n_incar , os.path.join(task_dir, "INCAR_%d"%n_incar))
            shutil.copyfile("POTCAR", os.path.join(task_dir, "POTCAR"))

        task = Task(
            command = command,
            task_work_path = task_dir,
            forward_files = ['POSCAR', 'INCAR_1', 'INCAR_2', 'INCAR_3', 'POTCAR'],
            backward_files = ['CONTCAR', 'OUTCAR', 'log', 'err']
        )
        task_list.append(task)
        submission = Submission(work_base = os.getcwd(), machine= machine, resources =resources, task_list = task_list)
        submission.run_submission()

        for pop in range(1, PopSize + 1):
            task_dir = "step%03d.pop%03d"% (step,pop)
            shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d"%pop )
            shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d"%pop)
                
os.system("./calypso.x | tee caly.log")

