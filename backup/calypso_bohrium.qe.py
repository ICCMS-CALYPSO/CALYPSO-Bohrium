#/usr/bin/env python
import dpdispatcher
from dpdispatcher import Machine, Resources, Task, Submission
import os
import sys
import shutil
import glob
import oss2
import sys
from pathlib import Path


def get_value(key):
    temp = (os.popen(f'grep {key} input.dat').read().split('=')[-1])
    return temp 

def get_pp(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    for idx, line in enumerate(lines):
        if 'pseudo_dir' in line:
            pp_dir = line.split('=')[-1].strip().strip(',').strip().strip("'").strip(',')
        if 'ntyp' in line:
            ntyp = int(line.split('=')[-1].strip().strip(',').strip())
        if 'ATOMIC_SPECIES' in line:
            pp_name = []
            for i in range(ntyp):
                pp_name.append(lines[idx+i+1].strip().strip(',').strip("'").strip('"').split()[-1].strip())
            # print(pp_name)
    return (pp_dir, pp_name)

qe='''Nat=`grep 'number of k points' -B 2 out.pw | head -n 1 | awk {'print($1)'}`;StruLine=`expr $Nat + 5 `; grep 'CELL_' -A $StruLine out.pw | tail -n `expr $StruLine + 1` > pwscf'''

machine = Machine.load_from_json('machine.json')
resources = Resources.load_from_json('resources.json')

_ncpu = machine.input_data['machine_type'].split('_')[0].strip('c')
ncpu = _ncpu if _ncpu != '' else 32

MaxStep = int(get_value('MaxStep'))
PopSize = int(get_value('PopSize'))
N_INCAR = int(get_value('NumberOfLocalOptim'))

if get_value('PickUp') == 'T':
    PickStep = int(get_value('PickStep'))
    os.system(f'echo {PickStep+1} > step')    # generate (pickstep+1)th structures
    os.system('touch pickup')
else:
    if os.path.exists('step'):
        PickStep = int(os.popen('cat step').read())
    else:
        PickStep = 1

os.system("echo 'Split = T' >> ./input.dat") if get_value('Split') == '' else ''

command_intel = 'source /opt/intel/oneapi/setvars.sh;'
# prepare a mapping of dft run command 
command_runvasp_list =  ['cp pwscf.ori pwscf']
command_runvasp_list += [f'cat pw_input_{idx} pwscf >> pw_input; mpirun -n {ncpu} pw.x < pw_input > out.pw 2>&1 ; {qe}' for idx in range(1, N_INCAR+1)]  # cpu number how to detect
command_runvasp = ';'.join(command_runvasp_list)

command = command_intel + command_runvasp

Path('./log_dir').mkdir(parents=True, exist_ok=True)

for step in range(PickStep, MaxStep + 1):

    if not os.path.exists('pickup') or (os.path.exists('pickup') and os.path.exists('restart')):
        os.system("./calypso.x | tee caly.log")

    task_list = []
    for pop in range(1, PopSize + 1):
        task_dir = "./data/step%03d.pop%03d"% (step,pop)
        Path(task_dir).mkdir(parents=True, exist_ok=True)

        if not os.path.exists('pickup') or (os.path.exists('pickup') and os.path.exists('restart')):
            # depand on the type of dft software
            shutil.copyfile("pwscf_%d"%pop, os.path.join(task_dir, "pwscf"))
            shutil.copyfile("pwscf_%d"%pop, os.path.join(task_dir, "pwscf.ori"))
            for n_incar in range(1, N_INCAR + 1):
                shutil.copyfile("pw_input_%d"%n_incar , os.path.join(task_dir, "pw_input_%d"%n_incar))
                pp_dir, _pp_name = get_pp('pw_input_1')
                for pp_name in _pp_name:
                    shutil.copyfile(pp_dir + '/' + pp_name, os.path.join(task_dir, pp_name))

        task = Task(
            command = command,
            task_work_path = task_dir,
            forward_files = ['pwscf', 'pwscf.ori'] + [f'pw_input_{idx}' for idx in range(1, N_INCAR + 1)] + [p_name for p_name in _pp_name],
            backward_files = ['out.pw', 'pwscf', 'log', 'err']
        )
        task_list.append(task)
    # print(['pwscf', 'pwscf.ori'] + [f'pw_input_{idx}' for idx in range(1, N_INCAR + 1)] + [p_name for p_name in _pp_name], )
    submission = Submission(work_base = os.getcwd(), machine= machine, resources =resources, task_list = task_list)
    submission.run_submission()
    if os.path.exists('pickup'):
        os.system('mv pickup pickup_done')
    if os.path.exists('restart'):
        os.system('mv restart restart_done')

    for pop in range(1, PopSize + 1):
        task_dir = "./data/step%03d.pop%03d"% (step,pop)
        shutil.copyfile(os.path.join(task_dir, "out.pw"), "out.pw_%d"%pop )
        shutil.copyfile(os.path.join(task_dir, "pwscf.ori"), "pwscf.%d"%pop )

    os.system('mv *.sub lbg-*.sh *_fail *_finished log_dir/')
                
os.system("./calypso.x | tee caly.log")
