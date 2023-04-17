#!/usr/bin/env python
import os
import glob
import shutil

import dpdata
from dpdispatcher import Task

from calypso_bohrium.write_outcar import write_files


def read_abacus(path):
    is_success = True
    try:
        atoms_list = dpdata.LabeledSystem(path, 'abacus/relax').to_ase_structure()
    except Exception as e:
        print('read_abacus', e)
        stru_path = os.path.join(path, 'STRU')
        atoms_list = dpdata.System(stru_path, 'abacus/stru').to_ase_structure()
        is_success = False

    return (atoms_list[-1], is_success)

def read_stress(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    stress_list = []
    for line in lines:
        if line[0] == '#' or len(line.strip('\n')) == 0 or 'INPUT_PARAMETERS' in line:
            continue
        key, value = line.strip('\n').split()
        if key in ['press1', 'press2', 'press3']: 
            stress_list.append(float(value))
    return stress_list

def sort_abacus():
    stru_list = glob.glob('./OUT.*/STRU_ION*_D')
    t = stru_list.sort(key=lambda x: int(x.split('_')[1].strip('ION')))

    return t[-1]

def to_stru(name):
    with open('pp', 'r') as f:
        pp_file = f.readlines()
    pp = []
    for temp in pp_file:
        pp.append(temp)

    data = dpdata.System(name, 'vasp/poscar')
    n_species = len(data.data['atom_names'])
    data.to_abacus_stru('STRU')
    with open('STRU', 'r') as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if 'ATOMIC_SPECIES' in line:
            for ii in range(n_species):
                lines[idx+ii+1] = pp[ii]
            break
    with open('STRU', 'w') as f:
        f.write(''.join(lines))

def get_pp(filename):
    p_name = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if len(line.strip('\n').strip()) == 0:
            continue
        if 'ATOMIC_SPECIES' in line:
            continue
        if 'upf' in line.strip('\n').strip().split()[-1].lower():
            p_name.append(line.strip('\n').strip().split()[-1])
    return p_name

def abacus_command(N_INCAR, ncpu):

    python_command = "python -c 'import glob, os; t=glob.glob('./OUT.*/STRU_ION*_D');t.sort(key=lambda x: int(x.split('_')[1].strip('ION')));os.system(cp t[-1] STRU)'"
    pre_command = "OMP_NUM_THREADS=1;"
    if N_INCAR == 1:
        command_runvasp_list = [
            f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1"
            ]  # cpu number how to detect
        command_runvasp = ";".join(command_runvasp_list)
        return command_runvasp

    elif N_INCAR == 2:
        
        string = f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1; python continue.py;mkdir -p old; mv ./OUT.* ./old;"
        string += f"cp INPUT_2 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1;"
        return string

    elif N_INCAR == 3:
        string = f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1; python continue.py;mkdir -p old; mv ./OUT.* ./old;"
        string += f"cp INPUT_2 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1;python continue.py;mkdir -p old; mv ./OUT.* ./old;"
        string += f"cp INPUT_3 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1;"
        return string

def abacus_task(pop, task_dir, N_INCAR, command):

    _pp_name = get_pp('pp')
    with open('continue.py', 'w') as f:
        f.write("import glob, os\nt=glob.glob('./OUT.*/STRU_ION*_D')\nt.sort(key=lambda x: int(x.split('_')[1].strip('ION')))\nos.system(f'cp {t[-1]} STRU')\n")
    if not os.path.exists('pickup') or (os.path.exists('pickup') and os.path.exists('restart')):
        to_stru("POSCAR_%d" % pop)
        # os.system("cat pp stru > STRU")
        shutil.copyfile("STRU" , os.path.join(task_dir, "STRU"))
        shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR.ori"))
        shutil.copyfile("continue.py" , os.path.join(task_dir, "continue.py"))
        for n_incar in range(1, N_INCAR + 1):
            shutil.copyfile(
                "INPUT_%d" % n_incar, os.path.join(task_dir, "INPUT_%d" % n_incar)
            )
            for pp_name in _pp_name:
                shutil.copyfile('./' + pp_name, os.path.join(task_dir, pp_name))
    
    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=["STRU"] + [p_name for p_name in _pp_name] + ['continue.py']
        + [f"INPUT_{idx}" for idx in range(1, N_INCAR + 1)],
        # backward_files=["STRU", "OUTCAR", "log", "err"],
        backward_files=[],
    )
    return task

def abacus_back(task_dir, pop):
    atoms, is_success = read_abacus(task_dir)

    stress_list = read_stress(os.path.join(task_dir, './INPUT'))
    pstress = sum(stress_list)/len(stress_list) if len(stress_list) != 0 else 0.00001
    
    write_files(atoms, pstress, is_success, task_dir)

    shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d" % pop)
    shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d" % pop)


