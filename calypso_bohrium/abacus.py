#!/usr/bin/env python
import os
import shutil

import dpdata
from dpdispatcher import Task

from write_outcar import write_files


def read_abacus(path):
    is_success = True
    try:
        atoms_list = dpdata.LabeledSystem(path, 'abacus/relax').to_ase_structure()
    except Exception as e:
        print(e)
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

def abacus_command(N_INCAR, ncpu):

    command_runvasp_list = [
        f"cp INCAR_{idx} INCAR; mpirun -n {ncpu} vasp_std > fp.log 2>&1"
        for idx in range(1, N_INCAR + 1)
    ]  # cpu number how to detect
    command_runvasp = ";".join(command_runvasp_list)

    return command_runvasp

def abacus_task(pop, task_dir, N_INCAR, command):

   if not os.path.exists('pickup') or (os.path.exists('pickup') and os.path.exists('restart')):
       shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR"))
       shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR.ori"))
       for n_incar in range(1, N_INCAR + 1):
           shutil.copyfile(
               "INCAR_%d" % n_incar, os.path.join(task_dir, "INCAR_%d" % n_incar)
           )
           shutil.copyfile("POTCAR", os.path.join(task_dir, "POTCAR"))
   
   task = Task(
       command=command,
       task_work_path=task_dir,
       forward_files=["POSCAR", "POTCAR"]
       + [f"INCAR_{idx}" for idx in range(1, N_INCAR + 1)],
       backward_files=["CONTCAR", "OUTCAR", "log", "err"],
   )
   return task

def abacus_back(task_dir, pop):
    atoms, is_success = read_abacus('./')

    stress_list = read_stress('./INPUT')
    pstress = sum(stress_list)/len(stress_list) if len(stress_list) != 0 else 0.00001
    
    write_files(atoms, pstress, is_success)

    shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d" % pop)
    shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d" % pop)


