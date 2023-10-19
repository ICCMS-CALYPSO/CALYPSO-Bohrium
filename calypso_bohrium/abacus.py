#!/usr/bin/env python
import os
import sys
import glob
import shutil

from dpdispatcher import Task


def abacus_command(N_INCAR, ncpu):
    if N_INCAR == 1:
        command_runabacus = f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1"

    elif N_INCAR == 2:
        command_runabacus = f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1;"
        command_runabacus += " if [ -f OUT.ABACUS/STRU_ION_D ]; then cp OUT.ABACUS/STRU_ION_D ./STRU;cp INPUT_2 INPUT; "
        command_runabacus += f"mpirun -n {ncpu} abacus > fp.log 2>&1; "
        command_runabacus += "else : ; fi "

    elif N_INCAR == 3:
        command_runabacus = f"cp INPUT_1 INPUT; mpirun -n {ncpu} abacus > fp.log 2>&1;"
        command_runabacus += " if [ -f OUT.ABACUS/STRU_ION_D ]; then cp OUT.ABACUS/STRU_ION_D ./STRU;cp INPUT_2 INPUT; "
        command_runabacus += f"mpirun -n {ncpu} abacus > fp.log 2>&1; "
        command_runabacus += "else : ; fi; "
        command_runabacus += " if [ -f OUT.ABACUS/STRU_ION_D ]; then cp OUT.ABACUS/STRU_ION_D ./STRU;cp INPUT_3 INPUT; "
        command_runabacus += f"mpirun -n {ncpu} abacus > fp.log 2>&1; "
        command_runabacus += "else : ; fi "

    return command_runabacus


def abacus_task(
    pop,
    task_dir,
    N_INCAR,
    command,
    backward_files=["STRU_ION_D", "running_cell-relax.log"],
    *args,
    **kwargs,
):
    if not os.path.exists('pickup') or (
        os.path.exists('pickup') and os.path.exists('restart')
    ):
        shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR.ori"))
        shutil.copyfile("STRU_%d" % pop, os.path.join(task_dir, "STRU"))
        shutil.copyfile("STRU_%d" % pop, os.path.join(task_dir, "STRU.ori"))
        for n_incar in range(1, N_INCAR + 1):
            shutil.copyfile(
                "INPUT_%d" % n_incar, os.path.join(task_dir, "INPUT_%d" % n_incar)
            )

        file_list = glob.glob("*.upf")
        file_list.extend(glob.glob("*.orb"))
        for temp_file in file_list:
            shutil.copyfile(temp_file, os.path.join(task_dir, temp_file))

    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=["STRU", "STRU.ori"]
        + [f"INPUT_{idx}" for idx in range(1, N_INCAR + 1)]
        + file_list,
        backward_files=backward_files,
    )
    return task


def abacus_back(task_dir, pop, *args, **kwargs):
    back_stru_file = os.path.join(task_dir, "OUT.ABACUS", "STRU_ION_D")

    if os.path.exists(back_stru_file):
        shutil.copyfile(back_stru_file, "STRU_ION_D_%d" % pop)
    else:
        back_stru_file = os.path.join(task_dir, "STRU.ori")
        shutil.copyfile(back_stru_file, "STRU_ION_D_%d" % pop)

    shutil.copyfile(
        os.path.join(task_dir, "OUT.ABACUS", "running_cell-relax.log"),
        "running_cell-relax.log_%d" % pop,
    )
