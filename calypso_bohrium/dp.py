#!/usr/bin/env python
import os
import shutil

from dpdispatcher import Task


def dp_command(N_INCAR, ncpu):
    pre_command = ""
    command_rundp_list = [
        "python calypso_run_opt.py >> opt.log 2>&1 || python calypso_check_outcar.py >> opt.log 2>&1"
    ]
    command_rundp = ";".join(command_rundp_list)

    command = pre_command + command_rundp

    return command


def dp_task(
    pop,
    task_dir,
    N_INCAR,
    command,
    backward_files=["CONTCAR", "OUTCAR", "opt.log", "traj.traj", "err"],
    lsurface=False,
):
    if not os.path.exists('pickup') or (
        os.path.exists('pickup') and os.path.exists('restart')
    ):
        shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR"))
        shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR.ori"))

    shutil.copyfile("input.dat", os.path.join(task_dir, "input.dat"))
    shutil.copyfile("frozen_model.pb", os.path.join(task_dir, "frozen_model.pb"))
    shutil.copyfile("calypso_run_opt.py", os.path.join(task_dir, "calypso_run_opt.py"))
    shutil.copyfile(
        "calypso_check_outcar.py", os.path.join(task_dir, "calypso_check_outcar.py")
    )

    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=[
            "POSCAR",
            "frozen_model.pb",
            "input.dat",
            "calypso_run_opt.py",
            "calypso_check_outcar.py",
        ],
        backward_files=backward_files,
    )
    return task


def dp_back(task_dir, pop, lsurface=False):
    shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d" % pop)
    shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d" % pop)
