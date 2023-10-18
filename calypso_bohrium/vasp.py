#!/usr/bin/env python
import os
import shutil

from dpdispatcher import Task


def vasp_command(N_INCAR, ncpu):
    pre_command = "cp POSCAR.ori CONTCAR;"
    command_runvasp_list = [
        f"cp CONTCAR POSCAR; cp INCAR_{idx} INCAR; mpirun -n {ncpu} vasp_std > fp.log 2>&1"
        for idx in range(1, N_INCAR + 1)
    ]  # cpu number how to detect
    command_runvasp = ";".join(command_runvasp_list)

    command = pre_command + command_runvasp

    return command


def vasp_task(
    pop,
    task_dir,
    N_INCAR,
    command,
    backward_files=["CONTCAR", "OUTCAR", "fp.log", "OSZICAR", "err"],
    lsurface="F",
):
    if not os.path.exists('pickup') or (
        os.path.exists('pickup') and os.path.exists('restart')
    ):
        if lsurface == "F":
            shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR"))
            shutil.copyfile("POSCAR_%d" % pop, os.path.join(task_dir, "POSCAR.ori"))
            for n_incar in range(1, N_INCAR + 1):
                shutil.copyfile(
                    "INCAR_%d" % n_incar, os.path.join(task_dir, "INCAR_%d" % n_incar)
                )
                shutil.copyfile("POTCAR", os.path.join(task_dir, "POTCAR"))
            forward_files = ["POSCAR", "POTCAR", "POSCAR.ori"] + [
                f"INCAR_{idx}" for idx in range(1, N_INCAR + 1)
            ]

        elif lsurface == "T":
            command = "python surface_run.py > log 2>&1"
            forward_files = [
                "INCAR-1",
                "INCAR-2",
                "KPOINTS-1",
                "KPOINTS-2",
                "POSCAR",
                "POTCAR",
                "SD-1",
                "SD-2",
                "submit.sh",
                "surface_run.py",
            ]
            backward_files = []

    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=forward_files,
        backward_files=backward_files,
    )
    return task


def vasp_back(task_dir, pop, lsurface="F"):
    if lsurface == "F":
        shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d" % pop)
        shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d" % pop)
    elif lsurface == "T":
        pass
