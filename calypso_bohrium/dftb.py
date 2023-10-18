#!/usr/bin/env python
import os
import glob
import shutil

from dpdispatcher import Task


def dftb_command(N_INCAR, ncpu):
    pre_command = "cp dftb.gen geom.out.gen; cp dftb.gen dftb.gen.ori;"
    command_rundftb_list = [
        f"cp geom.out.gen dftb.gen; cp dftb_in.hsd-{idx} dftb_in.hsd; mpirun -n {ncpu} dftb+ > fp.log 2>&1;"
        for idx in range(1, N_INCAR + 1)
    ]
    command_rundftb = ";".join(command_rundftb_list)

    command = pre_command + command_rundftb

    return command


def dftb_task(
    pop,
    task_dir,
    N_INCAR,
    command,
    backward_files=[
        "geom.out.gen",
        "geom.out.xyz",
        "fp.log",
        "detail.out",
        "dftb_pin.hsd",
        "err",
    ],
    lsurface="F",
):
    if not os.path.exists('pickup') or (
        os.path.exists('pickup') and os.path.exists('restart')
    ):
        if lsurface == "F":
            raise NotImplementedError("DFTB is only supported in surface!!!")
            # shutil.copyfile("dftb.gen-%d" % pop, os.path.join(task_dir, "dftb.gen"))
            # shutil.copyfile("dftb.gen-%d" % pop, os.path.join(task_dir, "dftb.gen.ori"))
            # for n_incar in range(1, N_INCAR + 1):
            #     shutil.copyfile(
            #         "dftb_in.hsd_%d" % n_incar, os.path.join(task_dir, "dftb_in.hsd_%d" % n_incar)
            #     )
            #     shutil.copyfile("POTCAR", os.path.join(task_dir, "POTCAR"))
            # forward_files = ["POSCAR", "POTCAR", "POSCAR.ori"] + [
            #     f"INCAR_{idx}" for idx in range(1, N_INCAR + 1)
            # ]

        elif lsurface == "T":
            command = "python surface_run.py > log 2>&1"
            _forward_files = [
                "dftb_in.hsd-1",
                "dftb_in.hsd-2",
                "dftb.gen",
                "submit.sh",
                "surface_run.py",
            ]
            skf_list = glob.glob(f"{task_dir}/*.skf")
            forward_files = _forward_files + skf_list
            backward_files = []

    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=forward_files,
        backward_files=backward_files,
    )
    return task


def dftb_back(task_dir, pop, lsurface="F"):
    if lsurface == "F":
        pass
        # shutil.copyfile(os.path.join(task_dir, "CONTCAR"), "CONTCAR_%d" % pop)
        # shutil.copyfile(os.path.join(task_dir, "OUTCAR"), "OUTCAR_%d" % pop)
    elif lsurface == "T":
        pass
