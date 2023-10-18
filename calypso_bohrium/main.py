#!/usr/bin/env python
import os
from pathlib import Path

import click
from dpdispatcher import Machine, Resources, Submission

from calypso_bohrium.utils import get_value

from calypso_bohrium.dp import dp_command, dp_task, dp_back
from calypso_bohrium.qe import qe_command, qe_task, qe_back
from calypso_bohrium.vasp import vasp_command, vasp_task, vasp_back
from calypso_bohrium.dftb import dftb_command, dftb_task, dftb_back
from calypso_bohrium.abacus import abacus_command, abacus_task, abacus_back


dft_task = {"vasp": vasp_task, "qe": qe_task, "abacus": abacus_task, "dp": dp_task, "dftb": dftb_task}
task_back = {"vasp": vasp_back, "qe": qe_back, "abacus": abacus_back, "dp": dp_back, "dftb": dftb_back}
dft_command = {
    "vasp": vasp_command,
    "qe": qe_command,
    "abacus": abacus_command,
    "dp": dp_command,
    "dftb": dftb_command,
}


@click.command()
@click.option(
    "--dft",
    default="vasp",
    type=click.Choice(["vasp", "qe", "abacus", "dp", "dftb"]),
    help="dft calculator selection, support vasp qe abacus and dp currently",
)
def main(dft):
    if get_value("Split") == "":
        os.system("echo 'Split = T' >> ./input.dat")

    lsurface = get_value("lsurface")  # "T"
    if lsurface == "T":
        lsurface = True
    else:
        lsurface = False

    machine = Machine.load_from_json("machine.json")
    resources = Resources.load_from_json("resources.json")

    out_files = machine.input_data.get('out_files', [])

    _ncpu = machine.input_data["machine_type"].split("_")[0].strip("c")
    ncpu = _ncpu if _ncpu != "" else 32

    MaxStep = int(get_value("MaxStep"))
    PopSize = int(get_value("PopSize"))
    number_of_local_optim = get_value("NumberOfLocalOptim")
    if number_of_local_optim is not "":
        N_INCAR = int(number_of_local_optim)
    else:
        if lsurface:
            N_INCAR = 2
        else:
            N_INCAR = 1

    if get_value("PickUp").lower().startswith("t"):
        PickStep = int(get_value("PickStep"))
        os.system(f"echo {PickStep+1} > step")
        os.system('touch pickup')
    else:
        if os.path.exists("step"):
            PickStep = int(os.popen("cat step").read())
        else:
            PickStep = 1

    command_intel = "source /opt/intel/oneapi/setvars.sh;"
    command_rundft = dft_command[dft](N_INCAR, ncpu)
    command = command_intel + command_rundft

    Path('./log_dir').mkdir(parents=True, exist_ok=True)

    for step in range(PickStep, MaxStep + 1):
        if not os.path.exists('pickup'):
            os.system("calypso.x | tee caly.log")
        elif os.path.exists('pickup') and not os.path.exists('restart'):
            pass
        elif os.path.exists('pickup') and os.path.exists('restart'):
            os.system("calypso.x | tee caly.log")
            step += 1

        task_list = []
        for pop in range(1, PopSize + 1):
            if lsurface :
                task_dir = "./results/Generation_%d/Indv_%d" % (step, pop)
            else:
                task_dir = "./data/step%03d.pop%03d" % (step, pop)
                Path(task_dir).mkdir(parents=True, exist_ok=True)

            task = dft_task[dft](
                pop, task_dir, N_INCAR, command, out_files, lsurface=lsurface
            )
            task_list.append(task)

        submission = Submission(
            work_base=os.getcwd(),
            machine=machine,
            resources=resources,
            task_list=task_list,
        )
        submission.run_submission()

        if os.path.exists('pickup'):
            os.system('mv pickup pickup_done')
        if os.path.exists('restart'):
            os.system('mv restart restart_done')

        for pop in range(1, PopSize + 1):
            if lsurface :
                task_dir = "./results/Generation_%d/Indv_%d" % (step, pop)
            else:
                task_dir = "./data/step%03d.pop%03d" % (step, pop)
            task_back[dft](task_dir, pop, lsurface=lsurface)

    os.system("calypso.x | tee caly.log")
