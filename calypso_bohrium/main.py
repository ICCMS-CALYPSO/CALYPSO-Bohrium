#!/usr/bin/env python
import os
import sys
import shutil
from pathlib import Path

import click
from dpdispatcher import Machine, Resources, Submission

from calypso_bohrium.utils import get_value
from calypso_bohrium.vasp import vasp_command, vasp_task, vasp_back
from calypso_bohrium.qe import qe_command, qe_task, qe_back


dft_task = {'vasp':vasp_task, 'qe':qe_task}
task_back = {'vasp':vasp_back, 'qe':qe_back}
dft_command = {'vasp':vasp_command, 'qe':qe_command}

@click.command()
@click.option("--dft", default="vasp", type=click.Choice(['vasp', 'qe']), help="dft calculator selection, support vasp and qe currently")
def main(dft):

    if get_value("Split") == "":
        os.system("echo 'Split = T' >> ./input.dat") 

    machine = Machine.load_from_json("machine.json")
    resources = Resources.load_from_json("resources.json")
    
    _ncpu = machine.input_data["machine_type"].split("_")[0].strip("c")
    ncpu = _ncpu if _ncpu != "" else 32
    
    MaxStep = int(get_value("MaxStep"))
    PopSize = int(get_value("PopSize"))
    N_INCAR = int(get_value("NumberOfLocalOptim"))
    
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
            os.system("./calypso.x | tee caly.log")
        elif os.path.exists('pickup') and not os.path.exists('restart'):
            pass
        elif (os.path.exists('pickup') and os.path.exists('restart')):
            os.system("./calypso.x | tee caly.log")
            step += 1
    
        task_list = []
        for pop in range(1, PopSize + 1):
            task_dir = "./data/step%03d.pop%03d" % (step, pop)
            Path(task_dir).mkdir(parents=True, exist_ok=True)
    
            task = dft_task[dft](pop, task_dir, N_INCAR, command)
            task_list.append(task)

        submission = Submission(
            work_base=os.getcwd(), machine=machine, resources=resources, task_list=task_list
        )
        submission.run_submission()

        if os.path.exists('pickup'):
            os.system('mv pickup pickup_done')
        if os.path.exists('restart'):
            os.system('mv restart restart_done')
    
        for pop in range(1, PopSize + 1):
            task_dir = "./data/step%03d.pop%03d" % (step, pop)
            task_back[dft](task_dir, pop)

        os.system('mv *.sub lbg-*.sh *_fail *_finished log_dir/')
    
    os.system("./calypso.x | tee caly.log")

