#!/usr/bin/env python
import os
import shutil

from dpdispatcher import Task

qe = '''Nat=`grep 'number of k points' -B 2 out.pw | head -n 1 | awk {'print($1)'}`;StruLine=`expr $Nat + 5 `; grep 'CELL_' -A $StruLine out.pw | tail -n `expr $StruLine + 1` > pwscf'''


def get_pp(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        if 'pseudo_dir' in line:
            pp_dir = (
                line.split('=')[-1].strip().strip(',').strip().strip("'").strip(',')
            )
        if 'ntyp' in line:
            ntyp = int(line.split('=')[-1].strip().strip(',').strip())
            for idx, line in enumerate(lines):
                if 'ATOMIC_SPECIES' in line:
                    pp_name = []
                    for i in range(ntyp):
                        pp_name.append(
                            lines[idx + i + 1]
                            .strip()
                            .strip(',')
                            .strip("'")
                            .strip('"')
                            .split()[-1]
                            .strip()
                        )
    return (pp_dir, pp_name)


def get_pwscf_natoms(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return len(lines) - 5


def update_input(filename, natoms):
    with open(filename, 'r') as fin:
        lines = fin.readlines()
        for idx, line in enumerate(lines):
            if 'nat' in line:
                lines[idx] = f'nat = {str(natoms)},\n'
                break
    with open(filename, 'w') as fout:
        for line in lines:
            fout.write(line)


def qe_command(N_INCAR, ncpu):
    command_runvasp_list = ['cp pwscf.ori pwscf']
    command_runvasp_list += [
        f'cat pw_input_{idx} pwscf > pw_input; mpirun -n {ncpu} pw.x -pd .true. < pw_input > out.pw_{idx} 2>&1 ; cp out.pw_{idx} out.pw; {qe}'
        for idx in range(1, N_INCAR + 1)
    ]  # cpu number how to detect
    command_runvasp = ';'.join(command_runvasp_list)
    return command_runvasp


def qe_task(
    pop, task_dir, N_INCAR, command, backward_files=['out.pw', 'pwscf', 'log', 'err']
):
    pp_dir, _pp_name = get_pp('pw_input_1')
    if not os.path.exists('pickup') or (
        os.path.exists('pickup') and os.path.exists('restart')
    ):
        # depand on the type of dft software
        shutil.copyfile("pwscf_%d" % pop, os.path.join(task_dir, "pwscf"))
        shutil.copyfile("pwscf_%d" % pop, os.path.join(task_dir, "pwscf.ori"))
        natoms = get_pwscf_natoms("pwscf_%d" % pop)
        for n_incar in range(1, N_INCAR + 1):
            update_input("pw_input_%d" % n_incar, natoms)
            shutil.copyfile(
                "pw_input_%d" % n_incar, os.path.join(task_dir, "pw_input_%d" % n_incar)
            )
            for pp_name in _pp_name:
                shutil.copyfile(pp_dir + '/' + pp_name, os.path.join(task_dir, pp_name))

    task = Task(
        command=command,
        task_work_path=task_dir,
        forward_files=['pwscf', 'pwscf.ori']
        + [f'pw_input_{idx}' for idx in range(1, N_INCAR + 1)]
        + [p_name for p_name in _pp_name],
        # backward_files = ['out.pw', 'pwscf', 'log', 'err'] + [f'out.pw_{idx}' for idx in range(1, N_INCAR + 1)]
        backward_files=backward_files,
    )
    return task


def qe_back(task_dir, pop):
    # natoms = get_pwscf_natoms(os.path.join(task_dir, "pwscf.ori"))
    # os.system(f'echo {str(natoms)} > natoms ')
    # os.system('cat natoms {os.path.join(task_dir, "pwscf.ori")} > {os.path.join(task_dir, "pwscf.ori")}')
    shutil.copyfile(os.path.join(task_dir, "out.pw"), "out.pw_%d" % pop)
    shutil.copyfile(os.path.join(task_dir, "pwscf.ori"), "pwscf_%d" % pop)
