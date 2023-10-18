#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np


def get_element_num(elements):
    '''Using the Atoms.symples to Know Element&Num'''
    element = []
    ele = {}
    element.append(elements[0])
    for x in elements:
        if x not in element:
            element.append(x)
    for x in element:
        ele[x] = elements.count(x)
    return element, ele


def write_contcar(element, ele, lat, pos, task_dir):
    '''Write CONTCAR'''
    f = open(os.path.join(task_dir, 'CONTCAR'), 'w')
    f.write('CALYPSO\n')
    f.write('1.0\n')
    for i in range(3):
        f.write('%15.10f %15.10f %15.10f\n' % tuple(lat[i]))
    for x in element:
        f.write(x + '  ')
    f.write('\n')
    for x in element:
        f.write(str(ele[x]) + '  ')
    f.write('\n')
    f.write('Direct\n')
    na = sum(ele.values())
    dpos = np.dot(pos, np.linalg.inv(lat))
    for i in range(na):
        f.write('%15.10f %15.10f %15.10f\n' % tuple(dpos[i]))


def write_outcar(element, ele, volume, lat, pos, ene, force, stress, pstress, task_dir):
    '''Write OUTCAR'''
    # f = open('OUTCAR','w')
    f = open(os.path.join(task_dir, 'OUTCAR'), 'w')
    for x in element:
        f.write('VRHFIN =' + str(x) + '\n')
    f.write('ions per type =')
    for x in element:
        f.write('%5d' % ele[x])
    f.write(
        '\nDirection     XX             YY             ZZ             XY             YZ             ZX\n'
    )
    f.write('in kB')
    f.write('%15.6f' % stress[0])
    f.write('%15.6f' % stress[1])
    f.write('%15.6f' % stress[2])
    f.write('%15.6f' % stress[3])
    f.write('%15.6f' % stress[4])
    f.write('%15.6f' % stress[5])
    f.write('\n')
    ext_pressure = np.sum(stress[0] + stress[1] + stress[2]) / 3.0 - pstress
    f.write(
        'external pressure = %20.6f kB    Pullay stress = %20.6f  kB\n'
        % (ext_pressure, pstress)
    )
    f.write('volume of cell : %20.6f\n' % volume)
    f.write('direct lattice vectors\n')
    for i in range(3):
        f.write('%10.6f %10.6f %10.6f\n' % tuple(lat[i]))
    f.write('POSITION                                       TOTAL-FORCE(eV/Angst)\n')
    f.write('-------------------------------------------------------------------\n')
    na = sum(ele.values())
    for i in range(na):
        f.write('%15.6f %15.6f %15.6f' % tuple(pos[i]))
        f.write('%15.6f %15.6f %15.6f\n' % tuple(force[i]))
    f.write('-------------------------------------------------------------------\n')
    f.write('energy  without entropy= %20.6f %20.6f\n' % (ene, ene / na))
    enthalpy = ene + pstress * volume / 1602.17733
    f.write('enthalpy is  TOTEN    = %20.6f %20.6f\n' % (enthalpy, enthalpy / na))


def write_files(atoms, pstress, is_success, task_dir):
    # pstress kbar
    # kBar to eV/A^3
    # 1 eV/A^3 = 160.21766028 GPa
    # 1 / 160.21766028 ~ 0.006242
    # print(pstress, is_success)

    atoms_lat = atoms.cell
    atoms_pos = atoms.positions
    atoms_force = (
        atoms.get_forces() if is_success else np.zeros((atoms_pos.shape[0], 3))
    )
    atoms_stress = atoms.get_stress() if is_success else np.array([0, 0, 0, 0, 0, 0])
    # eV/A^3 to GPa
    atoms_stress = atoms_stress / (0.01 * 0.6242)
    atoms_symbols = atoms.get_chemical_symbols()
    atoms_ene = atoms.get_potential_energy() if is_success else 610612509
    atoms_vol = atoms.get_volume()
    element, ele = get_element_num(atoms_symbols)

    write_contcar(element, ele, atoms_lat, atoms_pos, task_dir)
    write_outcar(
        element,
        ele,
        atoms_vol,
        atoms_lat,
        atoms_pos,
        atoms_ene,
        atoms_force,
        atoms_stress * -10.0,
        pstress,
        task_dir,
    )
