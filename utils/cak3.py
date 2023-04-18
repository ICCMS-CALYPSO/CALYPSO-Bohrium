#!/opt/mamba/bin/python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of CALYPSO PACKAGE.
# Spglib is used for finding and handling crystal symmetries.
#
# Author : Li Zhu < zl@calypso.cn >
# Version 2014.05.15
#
# Author : Xuecheng Shao
# Version 2016.10.10 add XRD and TS
# Version 2018.03.30 only for vsc (three elements)
# Version 2018.06.30 add mine and maxe
# Version 2018.07.02 remove same spacegroup and small energy difference structures
# Version 2018.10.19 add split vsc struct.dat --split
#
# Author : Zhenyu Wang <wzy@calypso.cn>
# Version 2022.07.01 update to python3 version with spglib-v1.16.0
# Version 2022.10.12 update to python3 version with spglib-v2.0.1
#
# Copyright CALYPSO Develop Group

import os
import math
import sys
import time
import glob

from optparse import OptionParser
try:
    import spglib as spg
except:
    print('ERROR: please install spglib.')
    exit(0)
try:
    import numpy as np
except:
    print('ERROR: please install numpy.')
    exit(0)


def readinput():
    if os.path.exists('input.dat'):
        finput = 'input.dat'
    elif os.path.exists('../input.dat'):
        finput = '../input.dat'
    else:
        print('Error: no input file: input.dat')
        exit(0)

    f = open(finput, 'r')
    indata = {}
    for line in f:
        if '=' in line:
            if line.strip()[0] != '#':
                litem = line.split('=')
                indata[litem[0].strip().lower()] = litem[1].strip()
    f.close()

    # NameOfAtoms = indata['nameofatoms'].split()
    Npop = int(indata['popsize'])
    try:
        SystemName = indata['systemname']
    except:
        SystemName = 'CALYPSO'
    # print indata['Mol'][0]

    try:
        Mtmp = indata['mol'][0]
        if Mtmp.upper() == 'T':
            Mol = True
        else:
            Mol = False
    except:
        # print 'exce'
        Mol = False

    try:
        H = indata['hardness'][0]
        if H.upper() == 'T':
            Hard = True
        else:
            Hard = False
    except:
        Hard = False

    try:
        vsc = indata['vsc'][0]
        if vsc.upper() == 'T':
            VSC = True
        else:
            VSC = False
    except:
        VSC = False


    try:
        d2 = indata['2d'][0]
        if d2.upper() == 'T':
            D2 = True
        else:
            D2 = False
    except:
        D2 = False

    try:
        cl = indata['cluster'][0]
        if cl.upper() == 'T':
            CL = True
        else:
            CL = False
    except:
        CL = False
    try:
        lsur = indata['lsurface'][0]
        if lsur.upper() == 'T':
            LSUR = True
        else:
            LSUR = False
    except:
        LSUR = False

    if LSUR:
        NameOfAtoms = []
    else:
        NameOfAtoms = indata['nameofatoms'].split()

    try:
        vsce = indata['vscenergy'].split()
        if len(vsce) < len(NameOfAtoms):
            print('VSCEnergy Error.')
            exit(0)
        else:
            VSCE = list(map(float, vsce))
    except:
        VSCE = [0]

    try:
        hm = indata['halfmetal'][0]
        if hm.upper() == 'T':
            HM = True
        else:
            HM = False
    except:
        HM = False

    try:
        bg = indata['bandgapdesign'][0]
        if bg.upper() == 'T':
            BG = True
        else:
            BG = False
    except:
        BG = False
    try:
        xrd = indata['lxrd'][0]
        if xrd.upper() == 'T':
            XRD = True
        else:
            XRD = False
    except:
        XRD = False
    try:
        ts = indata['ltranstate'][0]
        if ts.upper() == 'T':
            TS = True
            num_neb = int(indata['numberofimages'])
        else:
            TS = False
            num_neb=[0]
    except:
        TS = False
        num_neb=[0]
    # print num_neb
    # aaa=raw_input()



    return (SystemName, NameOfAtoms, Npop, Mol, Hard, VSC, VSCE, D2, CL, HM, LSUR, BG, XRD, TS, num_neb)

def parseStruct_old():
    """docstring for parseStruct"""
    try:
        f = open('struct.dat')
    except:
        print("Error: no struct.dat")
        exit(0)
    date = []
    try:
        for line in f:
            date.append(line)
    except:
        f.close()

    ldate = len(date)
    structData = []
    for i in range(ldate):
        if 'nstruct' in date[i]:
            hardness = 0.
            polar = 0.
            bandgap = 100.
            xrddiff=610612509.
            nstruct = int(date[i].split('=')[1].strip())
            # print 'ns', nstruct
            gtype = 1
            for j in range(i+1, ldate):
                if 'GenType' in date[j]:
                    gtype = date[j].split(':')[1].strip()
                if 'Optimized' in date[j]:
                    for k in range(j+1, ldate):
                        if 'Energy' in date[k]:
                            try:
                                energy = float(date[k].split('=')[1])
                            except:
                                energy = 610612509.
                        if 'Hardness' in date[k]:
                            try:
                                hardness = float(date[k].split('=')[1])
                            except:
                                hardness = 0.
                        if 'Polar' in date[k]:
                            try:
                                polar = abs(float(date[k].split('=')[1]))
                                if math.isnan(polar): polar = 0.
                            except:
                                polar = 0.
                        #if 'diff_gap' in date[k]:
                        if 'DiffGap' in date[k]:
                            try:
                                bandgap = float(date[k].split('=')[1])
                                if math.isnan(bandgap): bandgap = 100.
                            except:
                                bandgap = 100.
                        if 'XRDdiff' in date[k]:
                            try:
                                xrddiff = float(date[k].split('=')[1])
                                if math.isnan(xrddiff): xrddiff = 610612509.
                            except:
                                xrddiff = 610612509.
                        if 'Number Species' in date[k]:
                            ntype = int(date[k].split('=')[1])
                        if 'Ele_Num' in date[k]:
                            typt = list(map(int, date[k].split('=')[1].split()))
                            natom = sum(typt)
                            #print typt,natom
                        if 'lat_matrix' in date[k]:
                            # print  date[k+1 : k+4]
                            lat = np.array([item.split() for item in date[k+1 : k+4]], float)
                        if 'Atomic' in date[k]:
                            pos = np.array([item.split() for item in date[k+1 : k+1+natom]], float)

                        if 'nstruct' in date[k]:
                            break
                    break

            structData.append( [energy, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap, xrddiff] )
            # print '#', nstruct, gtype, energy, typt
            # print lat
            # print pos
    return structData

def parseStruct():
    """docstring for parseStruct"""
    try:
        f = open('struct.dat')
    except:
        print("Error: no struct.dat")
        exit(0)
    date = []
    try:
        for line in f:
            date.append(line)
    except:
        f.close()

    ldate = len(date)
    structData = []
    i = 0
    ldate = len(date)
    while True:
        if 'nstruct' in date[i]:
            hardness = 0.
            polar = 0.
            bandgap = 100.
            xrddiff=610612509.
            nstruct = int(date[i].split('=')[1].strip())
            # print 'ns', nstruct
            gtype = 1
            # for j in range(i+1, ldate):
            j = i
            lcurr = True
            while lcurr :
                j = j + 1
                if 'GenType' in date[j]:
                    gtype = date[j].split(':')[1].strip()
                if 'Optimized' in date[j]:
                    # for k in range(j+1, ldate):
                    k = j
                    while lcurr :
                        k = k + 1
                        # print k, date[k]
                        line = date[k]
                        if 'Energy' in line:
                            try:
                                energy = float(line.split('=')[1])
                            except:
                                energy = 610612509.
                        elif 'Hardness' in line:
                            try:
                                hardness = float(line.split('=')[1])
                            except:
                                hardness = 0.
                        elif 'Polar' in line:
                            try:
                                polar = abs(float(line.split('=')[1]))
                                if math.isnan(polar): polar = 0.
                            except:
                                polar = 0.
                        #if 'diff_gap' in line:
                        elif 'DiffGap' in line:
                            try:
                                bandgap = float(line.split('=')[1])
                                if math.isnan(bandgap): bandgap = 100.
                            except:
                                bandgap = 100.
                        elif 'XRDdiff' in line:
                            try:
                                xrddiff = float(line.split('=')[1])
                                if math.isnan(xrddiff): xrddiff = 610612509.
                            except:
                                xrddiff = 610612509.
                        elif 'Number Species' in line:
                            ntype = int(line.split('=')[1])
                        elif 'Ele_Num' in line:
                            typt = list(map(int, line.split('=')[1].split()))
                            natom = sum(typt)
                        elif 'lat_matrix' in line:
                            # print  date[k+1 : k+4]
                            lat = np.array([item.split() for item in date[k+1 : k+4]], float)
                            k = k + 3
                        elif 'Atomic' in line:
                            pos = np.array([item.split() for item in date[k+1 : k+1+natom]], float)
                            k = k + natom
                        elif 'nstruct' in line or k +3 > ldate :
                            i = k
                            lcurr = False

            structData.append( [energy, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap, xrddiff] )
            # print 'i', i, ldate
            if i +10 > ldate - 1  :
                break
            # print '#', nstruct, gtype, energy, typt
            # print lat
            # print pos
    return structData

def SplitVSCStruct():
    try:
        f = open('struct.dat')
    except:
        print("Error: no struct.dat")
        exit(0)
    date = []
    try:
        for line in f:
            date.append(line)
    except:
        f.close()

    ldate = len(date)
    structData = []
    AllData= {}
    i = 0
    ldate = len(date)
    while True:
        fstr = ''
        if 'nstruct' in date[i]:
            nstruct = int(date[i].split('=')[1].strip())
            j = i + 1
            lcurr = True
            while lcurr :
                fstr  += date[j]
                j = j + 1
                if 'Ele_Num' in date[j] :
                    typt = list(map(int, date[j].split('=')[1].split()))
                    typt = gcdnum(typt)
                    fus = '_'.join(map(str, typt))
                elif 'nstruct' in date[j] or j +3 > ldate :
                    if fus in AllData :
                        fstr = '            nstruct= ' + str(len(AllData[fus]) + 1)  + '\n'+ fstr
                        AllData[fus].append(fstr)
                    else :
                        fstr = '            nstruct= 1' + '\n' + fstr
                        AllData[fus] = [fstr]
                    i = j
                    lcurr = False

            if i +10 > ldate - 1  :
                break

    if fus in AllData :
        fstr = '            nstruct= ' + str(len(AllData[fus]) + 1)  + '\n'+ fstr
        AllData[fus].append(fstr)
    else :
        fstr = '            nstruct= 1' + '\n' + fstr
        AllData[fus] = [fstr]
    for fus in AllData :
        # print fus
        # print AllData[fus][2]
        fdir = 'data_' + fus + '/'
        os.system('mkdir -p ' + fdir)
        outf = fdir + 'struct.dat'
        with open(outf, 'w')  as fw :
            for line in AllData[fus] :
                fw.write(line)

    return 

def parseTS():
    """docstring for parseStruct"""
    try:
        f = open('TS.dat')
    except:
        print("Error: no TS.dat")
        exit(0)
    date = []
    try:
        for line in f:
            date.append(line)
    except:
        f.close()

    ldate = len(date)
    structData = []
    for i in range(ldate):
        if 'NTranState' in date[i]:
            nstruct = int(date[i].split('=')[1].strip())
        if 'Potential' in date[i]:
            pts=float(date[i].split('=')[1])
            for j in range(i+1, ldate):
                if 'Number Species' in date[j]:
                    ntype = int(date[j].split('=')[1])
                if 'Ele_num' in date[j]:
                    typt = list(map(int, date[j].split('=')[1].split()))
                    natom = sum(typt)
                    structure=[]
                    for k in range(j+1, ldate):
                        if 'structure' in date[k]:
                           # print date[k].split()[1].strip()
                           energy = float(date[k].split('=')[1].strip())
                        if 'lat_matrix' in date[k]:
                            lat = np.array([item.split() for item in date[k+1 : k+4]], float)
                        if 'Atomic' in date[k]:
                            pos = np.array([item.split() for item in date[k+1 : k+1+natom]], float)
                            structure.append((lat,pos,typt,natom,energy))
                        if 'NTranState' in date[k]:
                            break
                    break
            # print 'nstruct',nstruct
            # raw_input()
            structData.append( [pts, nstruct, structure] )
    return structData

def parseIni():
    """docstring for parseStruct"""
    try:
        f = open('struct.dat')
    except:
        print("Error: no struct.dat")
        exit(0)
    date = []
    try:
        for line in f:
            date.append(line)
    except:
        f.close()

    ldate = len(date)
    structData = []
    for i in range(ldate):
        if 'nstruct' in date[i]:
            hardness = 0.
            polar = 0.
            bandgap = 100.
            energy = 0.
            gtype = 'A'
            nstruct = int(date[i].split('=')[1].strip())
            for j in range(i+1, ldate):

                if 'Initial Structure' in date[j]:
                    for k in range(j+1, ldate):

                        if 'Number Species' in date[k]:
                            ntype = int(date[k].split('=')[1])
                        if 'Ele_Num' in date[k]:
                            typt = list(map(int, date[k].split('=')[1].split()))
                            natom = sum(typt)
                        if 'lat_matrix' in date[k]:
                            # print  date[k+1 : k+4]
                            lat = np.array([item.split() for item in date[k+1 : k+4]], float)
                        if 'Atomic' in date[k]:
                            pos = np.array([item.split() for item in date[k+1 : k+1+natom]], float)

                        if 'nstruct' in date[k]:
                            break
                    break

            structData.append( [energy, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap] )
            # print '#', nstruct, gtype, energy, typt
            # print lat
            # print pos
    return structData

def parseOpt(optf='opt', hm=False):
    if optf == 'opt':
        optfiles = glob.glob('pso_opt_*')
    else:
        optfiles = glob.glob('pso_dft_*')

    ofls = sorted(optfiles, key=lambda x: int(x.split('_')[2]))

    #print ofls
    hardness = 0.
    polar = 0.
    bandgap = 100.
    xrddiff=610612509.
    gtype = '(^3^)'
    nstruct = 0
    enth = []
    stru = []
    for optfile in ofls:
        fline = []
        f = open(optfile)
        try:
            for line in f:
                fline.append(line)
        finally:
            f.close()
        if hm:
            system = fline[2]
        else:
            system = fline[1]
        nstr = 0
        for item in fline:
            if item.strip() == system.strip(): nstr += 1
        num = 0
        for i in range(0, nstr):
            nstruct += 1
            if hm:
                polar = abs(float(fline[num].split()[0]))
                enthalpy = float(fline[num].split()[1])
                num += 1
            else:
                enthalpy = float(fline[num].split()[0])
                enth.append(enthalpy)
            natom = sum(map(int, fline[num+6].split()))
            #
            strutmp = fline[num+1 : num+natom+8]
            lat = []
            for i in range(2, 5):
                lat.append(list(map(float, strutmp[i].split())))
            lat = np.array(lat, float)
            typt = list(map(int, strutmp[5].split()))
            pos = []
            for item in strutmp[7:]:
                pos.append(list(map(float, item.split()[:3])))
            pos = np.array(pos, float)
            stru.append( [enthalpy, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap, xrddiff] )
            num = num + natom + 8

    return stru

def parseMol():

    psofiless = glob.glob('pso_struct_*')
    psofiles = sorted(psofiless, key=lambda x: int(x.split('_')[2]))
    d = []
    for psofile in psofiles:
        f = open(psofile)
        try:
            for line in f:
                d.append(line.split())
        finally:
            f.close()


    hardness = 0.
    bandgap = 100.
    xrddiff=610612509.
    polar = 0.
    gtype = '(^3^)'

    structData = []
    i = 0
    nstruct = 0
    while i < len(d):
        nstruct += 1
        e = float(d[i][0])
        lat = np.array(d[i+1:i+4], float)
        # print d[i+4]
        natom = int(d[i+4][0])
        # print natom
        atom = d[i+5:i+5+natom]

        atoms = sorted(atom, key=lambda x: int(x[0]))
        typ = [int(x[0]) for x in atoms]

        typt = []
        for x in set(typ):
            k = 0
            for item in typ:
                if item == x: k += 1
            typt.append(k)

        pos = np.array([x[2:] for x in atoms], float)

        for j in range(len(pos)):
            for k in range(3):
                pos[j][k] = pos[j][k] - np.floor(pos[j][k])

        structData.append( [e, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap,xrddiff ])
        i += (5+natom)
        # print 'i2',i

    # print len(structData)
    return structData

def parseSur():
    global name_ele
    try:
        f = open('CALYPSO.log')
    except:
        print('Error: no CALYPSO.log')
        sys.exit(0)

    hardness = 0.
    nstruct = 0
    polar = 0.
    bandgap = 100.
    xrddiff=610612509.
    data = []
    for line in f:
        data.append(line.split())
    f.close()

    structData = []

    for i in range(len(data)):
        if '/' not in data[i][0] and 'Generation' in data[i][0] and 'opt' in data[i][0]:
            for j in range(i+1, len(data)):
                if data[j][0] == 'Group_index':
                    for k in range(j+1, len(data)):
                        if data[k][0] == 'Gen':
                            gtype = data[k][1]
                        if data[k][0] == 'Fitness':
                            try:
                                energy = float(data[k][1])
                            except:
                                energy = 610612509.
                        if data[k][0] == 'Lattice_Matrix':
                            lat = np.array(data[k+1:k+4], float)
                        if data[k][0] == 'Atoms_Info':
                            natom = sum(map(int, data[k][1:]))
                            tpos = data[k+1:k+1+natom]
#                            sym = [ x[0] for x in tpos]
#                            pos = np.array([x[1:4] for x in tpos], float)
                            atominfo={}
                            for x in tpos:
                                atominfo.setdefault(x[0],[]).append([float(x[1]),float(x[2]),float(x[3]),x[7]])
                            NameOfAtoms=list(atominfo.keys())
                            typt=[len(atominfo[x]) for x in NameOfAtoms]
                            pos=[]
                            for x in NameOfAtoms:
                                for xx in range(len(atominfo[x])):
                                    pos.append(atominfo[x][xx])
#                            typt = []
#                            NameOfAtoms = list(set(sym))
#                            NameOfAtoms.sort(key = sym.index)
#                            NameOfAtoms=atominfo.keys()
#                            typt=[len(atominfo[x]) for x in NameOfAtoms]
#                            for x in NameOfAtoms:
#                                l = 0
#                                for item in sym:
#                                    if item == x: l += 1
#                                typt.append(l)
                        if data[k][0] == 'Group_index':
                            break
                    nstruct += 1
                    structData.append( [energy, hardness, nstruct, gtype, (lat, pos, typt, natom), polar, bandgap, xrddiff] )
                if '/' in data[j][0] and 'Generation' in data[j][0] and 'opt' in data[j][0]:
                    break
    name_ele = NameOfAtoms[:]
    return structData

def write_cif(pdir, cell, no, sgdata, wori, wprim):
    """convert the vasp poscar to cif file"""

    global name_ele
    tt = time.localtime()
    ttt = str(tt[0])+'-'+str(tt[1])+'-'+str(tt[2])+'\n'
    l = cell[0].copy()
    p = cell[1].copy()
    typt = cell[2][:]
    sg = sgdata[0]
    sgsymbol = sgdata[1]

    ra = math.sqrt(l[0][0]**2 + l[0][1]**2 + l[0][2]**2)
    rb = math.sqrt(l[1][0]**2 + l[1][1]**2 + l[1][2]**2)
    rc = math.sqrt(l[2][0]**2 + l[2][1]**2 + l[2][2]**2)

    cosinea = (l[1][0]*l[2][0] + l[1][1]*l[2][1] + l[1][2]*l[2][2])/rb/rc
    cosineb = (l[0][0]*l[2][0] + l[0][1]*l[2][1] + l[0][2]*l[2][2])/ra/rc
    cosinec = (l[0][0]*l[1][0] + l[0][1]*l[1][1] + l[0][2]*l[1][2])/ra/rb

    anglea = math.degrees(math.acos(cosinea))
    angleb = math.degrees(math.acos(cosineb))
    anglec = math.degrees(math.acos(cosinec))

    alpha = anglea
    beta = angleb
    gamma = anglec
    if wori:
        ofile = pdir + '/' + str(no) + ".cif"
        sgsymbol = 'P1'
    elif wprim:
        ofile = pdir + '/' + str(no) + '_' + str(sg) + "_p.cif"
    else:
        ofile = pdir +'/' + str(no) + '_' +str(sg)+".cif"
    f = open(ofile, 'w')
    f.write("data_" + sgsymbol + "\n"
            + "_audit_creation_date               " + ttt
            + "_audit_creation_method             'CALYPSO -> cif'\n"
            + "_symmetry_space_group_name_H-M     'P1'\n"
            + "_symmetry_Int_Tables_number        1\n"
            + "_symmetry_cell_setting             triclinic\n"
            + "loop_\n"
            + "_symmetry_equiv_pos_as_xyz\n"
            + "x,y,z\n"
            + "_cell_length_a" + ' '*8 + ("%8.4f" % ra)  + '\n'
            + "_cell_length_b" + ' '*8 + ("%8.4f" % rb) + '\n'
            + "_cell_length_c" + ' '*8 + ("%8.4f" % rc) + '\n'
            + "_cell_angle_alpha" + ' '*8 + ("%8.4f" % alpha) + '\n'
            + "_cell_angle_beta " + ' '*8 + ("%8.4f" % beta) + '\n'
            + "_cell_angle_gamma" + ' '*8 + ("%8.4f" % gamma) + '\n'
            + "loop_\n"
            + "_atom_site_label\n"
            + "_atom_site_type_symbol\n"
            + "_atom_site_fract_x\n"
            + "_atom_site_fract_y\n"
            + "_atom_site_fract_z\n"
            + "_atom_site_U_iso_or_equiv\n"
            + "_atom_site_adp_type\n"
            + "_atom_site_occupancy\n")
    k = 0
    num_ele = len(typt)
    for i in range(0, num_ele):
        for j in range(0, typt[i]):
            hao = str(j+1)
            f.write(name_ele[i].strip() + hao + " "*4 + name_ele[i] + ("%12.5f" % p[k][0]) +
                    ("%12.5f" % p[k][1]) + ("%12.5f" % p[k][2]) + "   0.00  Uiso   1.00\n")
            k += 1
    f.close()

def write_xyz(pdir, cell, no, sgdata, wori, wprim):

    global name_ele

    l = cell[0].copy()
    p = cell[1].copy()
    typt = cell[2][:]
    natom = sum(typt)
    sg = sgdata[0]
    sgsymbol = sgdata[1]
    ofile = pdir + '/' + str(no) + '.xyz'
    f = open(ofile, 'w')
    f.write(str(natom) + "\n")
    f.write("xyz file : generated by cak.py\n")

    k = 0

    for i in range(0, len(typt)):
        for j in range(0, typt[i]):
            x = l[0][0] * p[k][0]
            y = l[1][1] * p[k][1]
            z = l[2][2] * p[k][2]
            f.write(str(name_ele[i]) +
                    ("%12.7f " % x ) + ("%12.7f " % y ) + ("%12.7f " % z )
                    + "\n")
            k += 1

def write_vasp(pdir, cell, no, sgdata, wori, wprim):

    global name_ele

    l = cell[0].copy()
    p = cell[1].copy()
    typt = cell[2][:]
    sg = sgdata[0]
    sgsymbol = sgdata[1]

    if wori:
        ofile = pdir + '/OCell_' + str(no) + '.vasp'
        sgsymbol = 'P1'
    elif wprim:
        ofile = pdir + '/PCell_' + str(no) + '_' + str(sg) + '.vasp'
    else:
        ofile = pdir + '/UCell_' + str(no) + '_' + str(sg) + '.vasp'

    f = open(ofile, 'w')
    f.write(str(no) + ':' + sgsymbol + '(' + str(sg) + ')')
    if wprim:
        f.write(' PRIMITIVE\n')
    else:
        f.write('\n')
    f.write('1.0\n')
    for item in l:
        f.write("%12.7f %12.7f %12.7f\n" % tuple(item))
    for item in name_ele:
        f.write("%2s " % item)
    f.write('\n')
    for item in typt:
        f.write("%3d " % item)
    f.write('\n')
    f.write('Direct\n')
    for item in p:
        f.write("%8.5f %8.5f %8.5f\n" % tuple(item))
    f.close()

def findsym(cell, prec, is_refine, is_prim):
    '''
    cell[ lattice,
          positions,
          typt,
          num_atoms
    ]
    '''

    aprec = -1

    # sl = cell[0].T.copy()
    sl = cell[0].copy()
    sp = cell[1].copy()
    stypt = cell[2][:]
    num_atoms = cell[3]

    snumbers = []
    for i in range(len(stypt)):
        snumbers += [i+1] * stypt[i]
    snumbers = np.array(snumbers, int)

    l = sl[:]
    p = sp[:]
    typt = stypt[:]
    numbers = snumbers[:]

    spacegroup_info = spg.get_spacegroup((l, p, numbers), prec, aprec)
    if spacegroup_info is None:
        num_spg = 0
        symbol_spg = 'NULL'
    else:
        num_spg = spacegroup_info.split()[1].strip('(').strip(')')
        symbol_spg = spacegroup_info.split()[0]

    if is_refine:
        # pp = np.zeros((num_atoms*4, 3), float)
        # numn = np.array([0]*(num_atoms*4), int)

        # for i in range(0, num_atoms):
        #     pp[i] = p[i]
        #     numn[i] = numbers[i]
        if num_spg is not None:
            refine_info = spg.refine_cell((l, p, numbers), prec, aprec)
            if refine_info is not None:
                ll, pp, numn = refine_info 
            else:
                ll, pp, numn = (l, p, numbers)
        else:
            ll, pp, numn = (l, p, numbers)

        num_atoms_brv = len(numn)

        refine_lat = ll
        refine_pos = []
        refine_numbers = []
        refine_typt = []


        for inums in set(numbers):
            itypt = 0
            for i in range(num_atoms_brv):
                if inums == numn[i]:
                    refine_pos.append(pp[i])
                    refine_numbers.append(numn[i])
                    itypt += 1
            refine_typt.append(itypt)

        refine_lat = np.array(refine_lat, float)
        refine_pos = np.array(refine_pos, float)
        refine_numbers = np.array(refine_numbers, int)
    else:
        refine_lat = None
        refine_pos = None
        refine_typt = None
        num_atoms_brv = None
        refine_numbers = None

    if is_prim and refine_lat is not None:
        _p_l = refine_lat.copy()
        _p_p = refine_pos.copy()
        _p_typt = refine_typt[:]
        _p_numbers = refine_numbers[:]
        if num_spg is not None:
            prim_info = spg.find_primitive((_p_l, _p_p, _p_numbers), 1e-5, -1)
            if prim_info is not None:
                p_l, p_p, p_numbers = prim_info
            else:
                p_l, p_p, p_numbers = (_p_l, _p_p, _p_numbers)
        else:
            p_l, p_p, p_numbers = (_p_l, _p_p, _p_numbers)
        num_atom_prim = len(p_numbers)
        if num_atom_prim > 0:
            prim_lat = p_l.copy()
            prim_positions = p_p[:]
            prim_numbers = p_numbers[:]
            # prim_positions = []
            # prim_numbers = []
            prim_typt = []
            for inums in set(prim_numbers):
                itypt = 0
                for i in range(num_atom_prim):
                    if inums == prim_numbers[i]:
                        # prim_positions.append(p_p[i])
                        # prim_numbers.append(p_numbers[i])
                        itypt += 1
                prim_typt.append(itypt)
        else:
            prim_lat = None
            prim_positions = None
            prim_numbers = None
            prim_typt = None
            num_atom_prim = None
    else:
        prim_lat = None
        prim_positions = None
        prim_numbers = None
        prim_typt = None
        num_atom_prim = None

    return ((num_spg, symbol_spg),
            (refine_lat, refine_pos, refine_typt, num_atoms_brv, refine_numbers),
            (prim_lat, prim_positions, prim_typt, num_atom_prim))

def plot(struct, npop):
    enth = [item[0] for item in struct]
    t = len(enth)
    nge = t / npop
    mod = t % npop

    Genth = []
    for i in range(int(nge)):
        Genth.append(enth[i*npop : (i+1)*npop])
    if mod != 0:
        Genth.append(enth[nge*npop:])

    dp = []
    for item in Genth:
        dp.append(min(item))

    Dp = []
    mine = dp[0]
    for item in dp:
        if item < mine:
            mine = item
        Dp.append(mine)

    f = open('plot.dat', 'w')
    for i in range(len(Dp)):
        f.write("%4d    %12.5f\n" % (i+1, Dp[i]))
    f.close()

    try:
        import matplotlib.pyplot as plt
    except:
        print("No matplotlib")
        if os.system('which gnuplot >/dev/null 2>/dev/null') == 0:
            os.system('''gnuplot << EOF\n set terminal dumb;\n plot 'plot.dat' w l;\n ''')
        return 0

    x = [i+1 for i in range(len(Dp))]
    y = Dp[:]

    plt.plot(x, y, 'bo-')
    plt.xlabel('Generation')
    plt.ylabel('Enthalpy (eV)')
    plt.show()


    return 0

def Zoutput(structure, options, num_proce, prec_pool, is_refine, is_prim, hard, fdir, d2, cl, norefine, hm, bg, xrd, lsur):
    global name_ele
    global prec

    OutputData = []
    # for i in range(num_proce):
    spg_old = 0
    ene_old = 610612509
    ip = 1
    for i in range(len(structure)):
        if ip > num_proce  :
            break
        # print i
        outputdata = []
        outputdata.append(structure[i][0]) # 0 enthalpy
        outputdata.append(structure[i][2]) # 1 nstruct
        outputdata.append(structure[i][3]) # 2 gen type
        outputdata.append(structure[i][1]) # 3 hardness
        outputdata.append(structure[i][5]) # 4 polar
        outputdata.append(structure[i][6]) # 5 bandgap
        outputdata.append(structure[i][7]) # 6 xrddiff

        if norefine:
            pass
        else:
            for prec in prec_pool:
                pdir = fdir + 'dir_' + str(prec)
                if options.is_wcif or options.is_wvasp:
                    os.system('mkdir -p ' + pdir)
                if abs(structure[i][0] - 610612509.) > 0.1:
                    (spgdata, recell, primcell) = findsym(structure[i][4], prec, is_refine, is_prim)
                    if spgdata[0] == spg_old and abs(outputdata[0]-ene_old) < options.rme :
                        recell = None
                        break
                        
                else:
                    spgdata = (0, 'NULL')
                    recell = None
                    break

                outputdata.append(spgdata)
                if options.is_wcif and recell is not None:
                    write_cif(pdir, recell, ip, spgdata, False, False)
                    if options.is_prim and primcell[0] is not None:
                        write_cif(pdir, primcell, i+1, spgdata, False, True)
                if options.is_wvasp and recell is not None:
                    write_vasp(pdir, recell, ip, spgdata, False, False)
                    if options.is_prim and primcell[0] is not None:
                        write_vasp(pdir, primcell, ip, spgdata, False, True)
            if not recell :
                continue

        if norefine:
            options.is_origin = True
            spgdata = (1, 'P1')

        if options.is_origin:
            if options.is_wcif or options.is_wvasp or options.is_xyz:
                pdir = fdir + 'dir_origin'
                os.system('mkdir -p ' + pdir)
            if options.is_xyz:
                try:
                    write_xyz(pdir, structure[i][4], ip, spgdata, True, False)
                except:
                    pass
            if options.is_wcif:
                try:
                    write_cif(pdir, structure[i][4], ip, spgdata, True, False)
                except:
                    pass
            if options.is_wvasp:
                try:
                    write_vasp(pdir, structure[i][4], ip, spgdata, True, False)
                except:
                    pass


        ip += 1
        spg_old = spgdata[0]
        ene_old = outputdata[0]
        OutputData.append(outputdata)


    if options.is_dft:
        fw = open(fdir + 'Analysis_Dft_Output.dat', 'w')
    else:
        fw = open(fdir + 'Analysis_Output.dat', 'w')
    if options.is_gt:
        fw.write("        No.     GenType    Enthalpy")
    else:
        fw.write("        No.      Enthalpy")
    if hard: fw.write("    Hardness")
    if hm:   fw.write("       Polar")
    if bg:   fw.write("    Diff_gap")
    if xrd:   fw.write("    XRD_Diff")
    if norefine:
        pass
    else:
        for i in range(len(prec_pool)):
            fw.write("%15g" % prec_pool[i])
    fw.write("\n")

    for i in range(len(OutputData)):
        if abs(OutputData[i][0] - 610612509.) > 0.1:
            if options.is_gt :
                fw.write("%4d (%4d)  %10s%12.5f" % (i+1, OutputData[i][1], OutputData[i][2], OutputData[i][0]))
            else:
                fw.write("%4d (%4d)  %12.5f" % (i+1, OutputData[i][1], OutputData[i][0]))
        else:
            if options.is_gt:
                fw.write("%4d (%4d)  %10s%12s" % (i+1, OutputData[i][1], OutputData[i][2], 'NULL'))
            else:
                fw.write("%4d (%4d)  %12s" % (i+1, OutputData[i][1], 'NULL'))
        if hard:
            fw.write("%12.5f" % (OutputData[i][3] * -1))
        if hm:
            fw.write("%12.5f" % (OutputData[i][4]))
        if bg:
            fw.write("%12.5f" % (OutputData[i][5]))
        if xrd:
            fw.write("%12.5f" % (OutputData[i][6]))

        if norefine:
            pass
        else:
            for sspg in OutputData[i][7:]:
                fw.write("%15s" % (sspg[1].strip() + "(" + str(sspg[0]) + ")"))
        fw.write("\n")
    fw.close()

def output(structure, options, num_proce, num_neb):
    global name_ele
    global prec

    spgdata = (1, 'P1')
    OutputData = []
    fdir = './'
    pdir = fdir + 'dir_' + 'ts'
    os.system('mkdir -p ' + pdir)
    for  i in range(num_proce):
        outputdata = []
        outputdata.append(structure[i][1]) # 0 nstruct
        outputdata.append(structure[i][0]) # 1 pts
        energy=[]
        for j in range(num_neb):
            energy.append(structure[i][2][j][4])

        ids=energy.index(max(energy))
        ts_struct=structure[i][2][ids][0:4]
        if options.is_xyz:
            try:
                write_xyz(pdir, ts_struct, i+1, spgdata, True, False)
            except:
                pass
        if options.is_wcif:
            try:
                write_cif(pdir, ts_struct, i+1, spgdata, True, False)
            except:
                pass
        if options.is_wvasp:
            try:
                write_vasp(pdir, ts_struct, i+1, spgdata, True, False)
            except:
                pass
        tsdir = pdir + '/' + str(i+1)+ '_path'
        os.system('mkdir -p ' + tsdir)
        for  j in range(num_neb):
            ts_struct=structure[i][2][j][0:4]
            if options.is_xyz:
                try:
                    write_xyz(tsdir, ts_struct, j+1, spgdata, True, False)
                except:
                    pass
            if options.is_wcif:
                try:
                    write_cif(tsdir, ts_struct, j+1, spgdata, True, False)
                except:
                    pass
            if options.is_wvasp:
                try:
                    write_vasp(tsdir, ts_struct, j+1, spgdata, True, False)
                except:
                    pass

        OutputData.append(outputdata)


    if options.is_dft:
        fw = open(fdir + 'Analysis_Dft_Output.dat', 'w')
    else:
        fw = open(fdir + 'Analysis_Output.dat', 'w')

    fw.write("        No.     Potential_Barrier")
    fw.write("\n")

    for i in range(len(OutputData)):
        if abs(OutputData[i][1] - 610612509.) > 0.1:
            fw.write("%4d (%4d)  %12.5f" % (i+1, OutputData[i][0], OutputData[i][1]))
        else:
            fw.write("%4d (%4d)  %12s" % (i+1, OutputData[i][0], 'NULL'))
        fw.write("\n")
    fw.close()

def vsckit_old(structure, vsce, name_ele, options, prec_pool, is_refine, is_prim, hard, cl, norefine, hm, bg, xrd):

    os.system('rm -rf dir_* > /dev/null 2> /dev/null')
    #structure
    #Bstruct = sorted(structure, key = lambda x : x[4][2]) # sorted with typt
    Bstruct = sorted(structure, key = lambda x: float(x[4][2][0])/float(sum(x[4][2][:])))
    #TYPT = [x[4][2] for x in Bstruct]
    TYPT = [float(x[4][2][0])/float(sum(x[4][2][:])) for x in Bstruct]
    TYPT2 = [float(x[4][2][1])/float(sum(x[4][2][:])) for x in Bstruct]
    print('TYPT', TYPT, len(TYPT))
    print('TYPT2', TYPT2, len(TYPT2))
    x = TYPT[0]
    Typt = []
    Typt.append(x)
    x2 = TYPT2[0]
    Typt2 = []
    Typt2.append(x2)
    for item, item2 in zip(TYPT, TYPT2):
        if abs(item - x) > 1e-4 and abs(item2 - x2) > 1e-4 :
            x = item
            Typt.append(x)
            x2 = item2
            Typt2.append(x2)
    print(Typt, len(Typt))
    print(Typt2, len(Typt))
    sys.exit(0)
    BStruct = []
    for item, item2 in zip(Typt, Typt2):
        bst = []
        for x in Bstruct:
            if abs( float(x[4][2][0])/float(sum(x[4][2][:])) - item ) < 1e-4 and \
                    abs( float(x[4][2][1])/float(sum(x[4][2][:])) - item2 ) < 1e-4:
                # print x[4][2][:]
                # print float(x[4][2][0])/float(sum(x[4][2][:])), item
                #print float(x[4][2][0])/float(sum(x[4][2][:]))-item
                bst.append(x)
        BStruct.append(bst)

    # print Typt[1]
    # print BStruct[1][0]
    # sys.exit(0)

    for i in range(len(BStruct)):
        for j in range(len(BStruct[i])):
            e1 = 0.
            for k in range(len(vsce)):
                e1 += vsce[k]*BStruct[i][j][4][2][k]
            forme = (BStruct[i][j][0] * BStruct[i][j][4][3] - e1) / BStruct[i][j][4][3]
            BStruct[i][j].append(forme)

    # print Typt[3]
    # print BStruct[3][1]
    # print len(BStruct[3][1])
    # sys.exit(0)

    Ydata = []
    Xdata = []
    xlabel = []
    Xlab = []
    for i in range(len(BStruct)):

        SBS = sorted(BStruct[i], key = lambda x : x[-1])
        # xlabel = float(Typt[i][-1]) / sum(Typt[i])
        # pdata.append(xlabel)
        #print Typt[i], SBS[0][5]
        #if len(name_ele) > 2:
        #    Xdata.append(i+1)
        #else:
        #    Xdata.append(float(Typt[i][-1]) / sum(Typt[i]))
        if  abs(SBS[0][-1] - 610612509) > 100.:
            Xdata.append(Typt[i])
            Ydata.append(SBS[0][-1])
        # print Ydata
        vdir = '%.4f' % Typt[i]
        vdir2 = '%.4f' % Typt2[i]
        xt = 'x'
        # for j in range(len(Typt[i])):
        #     vdir += name_ele[j] + str(Typt[i][j])
        #     xt += name_ele[j] + '_{' + str(Typt[i][j]) + '}'

        # print vdir
        xxt = r'$' + xt + r'$'
        if  abs(SBS[0][-1] - 610612509) > 100.:
            xlabel.append(xxt)
            Xlab.append(vdir)
        if options.is_all:
            num_proce = len(SBS)
        else:
            if options.num_proce < len(SBS):
                num_proce = options.num_proce
            else :
                num_proce = len(SBS)
        fdir = 'dir_' + vdir + '_' + vdir2 + '/'
        os.system('mkdir -p ' + fdir)
        lsur='F'
        Zoutput(SBS, options, num_proce, prec_pool, is_refine, is_prim, hard, fdir, False, cl, norefine, hm, bg, xrd, lsur)
    if len(name_ele) == 2:
        if Xdata[0] < Xdata[-1]:
            Xdata = [0.] + Xdata + [1.]
        else:
            Xdata = [1.] + Xdata + [0.]
        Ydata = [0] + Ydata + [0]
        xlabel = [r'$'+name_ele[0]+r'$'] + xlabel + [r'$'+name_ele[1]+r'$']
        Xlab = [name_ele[0]] + Xlab + [name_ele[1]]

    f = open('Convexhull.dat', 'w')
    for i in range(len(Xdata)):
        if len(name_ele) == 2:
            f.write('%5.4f  %12.5f  %15s\n' % (Xdata[i], Ydata[i], Xlab[i]))
        else:
            f.write('%4d  %12.5f  %15s\n' % (Xdata[i], Ydata[i], Xlab[i]))
    f.close()

    if options.is_pch:

        try:
            import matplotlib.pyplot as plt
        except:
            print("No matplotlib")
            return 0


        plt.plot(Xdata, Ydata, 'bo--')
        plt.ylabel('Enthalpy (eV)')
        plt.xticks(Xdata,tuple(xlabel))
        plt.title('Convex hull')
        if len(name_ele) == 2:
            plt.axis([0, 1, min(Ydata)-abs(min(Ydata)*.05), max(Ydata)+abs(max(Ydata)*.05)])
        else:
            plt.axis([0, len(Xdata)+1, min(Ydata)-abs(min(Ydata)*.05), max(Ydata)+abs(max(Ydata)*.05)])
        plt.show()

def vsckit(structure, vsce, name_ele, options, prec_pool, is_refine, is_prim, hard, cl, norefine, hm, bg, xrd):
    import itertools
    import pandas as pd

    os.system('rm -rf dir_* > /dev/null 2> /dev/null')
    Bstruct = structure
    #structure
    #Bstruct = sorted(structure, key = lambda x : x[4][2]) # sorted with typt
    #Bstruct = sorted(structure, key = lambda x: float(x[4][2][0])/float(sum(x[4][2][:])))
    # Bstruct = sorted(structure, key = lambda x: (float(x[4][2][0])/float(sum(x[4][2][:])),float(x[4][2][0])/float(sum(x[4][2][:]))))
    TYPT = [x[4][2][:] for x in Bstruct]
    ENTH = [x[0] for x in Bstruct]
    StruNum = [x[2] for x in Bstruct]
    tmp_list = []
    for one in TYPT :
        if one not in tmp_list :
            tmp_list.append(one)
    # print len(tmp_list)
    typt_list = []
    for one in tmp_list :
        one = gcdnum(one)
        if one not in typt_list :
            typt_list.append(one)

    print(len(typt_list))
    typt_dict = {}
    BStruct = []
    for i in range(len(typt_list)):
        temp_int = list(map(int, typt_list[i]))
        str1 = '_'.join(map(str, temp_int))
        typt_dict[str1] = i
        BStruct.append([])
    # print BStruct
    for x in Bstruct:
        typt = x[4][2][:]
        typt = gcdnum(typt)
        typt = list(map(int, typt))
        str1 = '_'.join(map(str, typt))
        BStruct[int(typt_dict[str1])].append(x)


    for i in range(len(BStruct)):
        # print len(BStruct[i])
        for j in range(len(BStruct[i])):
            e1 = 0.
            for k in range(len(vsce)):
                e1 += vsce[k]*BStruct[i][j][4][2][k]
            # forme = (BStruct[i][j][0] * BStruct[i][j][4][3] - e1) / BStruct[i][j][4][3]
            forme = BStruct[i][j][0]
            BStruct[i][j].append(forme)

    # for i in range(len(BStruct)):
    for key in typt_dict :
        i = typt_dict[key]
        SBS = sorted(BStruct[i], key = lambda x : x[-1])
        if options.is_all:
            num_proce = len(SBS)
        else:
            if options.num_proce < len(SBS):
                num_proce = options.num_proce
            else :
                num_proce = len(SBS)
        fdir = 'dir_' + key + '/'
        print(fdir, len(SBS))
        os.system('mkdir -p ' + fdir)
        lsur='F'
        Zoutput(SBS, options, num_proce, prec_pool, is_refine, is_prim, hard, fdir, False, cl, norefine, hm, bg, xrd, lsur)

    # convexhull.csv
    kformula = []
    kenth = []
    knum = []
    cdata = zip(TYPT, ENTH, StruNum)
    for tformula_num, tenth, Num in cdata:
        tformula_num_float = list(map(int, tformula_num))
        tformula_num_str   = list(map(str, tformula_num_float))
        tformula =  ''.join([a + b for a,b in zip(name_ele, tformula_num_str)])
        kformula.append(tformula)
        kenth.append(tenth)
        knum.append(Num)

    df = pd.DataFrame({'Number':knum,'formula':kformula,'enthalpy':kenth,})
    df.to_csv('convexhull.csv',index=False)

def gcdnum(xl):
    m = min(xl)
    for i in range(m,1,-1):
        if not any([x%i for x in xl]):
            return [x/i for x in xl]
    return xl

def run():
    global name_ele
    global prec

    pinfinty = float(1e300)
    einfinty = float(-1e300)

    evolution = False

    parser = OptionParser()
    parser.set_defaults( is_all = False,
                         is_wcif = False,
                         is_wvasp = False,
                         is_xyz   = False,
                         is_prim = False,
                         is_origin = False,
                         is_hard = False,
                         is_gt = False,
                         is_v = False,
                         is_pch = False,
                         is_popt = False,
                         is_dft = False,
                         is_nosym = False,
                         is_ini = False,
                         plotch = False,
                         plotevo = False,
                         num_proce = 50,
                         prec = 0.1,
                         mprec = None,
                         minh = einfinty,
                         maxh = pinfinty,
                         minp = einfinty,
                         maxp = pinfinty,
                         mine = einfinty,
                         is_vsc_split = False, 
                         rme = 0.001, 
                         maxe = pinfinty )
    parser.add_option("-a", dest="is_all", action="store_true")
    parser.add_option("--cif", dest="is_wcif", action="store_true")
    parser.add_option("--pos", "--vasp", dest="is_wvasp", action="store_true")
    parser.add_option("--xyz", dest="is_xyz", action="store_true")
    parser.add_option("--pri", "--prim", "--primitive", dest="is_prim", action="store_true")
    parser.add_option("--ori", "--orig",  "--origin", dest="is_origin", action="store_true")
    parser.add_option("--gt", dest="is_gt", action="store_true")
    parser.add_option("--nosym", dest="is_nosym", action="store_true")
    parser.add_option("-n", dest="num_proce", type="int")
    parser.add_option("-t", "--tolerance", dest="prec", type="float")
    parser.add_option("-m", "--multi-tolerance", dest="mprec", action="store", type="string")
    parser.add_option("-p", dest="is_plot", action="store_true")
    parser.add_option("--hard", dest="is_hard", action="store_true")
    parser.add_option("--minh", dest="minh", type="float")
    parser.add_option("--maxh", dest="maxh", type="float")
    parser.add_option("--mine", dest="mine", type="float")
    parser.add_option("--maxe", dest="maxe", type="float")
    parser.add_option("--minp", dest="minp", type="float")
    parser.add_option("--maxp", dest="maxp", type="float")
    parser.add_option("--pch", dest="is_pch",  action="store_true")
    parser.add_option("--popt", dest="is_popt",  action="store_true")
    parser.add_option("--dft", dest="is_dft",  action="store_true")
    parser.add_option("--ini", dest="is_ini", action="store_true")
    parser.add_option("--ts", dest="is_ts", action="store_true")
    parser.add_option("--split", dest="is_vsc_split", action="store_true")
    parser.add_option("--rme","--remove-energy",dest="rme", type="float")
    parser.add_option("--plotch","--plot_convexhull",dest="plotch", action="store_true")
    parser.add_option("--plotevo","--plot_evo",dest="plotevo", action="store_true")


    parser.add_option("-v", "--version", dest="is_v", action="store_true")

    (options, args) = parser.parse_args()

    if options.is_v:
        print("Version: 2016.10.10")
        exit(0)

    if options.is_vsc_split :
        SplitVSCStruct()
        exit(0)


    (sysname, name_ele, npop, mol, hard, vsc, vsce, d2, cl, hm, lsur, bg, xrd, ts, num_neb) = readinput()

    if options.plotch:
        from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
        from pymatgen.core.composition import Composition
        import pandas as pd
        
        ini_entries = []
        convexhull_data = pd.read_csv('./convexhull.csv', header=0, sep=r',')
        
        for idx, row in convexhull_data.iterrows():
            #if row['form_energy'] <= 0.02:
            comp = Composition(row['formula'])
            num_atoms = comp.num_atoms
            enth_per_atom = row['enthalpy']
            enth = enth_per_atom * num_atoms
            entry_id = row['Number']
            _entry = PDEntry(comp, enth)
            _entry.entry_id = entry_id
            # _entry.spg = row['spg']
            # _entry.p = row['P']
            # _entry.nsw = row['NSW']
            # _entry.epa = enth_per_atom
            ini_entries.append(_entry)              

        for i in range(len(name_ele)):
            comp = name_ele[i]
            enth = vsce[i]
            entry_id = name_ele[i]
            _entry = PDEntry(comp, enth)
            _entry.entry_id = entry_id
            ini_entries.append(_entry)              
            
        ini_pd = PhaseDiagram(ini_entries)
        plotter = PDPlotter(ini_pd, show_unstable=0.050, backend='plotly')
        plotter.get_plot().write_image('convexhull.png')     
        
        
        e_above_hull = []
        form_energy = []
        names = []
        nsws = []
        formulas = []
        enths = []
        spgs = []
        ps = []
        threshold_start = 0.00
        #threshold = 0.05
        
        for entry in ini_entries:
            ebh = ini_pd.get_e_above_hull(entry, )
            fme = ini_pd.get_form_energy_per_atom(entry)
            name = entry.entry_id
            # spg = entry.spg
            # p = entry.p
            formula = entry.composition.formula.replace(".", "")
            #if threshold_start <= ebh <= threshold:
            e_above_hull.append(ebh)
            form_energy.append(fme)
            names.append(name)
            formulas.append(''.join(  (''.join(str(formula))).split()  ))
            # spgs.append(spg)
            # ps.append(p)
            # enths.append(enth)
        
        
        df = pd.DataFrame(
            {
                'formula':formulas,
                'e_above_hull':e_above_hull,
                # 'spg':spgs,
                # 'p':ps,
                'form_energy':form_energy,
                'name':names,
                }
            )
        df.to_csv('./e_above_hull_50meV.csv', index=False, sep=' ')
        sys.exit()

    if options.plotevo:
        import matplotlib.pyplot as plt
        import scienceplots
        
        plt.style.use(['science', 'ieee'])
        
        
        # data
        pso_sor_list = glob.glob('./pso_sor_*')
        x = []
        yy = []
        y = [np.loadtxt('pso_sor_%d'%(t)).tolist() for t in range(1, len(pso_sor_list)+1)]
        for evo_num, ene_list in enumerate(y):
            for ene in ene_list:
                yy.append(ene)
                x.append(evo_num+1)

        temp = list(set(yy))
        temp.sort()
        max_y = temp[-1] if temp[-1] < 1000 else temp[-2]
        min_y = temp[0] if temp[0] > -1000 else temp[1]

        # colormap
        cm = plt.cm.get_cmap('rainbow')
        
        # fig, ax
        fig, ax = plt.subplots(figsize=(10, 4))
        
        # ax.scatter(x, yy, c=yy, cmap=cm, vmin=2.3, vmax=2.6, s=5)
        sc = ax.scatter(x, yy, c=yy, cmap=cm, vmin=min_y, vmax=max_y, s=10)
        
        # beauty
        ax.set_xlabel('Structure Evolution Step', fontsize=25, )
        ax.set_ylabel('Energy (eV/atom)', fontsize=25, )
        # ax.set_xticks([0, 10, 20, 30, 40, 50],)
        # ax.set_yticks([-3.6, -3.5, -3.4, -3.3, -3.2],)

        # ax.tick_params(width=5, labelsize=10)
        ax.set_ylim((min_y, max_y))
        
        plt.tick_params(labelsize=20)
        
        # colorbar
        cbar = fig.colorbar(sc)
        # cbar.set_ticks(ticks=[-3.6, -3.5, -3.4, -3.3, -3.2])
        # cbar.ax.tick_params(labelsize=20)
        cbar.ax.minorticks_off()
        cbar.outline.set_visible(False)
        
        
        # spines width
        width = 2
        ca = plt.gca()
        
        ca.spines['bottom'].set_linewidth(width)
        ca.spines['left'].set_linewidth(width)
        ca.spines['top'].set_linewidth(width)
        ca.spines['right'].set_linewidth(width)
        
        plt.tight_layout()
        
        fig.savefig('evo.png', dpi=350)
        sys.exit()

    # struct = parseStruct()
    if options.is_ini:
        struct = parseIni()
    elif lsur:
        struct = parseSur()
    elif ts:
        struct = parseTS()
    elif options.is_dft:
        struct = parseOpt('dft')
    elif options.is_popt:
        struct = parseOpt('opt', hm)
    else:
        if os.path.exists('struct.dat'):
            try:
                struct = parseStruct()
            except:
                struct = parseOpt('opt', hm)
        else:
            struct = parseOpt('opt', hm)

    # struct [0] enthalpy; [1] hardness; [2] nstruct; [3] gtype; [4] cell; [5] polar; [6] bandgap

    structure = []

    if cl or lsur or options.is_nosym or ts:
        norefine = True
    else:
        norefine = False

    if hard:
        for x in struct:
            if x[0] < options.maxe and x[0] > options.mine and \
               x[1] < -1*options.minh and x[1] > -1*options.maxh:
                structure.append(x)
    elif hm:
        for x in struct:
            if x[0] <= options.maxe and x[0] >= options.mine and \
               x[5] <= options.maxp and x[5] >= options.minp:
                structure.append(x)
    else:
        structure = struct[:]


    if options.is_all:
        num_proce = len(structure)
    else:
        if options.num_proce < len(structure):
            num_proce = options.num_proce
        else :
            num_proce = len(structure)

    if options.mprec is None:
        prec_pool = [options.prec]
    else:
        try:
            prec_pool = list(map(float, options.mprec.split()))
        except:
            print('warning: parse multi tolerance error')
            prec_pool = [options.prec]

    if options.is_wvasp or options.is_wcif:
        is_refine = True
        is_prim = options.is_prim
        os.system('rm -rf dir_* > /dev/null 2> /dev/null')
    else:
        is_refine = False
        is_prim = False

    if options.is_plot:
        plot(struct, npop)

    for item in structure :
        if item[0] > options.maxe or item[0] < options.mine:
            item[0] = 610612509

    if vsc:
        vsckit(structure, vsce, name_ele, options, prec_pool, is_refine, is_prim, hard, cl, norefine, hm, bg, xrd)
    elif ts:
        structure.sort(key=lambda x:x[0])
        output(structure, options, num_proce, num_neb)
    else:
        if options.is_hard:
            structure.sort(key=lambda x:x[1])
        else:
            structure.sort(key=lambda x:x[0])
        if hm:
            structure.sort(key=lambda x:x[5], reverse=True)
        if bg:
            structure.sort(key=lambda x:x[6])
        if xrd:
            structure.sort(key=lambda x:x[7])

        fdir = './'
        Zoutput(structure, options, num_proce, prec_pool, is_refine, is_prim, hard, fdir, d2, cl, norefine, hm, bg, xrd, lsur)



if __name__ == '__main__':
    run()
    # print gcdnum([4, 4, 32])
    # print gcdnum([4, 6])


