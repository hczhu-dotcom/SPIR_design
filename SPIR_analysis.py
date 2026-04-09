#这个的脚本的目的，是从不同数据库中收集nanobody,并进行汇总，准备进行之后的虚拟筛选
import pandas as pd
import os
import numpy as np
import datetime
import re
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.Seq import Seq
import Bio
import random
#from pymol import cmd
import pymol
import json
import argparse
import shutil
from pymol import stored   #这里不会报错
from pymol.cmd import util
from Bio import SeqIO
from pymol import stored   #这里不会报错
from pymol.cmd import util
import sys
import math
import plotnine
from plotnine import ggplot, aes, ggsave
import subprocess
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB import PDBParser, Superimposer, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Polypeptide import is_aa
from Bio.Data import CodonTable
import random
from lxml import etree
import warnings


import traceback

fr1_seq = 'QVQLVESGGGLVQAGGSLRLSCAAS'
fr2_seq = 'MGWFRQAPGKEREFVAA'
fr3_seq = 'YYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC'
fr4_seq = 'WGQGTQVTVS'
#‘AQVQLVESGGGLVQAGGSLRLSCAASGRTFSEYNMGWFRQAPGKEREFVAAIRSSGTTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAMSRVDTDSPAFYDYWGQGTQVTVSK’
#创建一个nano body的class.
class My_nanobody():

    def __init__(self,sequence:str = 'AQVQLVESGGGLVQAGGSLRLSCAASGRTFSEYNMGWFRQAPGKEREFVAAIRSSGTTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAMSRVDTDSPAFYDYWGQGTQVTVSK',
                 stringent = 1,silent = 0,chain_name = 'H'):
        self.sequence = sequence.upper()
        self.chain_name = chain_name
        #QVQLQESGGGSVQVGGSLRVACAASGDTFSGYLAAWFRQAPGKGREGVAAINSKRHTTSYADSVKGRFTISKDNADNIMYLEMNSLKPEDTAIYYCAAADAIGLAEYWSTPTLSAARYKYWGQGTQVTVSS
        if stringent == 0:
            key_fr1 = re.search(r'\w{8}[PAVST]G[GED]SL[RTKI][LVAT][SAT]C',self.sequence)
            key_fr2 = re.search(r'\w\wRQ[AVG]PG[KQ]\w{3}[GWLRYELFA][VIL]\w',self.sequence)
            key_fr3 = re.search(r'Y\w[DGNVEA][SAP][VWAMILG][KRE][GSA]R[FSA][TSIVA][IVTLN][STAYNP]\w[DESVHN]\w\w[KRQEHANLGD][NKSRDGTAI][TLKVPMRIMA][VLAIGM][YFSHAT][LVRFM][QHYKRE][MLVIDR][NDSYVERT]\w[LVWPR][KERNQ][PSTAIL][ED][DE][TKSAP][AGT][VMLITLR]Y[YFHWRS]\w', self.sequence)
            key_fr4 = re.search(r'[GSAD][QKRDE]G[TA][QLP]V', self.sequence)
        elif stringent == 1:
            key_fr1 = re.search(r'[VL]QL[VQ]E[ST]GGG[LS]VQ[PAV]G[GE]SLRLSC',self.sequence)
            key_fr2 = re.search(r'W[FYV]RQ[AVG]PG[KQ][EGQ][RL]E[GWLFA][VIL][AS]',self.sequence)
            key_fr3 = re.search(r'Y[AS][DEA]SV[KRE]GRF[TA][IV]S[RK]D[ND][AS]K[NKS][TIMA][VLAM][YFS]L[QE][MLI][NDST][SNDR]L[KERNQ][PSTA][ED]DT[AG][VMLI]Y[YFS]C', self.sequence)
            key_fr4 = re.search(r'G[QKR]GT[QLP]VTVS', self.sequence)

        elif stringent == 2:
            key_fr1 = re.search(r'\w{6}GGG\w{6}SL\wLSC',self.sequence)
            key_fr2 = re.search(r'W[FYV]RQ\w\w\w\w\w[RLG][EDRT]\w[VIL][AS]',self.sequence)
            key_fr3 = re.search(r'[YL]\w\w[ASF][VAG][KRE][GD]RF[TSA][IV]S\w\w\w[ATS]\w{5}L\w{12}Y\wC', self.sequence)
            key_fr4 = re.search(r'G[QKR]GT[QLRP]VTVS', self.sequence)
        elif stringent == 3:
            key_fr1 = re.search(r'\w{6}GGG\w{6}[SR]L\wL[SA]C',self.sequence)
            key_fr2 = re.search(r'W[FYV]RQ\w\w\w\w\w[RLGK][EDRT]\w[VIL][AST]',self.sequence)
            key_fr3 = re.search(r'[YL]\w\w[ASF][VAG][KRE][GD]R[FT][TSA][IVA]S\w\w\w[ATSD]\w{5}L\w{12}Y\wC', self.sequence)
            key_fr4 = re.search(r'G[QKR]GT[QLRPIK]VTV[SN]', self.sequence)



        if key_fr1 == None or key_fr2 == None or key_fr3 == None or key_fr4 == None:
            if silent == 0:
                if key_fr1 == None:
                    print('no_fr1_sequence')
                if key_fr2 == None:
                    print('no_fr2_sequence')
                if key_fr3 == None:
                    print('no_fr3_sequence')
                if key_fr4 == None:
                    print('no_fr4_sequence')
            self.itis_nanobody = False
            #return 'the sequence my not a nanobody'


        elif key_fr1 != None and key_fr2 != None and key_fr3 != None and key_fr4 != None:
            self.itis_nanobody = True

            self.cdr1 = self.sequence[key_fr1.span()[1] + 3:key_fr2.span()[0] - 2]
            self.cdr2 = self.sequence[key_fr2.span()[1] + 1:key_fr3.span()[0] - 1]
            self.cdr3 = self.sequence[key_fr3.span()[1] + 0:key_fr4.span()[0] - 1]
            if silent == 0:
                print('cdr1',self.cdr1,'cdr2',self.cdr2,'cdr3',self.cdr3)

            self.len_cdr1 = len(self.cdr1)
            self.len_cdr2 = len(self.cdr2)
            self.len_cdr3 = len(self.cdr3)

    def is_nanobody(self):
        #判断输入序列是否可能是一个nanobody

        return self.itis_nanobody


    def standized_nb(self):
        if self.itis_nanobody == 1:
            std_seq = fr1_seq+self.cdr1+fr2_seq+self.cdr2+fr3_seq+self.cdr3+fr4_seq
            self.sequence = std_seq
            return std_seq
        else:
            return ''

    def mut_cdr1(self):
        #将cdr1区域变成X
        self.mut = ''
        for i in range(len(self.cdr1)):
            self.mut += 'X'
        return fr1_seq+self.mut+fr2_seq+self.cdr2+fr3_seq+self.cdr3+fr4_seq

    def mut_cdr2(self):
        #将cdr1区域变成X
        self.mut = ''
        for i in range(len(self.cdr2)):
            self.mut += 'X'
        return fr1_seq+self.cdr1+fr2_seq+self.mut+fr3_seq+self.cdr3+fr4_seq

    def mut_cdr3(self):
        #将cdr1区域变成X
        self.mut = ''
        for i in range(len(self.cdr3)):
            self.mut += 'X'
        return fr1_seq+self.cdr1+fr2_seq+self.cdr2+fr3_seq+self.mut+fr4_seq

    def mut_cdr123(self):
        self.mut1 = ''
        self.mut2 = ''
        self.mut3 = ''
        for i in range(len(self.cdr1)):
            self.mut1 += 'X'
        self.mut2 = ''
        for i in range(len(self.cdr2)):
            self.mut2 += 'X'
        self.mut3 = ''
        for i in range(len(self.cdr3)):
            self.mut3 += 'X'
        self.mut_cdr123_sequence = fr1_seq + self.mut1 + fr2_seq + self.mut2 + fr3_seq + self.mut3 + fr4_seq
        #return fr1_seq + self.mut1 + fr2_seq + self.mut2 + fr3_seq + self.mut3 + fr4_seq
        return self.mut_cdr123_sequence

    def obtain_cdr_border(self,cdr = 'cdr3'):
        cdr_begin = 0
        cdr_end = 0
        if cdr == 'cdr1':
            mut_cdr_seq = self.mut_cdr1()
        if cdr == 'cdr2':
            mut_cdr_seq = self.mut_cdr2()
        if cdr == 'cdr3':
            mut_cdr_seq = self.mut_cdr3()

        for i in range(len(mut_cdr_seq)):
            if mut_cdr_seq[i] == 'X':
                cdr_begin = i
                break

        for j in range(len(mut_cdr_seq)):
            if mut_cdr_seq[j] != 'X' and j> cdr_begin:
                cdr_end = j-1
                break

        return cdr_begin, cdr_end


    def obtain_the_fixed_residues_cdr123(self,with_chain = True):
        #用于提供参数 给proteinmpnn和ligandmpnn
        fixed_residues = ''
        self.mut_cdr123_sequence = self.mut_cdr123()
        for i in range(len(self.mut_cdr123_sequence)):
            if self.mut_cdr123_sequence[i] != 'X':
                if self.sequence.startswith('M') or self.sequence.startswith('A'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 2) + ' '
                    else:
                        fixed_residues += str(i + 2) + ' '
                elif self.sequence.startswith('Q') or self.sequence.startswith('E'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 1) + ' '
                    else:
                        fixed_residues += str(i + 1) + ' '
        fixed_residues = fixed_residues[:-1]
        #print('fixed_residues',fixed_residues)
        return fixed_residues

    def obtain_the_fixed_residues_cdr3(self,with_chain = True):
        # 用于提供参数 给proteinmpnn和ligandmpnn
        fixed_residues = ''
        self.mut_cdr3_sequence = self.mut_cdr3()
        for i in range(len(self.mut_cdr3_sequence)):
            if self.mut_cdr3_sequence[i] != 'X':
                if self.sequence.startswith('M') or self.sequence.startswith('A'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 2) + ' '
                    else:
                        fixed_residues += str(i + 2) + ' '
                elif self.sequence.startswith('Q') or self.sequence.startswith('E'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 1) + ' '
                    else:
                        fixed_residues += str(i + 1) + ' '
        fixed_residues = fixed_residues[:-1]
        #print('fixed_residues', fixed_residues)
        return fixed_residues

    def obtain_the_fixed_residues_cdr1(self,with_chain = True):
        # 用于提供参数 给proteinmpnn和ligandmpnn
        fixed_residues = ''
        self.mut_cdr1_sequence = self.mut_cdr1()
        for i in range(len(self.mut_cdr1_sequence)):
            if self.mut_cdr1_sequence[i] != 'X':
                if self.sequence.startswith('M') or self.sequence.startswith('A'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 2) + ' '
                    else:
                        fixed_residues += str(i + 2) + ' '
                elif self.sequence.startswith('Q') or self.sequence.startswith('E'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 1) + ' '
                    else:
                        fixed_residues += str(i + 1) + ' '
        fixed_residues = fixed_residues[:-1]
        #print('fixed_residues', fixed_residues)
        return fixed_residues

    def obtain_the_fixed_residues_cdr2(self,with_chain = True):
        # 用于提供参数 给proteinmpnn和ligandmpnn
        fixed_residues = ''
        self.mut_cdr2_sequence = self.mut_cdr2()
        for i in range(len(self.mut_cdr2_sequence)):
            if self.mut_cdr2_sequence[i] != 'X':
                if self.sequence.startswith('M') or self.sequence.startswith('A'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 2) + ' '
                    else:
                        fixed_residues += str(i + 2) + ' '
                elif self.sequence.startswith('Q') or self.sequence.startswith('E'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 1) + ' '
                    else:
                        fixed_residues += str(i + 1) + ' '

        fixed_residues = fixed_residues[:-1]
        #print('fixed_residues', fixed_residues)
        return fixed_residues

        #for each##
###################################
    def obtain_the_designed_residues_cdr123(self,with_chain = True):
        #用于提供参数 给proteinmpnn和ligandmpnn
        fixed_residues = ''
        self.mut_cdr123_sequence = self.mut_cdr123()
        for i in range(len(self.mut_cdr123_sequence)):
            if self.mut_cdr123_sequence[i] == 'X':
                if self.sequence.startswith('M') or self.sequence.startswith('A'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 2) + ' '
                    else:
                        fixed_residues += str(i + 2) + ' '
                elif self.sequence.startswith('Q') or self.sequence.startswith('E'):
                    if with_chain:
                        fixed_residues += self.chain_name + str(i + 1) + ' '
                    else:
                        fixed_residues += str(i + 1) + ' '

        fixed_residues = fixed_residues[:-1]
        #print('fixed_residues',fixed_residues)
        return fixed_residues



def get_current_date_hczhu():

    #获取当前日期。
    #输出格式类似于‘20250613’
    current_date = datetime.datetime.now().date()
    if current_date.month < 10:
        current_month = '0' + str(current_date.month)
    else:
        current_month = str(current_date.month)
    if current_date.day < 10:
        current_day = '0' + str(current_date.day)
    else:
        current_day = str(current_date.day)
    current_data2 = str(current_date.year)+current_month+current_day
    return current_data2
    #‘20250613’

def pi_pi_angle(x1,x2,x3,y1,y2,y3):
    import numpy as np
    #print(x1,x2,x3,y1,y2,y3)

    B1, B2, B3 = [x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]]
    C1, C2, C3 = [x1[0] - x3[0], x1[1] - x3[1], x1[2] - x3[2]]
    #print([x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]])
    n1 = [B2 * C3 - C2 * B3, B3 * C1 - C3 * B1, B1 * C2 - C1 * B2]
    D1, D2, D3 = [y1[0] - y2[0], y1[1] - y2[1], y1[2] - y2[2]]
    E1, E2, E3 = [y1[0] - y3[0], y1[1] - y3[1], y1[2] - y3[2]]
    n2 = [D2 * E3 - E2 * D3, D3 * E1 - E3 * D1, D1 * E2 - E1 * D2]
    dot_product = np.dot(n1, n2)
    magnitude1 = np.linalg.norm(n1)  # ????????????????
    magnitude2 = np.linalg.norm(n2)
    #print(dot_product,magnitude1,magnitude2,dot_product / (magnitude1 * magnitude2))

    xx = math.acos(dot_product / (magnitude1 * magnitude2))
    degree = math.degrees(xx)
    if degree > 90:
        degree = 180 - degree
    return round(degree,2)

def pi_pi(aro_ring_info,cutoff=5.0):
    aro_resn =['HIS','PHE','TRP','TYR']
    aro_resn_list = '+'.join(aro_resn)
    chA_aro_res = f'inter_1 and resn {aro_resn_list}'
    chB_aro_res = f'inter_2 and resn {aro_resn_list}'
    chA_aro_res_list = list(set([f'resi {a.resi} and resn {a.resn}' for a in pymol.cmd.get_model(chA_aro_res).atom]))
    chB_aro_res_list = list(set([f'resi {a.resi} and resn {a.resn}' for a in pymol.cmd.get_model(chB_aro_res).atom]))
    for A in chA_aro_res_list:
        for B in chB_aro_res_list:
            A_type = A.split()[-1]
            B_type = B.split()[-1]
            #print(f'{A} and name {aro_ring_info[A_type]}')
            x = pymol.cmd.centerofmass(f'inter_1 and {A} and name {aro_ring_info[A_type]}')
            y = pymol.cmd.centerofmass(f'inter_2 and {B} and name {aro_ring_info[B_type]}')
            dist = math.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)
            if dist < cutoff:
                pymol.cmd.show('sticks',f'inter_1 and {A}')
                pymol.cmd.show('sticks', f'inter_2 and {B}')
                pymol.cmd.label(f'name CA and  inter_1 and {A}', 'oneletter+resi')
                pymol.cmd.label(f'name CA and  inter_2 and {B}', 'oneletter+resi')
                x1 = pymol.cmd.get_coords(f"{A} and name {aro_ring_info[A_type].split('+')[0]}")[0]
                x2 = pymol.cmd.get_coords(f"{A} and name {aro_ring_info[A_type].split('+')[1]}")[0]
                x3 = pymol.cmd.get_coords(f"{A} and name {aro_ring_info[A_type].split('+')[2]}")[0]
                y1 = pymol.cmd.get_coords(f"{B} and name {aro_ring_info[B_type].split('+')[0]}")[0]
                y2 = pymol.cmd.get_coords(f"{B} and name {aro_ring_info[B_type].split('+')[1]}")[0]
                y3 = pymol.cmd.get_coords(f"{B} and name {aro_ring_info[B_type].split('+')[2]}")[0]
                print(A,B)
                pi_angle = pi_pi_angle(x1, x2, x3, y1, y2, y3)
                pi_name = f'pi_{A.split()[1]}_{B.split()[1]}_{str(pi_angle)}'
                print(A, B,pi_name)


                pymol.cmd.distance(pi_name, f"inter_1 and resi {A.split()[1]} and name {aro_ring_info[A_type]}", f'inter_2 and resi {B.split()[1]} and name {aro_ring_info[B_type]}', mode=4, label=0)


def Salt_bridge(ChA,ChB,cutoff=6.0):
    selection_1 = f'{ChB} and ((resn LYS and name NZ) or (resn ARG and name NE+NH*))'
    a_charge_res = set([f'chain {a.chain} and  resi {a.resi} and resn {a.resn}' for a in pymol.cmd.get_model(f'inter_1 and resn LYS+ARG+ASP+GLU').atom])
    b_charge_res = set([f'chain {a.chain} and resi {a.resi} and resn {a.resn}' for a in pymol.cmd.get_model(f'inter_2  and resn LYS+ARG+ASP+GLU').atom])
    #print(a_charge_res)
    p_info = {'LYS':'NZ','ARG':'NE+NH*'}
    n_info = {'ASP':'OD*+OE*','GLU':'OD*+OE*'}
    for a  in a_charge_res:
        a_resn = a.split()[-1]
        for b in b_charge_res:
            b_resn = b.split()[-1]
            name = f'SB_{a.split()[4]}_{b.split()[4]}'

            if (a_resn in ['LYS','ARG']) and (b_resn in ['ASP','GLU']):
                a_resn_atom = p_info[a_resn]
                b_resn_atom = n_info[b_resn]
                #print(a,b)
                #print(a_resn_atom,b_resn_atom)
                selection_1 = f'{a} and name {a_resn_atom}'
                selection_2 = f'{b} and name {b_resn_atom}'
                cal_saltbrige(selection_1, selection_2, name, a,b,cutoff)


            elif (a_resn in ['ASP','GLU'] ) and (b_resn in ['LYS','ARG']):
                a_resn_atom = n_info[a_resn]
                b_resn_atom = p_info[b_resn]


                selection_1 = f'{a} and name {a_resn_atom}'
                selection_2 = f'{b} and name {b_resn_atom}'
                cal_saltbrige(selection_1,selection_2,name, a,b,cutoff)


def cal_saltbrige(selection1,selection2,name,a,b,cutoff):
    x = pymol.cmd.centerofmass(selection1)
    y = pymol.cmd.centerofmass(selection2)
    dist = math.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)

    if dist < cutoff:
        pymol.cmd.show('stick', a)
        pymol.cmd.show('stick', b)
        print(a,b)
        pymol.cmd.label(f'name CA and  {a}', 'oneletter+resi')
        pymol.cmd.label(f'name CA and  {b}', 'oneletter+resi')

        pymol.cmd.distance(name, selection1, selection2, mode=4, label=0)
    return name

def DA_Hbond(DA_type,select1,select2,cutoff='3.5',angle='150'):
    hbond_pairs = []
    pairs = pymol.cmd.find_pairs(select1, select2, mode=0, cutoff=cutoff)
    for pair in pairs:

        D_atom = f'index {pair[0][1]}'
        A_atom = f'index {pair[1][1]}'
        HD_atom = f'e. h and neighbor({D_atom})'
        index_1 = [a.index for a in pymol.cmd.get_model(A_atom).atom]
        index_2 = [a.index for a in pymol.cmd.get_model(HD_atom).atom]
        index_3 = [a.index for a in pymol.cmd.get_model(D_atom).atom]

        for h in index_2:
            h_angle = pymol.cmd.angle('angle', A_atom, f'index {h}', D_atom)
            #print(type, pair, h_angle)
            if h_angle > float(angle):
                D_atom_name = [f'{a.resi}-{a.name}' for a  in pymol.cmd.get_model(f'index {pair[0][1]} ').atom][0]

                A_atom_name = [f'{a.resi}-{a.name}' for a  in pymol.cmd.get_model(f'index {pair[1][1]} ').atom][0]
                #name = [f'HB-{a.resi}{a.resn}{a.chain}-{base_atom_name}-lig-{lig_atom_name}' for a in cmd.get_model(f'index {pair[0][1]} and polymer').atom][0]
                name = f'{DA_type}_{D_atom_name}_{A_atom_name}'
                pymol.cmd.show('stick', f'byres(index {pair[0][1]})')
                pymol.cmd.show('stick', f'byres(index {pair[1][1]})')
                pymol.cmd.label(f'name CA and byres(index {pair[0][1]})', 'oneletter+resi')
                pymol.cmd.label(f'name CA and  byres(index {pair[1][1]})', 'oneletter+resi')


                pymol.cmd.distance(name, f'index {index_1[0]} ', f'index {h} ', label=0)

                break
    return hbond_pairs
def ppi(pdb_name,ChA="A",ChB="B",dist=3):
    #cmd.h_add('all')
    #cmd.remove('solvent and resn SO4, NA,CL,GOL')
    pymol.cmd.set_color("c1",[1.000, 0.675, 0.718])
    #pymol.cmd.set_color("c2",[0.8353,0.6039,0.7098])
    pymol.cmd.remove('not polymer')
    #pymol.cmd.util.cbc('all', first_color=7, quiet=1, legacy=0, _self=cmd)
    pymol.cmd.select('inter_1',f'{pdb_name} and chain {ChA} and byres(chain {ChB}) around {dist}')
    pymol.cmd.select('inter_2', f'{pdb_name} and chain {ChB} and byres(chain {ChA}) around {dist}')
    pymol.cmd.util.cba('c1',f'chain {ChA}')
    pymol.cmd.util.cba('slate', f'chain {ChB}')
    else_part  = pymol.cmd.select('else',f'{pdb_name} and  not (chain {ChA}+{ChB})')
    if else_part > 0:
        pymol.cmd.util.cba('lime', 'else')

    chA_D = f'inter_1 and  e. o+n and (neighbor e. h)'
    chB_A =  f'inter_2 and  e. o+n '
    DA_Hbond('HDA',chA_D,chB_A)
    chB_D = f'inter_2 and  e. o+n and (neighbor e. h)'
    chA_A = f'inter_1 and  e. o+n '
    DA_Hbond('HAD',chB_D,chA_A)
    aro_ring_name = {'HIS':'CG+ND1+CE1+NE2+CD2',
                     'TRP':'CH2+CZ3+CE3+CZ2+CE2+CD2',
                     'TYR':'CG+CD1+CE1+CZ+CE2+CD2',
                     'PHE':'CG+CD1+CE1+CZ+CE2+CD2'}
    pi_pi(aro_ring_name)
    Salt_bridge(ChA,ChB)
    pymol.cmd.do('set cartoon_transparency, 0.6')
    pymol.cmd.do('set dash_gap, 0.2')
    pymol.cmd.do('set dash_round_ends, 0')
    pymol.cmd.do('set label_size, 20')
    pymol.cmd.do('set dash_color, cyan,pi*')
    pymol.cmd.do('set dash_color, violet,SB*')
    pymol.cmd.do('set dash_radius, 0.1')
    pymol.cmd.do('set stick_radius, 0.15')
    pymol.cmd.delete('angle')
    pymol.cmd.disable('angle')
    pymol.cmd.set('bg_rgb', 'white')
    pymol.cmd.remove("all & hydro & not nbr. (don.|acc.)")
    pymol.cmd.do('set ray_trace_mode, 1')
    #cmd.do('set label_bg_color, white')
    #cmd.do('set label_bg_transparency, 0.7')
    pymol.cmd.do('set label_connector, on')
    pymol.cmd.do(f'set label_color, black ,chain {ChA}')
#    cmd.do('set stick_color, c1, chain A,elemC')
    pymol.cmd.do(f'set label_color, black ,chain {ChB}')
#    cmd.do('set stick_color, slate, chain B')
    pymol.cmd.do('set cartoon_loop_radius, 0.3')
    pymol.cmd.do('set ray_trace_mode, 1')
    pymol.cmd.do('set ray_shadows,0')
    pymol.cmd.do('set specular, 0')
    pymol.cmd.do('space cmyk')
    pymol.cmd.do('set ray_trace_color, [0,0,0]')
    pymol.cmd.do('set stick_h_scale,1')
    #set ray_shadow, off
    #ray_trace_mode set to 3
    #cmd.do('set ray_shadow,on')
    #cmd.save(f'{pdb_name}_interaction.pse')
pymol.cmd.extend("ppi", ppi)


def interfaceResidues(cmpx = '*', cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    """
    interfaceResidues -- finds 'interface' residues between two chains in a complex.

    PARAMS
        cmpx
            The complex containing cA and cB ——GHZ

        cA
            The first chain in which we search for residues at an interface
            with cB  ——GHZ

        cB
            The second chain in which we search for residues at an interface
            with cA  ——GHZ

        cutoff
            The difference in area OVER which residues are considered
            interface residues.  Residues whose dASA from the complex to
            a single chain is greater than this cutoff are kept.  Zero
            keeps all residues.

        selName
            The name of the selection to return.

    RETURNS
        * A selection of interface residues is created and named
            depending on what you passed into selName
        * An array of values is returned where each value is:
            ( modelName, residueNumber, dASA )

    NOTES
        If you have two chains that are not from the same PDB that you want
        to complex together, use the create command like:
            create myComplex, pdb1WithChainA or pdb2withChainX
        then pass myComplex to this script like:
            interfaceResidues myComlpex, c. A, c. X

        This script calculates the area of the complex as a whole.  Then,
        it separates the two chains that you pass in through the arguments
        cA and cB, alone.  Once it has this, it calculates the difference
        and any residues ABOVE the cutoff are called interface residues.

    AUTHOR:
        Jason Vertrees, 2009.
    """
    # Save user's settings, before setting dot_solvent
    oldDS = pymol.cmd.get("dot_solvent")
    pymol.cmd.set("dot_solvent", 1)

    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName + "1"
    chA, chB = "chA", "chB"

    # operate on a new object & turn off the original
    pymol.cmd.create(tempC, cmpx)
    pymol.cmd.disable(cmpx)

    # remove cruft and inrrelevant chains
    pymol.cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))

    # get the area of the complete complex
    pymol.cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    pymol.cmd.alter(tempC, 'q=b')

    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    pymol.cmd.extract(chA, tempC + " and (" + cA + ")")
    pymol.cmd.extract(chB, tempC + " and (" + cB + ")")
    pymol.cmd.get_area(chA, load_b=1)
    pymol.cmd.get_area(chB, load_b=1)

    # update the chain-only objects w/the difference
    pymol.cmd.alter("%s or %s" % (chA, chB), "b=b-q")

    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    pymol.cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

    pymol.cmd.enable(cmpx)
    pymol.cmd.select(selName1, 'none')
    for (model, resi, diff) in stored.r:
        key = resi + "-" + model
        if abs(diff) >= float(cutoff):
            if key in seen:
                continue
            else:
                seen.append(key)
            rVal.append((model, resi, diff))
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            pymol.cmd.select(selName1, selName1 + " or (%s and i. %s)" % (model, resi))

    # this is how you transfer a selection to another object.
    pymol.cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    pymol.cmd.delete(selName1)
    pymol.cmd.delete(chA)
    pymol.cmd.delete(chB)
    pymol.cmd.delete(tempC)
    # show the selection
    pymol.cmd.enable(selName)
    #cmd.show('sticks',selection=selName)#####################################
    pymol.cmd.color('blue', selection=selName)

    pymol.cmd.util.cbc()
    pymol.cmd.util.cbak(selection=selName)
    pymol.cmd.distance(name = 'polar_contact', selection1= 'c. A', selection2='c. B',mode = 2)
    pymol.cmd.color('yellow',selection= 'polar_contact')
    pymol.cmd.distance(name = 'pi-pi', selection1= 'c. A', selection2='c. B',mode = 1)
    pymol.cmd.color('pink',selection= 'pi-pi')

    # reset users settings
    pymol.cmd.set("dot_solvent", oldDS)
    #cmd.save(out_put_name)
    #return rVal

pymol.cmd.extend("interfaceResidues", interfaceResidues)

def cif_to_pdb(input_folder, output_folder='',):
    """
    将指定文件夹下所有.cif文件转换为.pdb格式
    :param input_folder: 输入文件夹路径
    :param output_folder: 输出文件夹路径（默认同输入文件夹）
    """
    # 设置输出路径
    if output_folder == '':
        output_folder = input_folder
    elif output_folder != '':
        os.makedirs(output_folder, exist_ok=True)

    # 初始化解析器
    parser = MMCIFParser()
    io = PDBIO()

    # 遍历文件夹
    converted = 0
    errors = []

    for filename in os.listdir(input_folder):
        if not filename.lower().endswith(".cif"):
            continue

        # 构造文件路径
        input_path = os.path.join(input_folder, filename)
        base_name = os.path.splitext(filename)[0]
        output_path = os.path.join(output_folder, f"{base_name}.pdb")

        try:
            # 解析CIF文件
            structure = parser.get_structure(base_name, input_path)

            # 保存为PDB
            io.set_structure(structure)
            io.save(output_path)

            print(f"转换成功: {filename} -> {base_name}.pdb")
            converted += 1

        except Exception as e:
            errors.append((filename, str(e)))
            print(f"转换失败: {filename} - {str(e)}")

    # 输出统计信息
    print(f"\n转换完成！成功: {converted}, 失败: {len(errors)}")
    if errors:
        print("\n错误详情:")
        for filename, error in errors:
            print(f"- {filename}: {error}")

def cif_to_pdb_v2(input_folder, output_folder=None):
    """
    将输入文件夹中的所有 CIF 文件转换为 PDB 格式

    参数:
        input_folder: 包含 CIF 文件的输入文件夹路径
        output_folder: 输出 PDB 文件的文件夹路径 (默认: 输入文件夹下创建 'pdb_output' 文件夹)
    """
    # 设置默认输出文件夹
    if output_folder is None:
        #output_folder = os.path.join(input_folder, 'pdb_output')
        output_folder = input_folder

    # 创建输出文件夹（如果不存在）
    os.makedirs(output_folder, exist_ok=True)

    # 获取所有 CIF 文件
    cif_files = [f for f in os.listdir(input_folder)
                 if f.lower().endswith(('.cif', '.cif.gz'))]

    if not cif_files:
        print(f"在 {input_folder} 中没有找到 CIF 文件")
        return

    print(f"找到 {len(cif_files)} 个 CIF 文件，开始转换...")

    # 初始化解析器和写入器
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    success_count = 0
    failed_files = []

    for cif_file in cif_files:
        input_path = os.path.join(input_folder, cif_file)
        pdb_file = os.path.splitext(cif_file)[0] + '.pdb'
        if cif_file.lower().endswith('.gz'):
            pdb_file = os.path.splitext(pdb_file)[0] + '.pdb'  # 处理 .cif.gz 文件

        output_path = os.path.join(output_folder, pdb_file)

        try:
            # 读取结构
            structure = parser.get_structure(os.path.splitext(cif_file)[0], input_path)

            # 保存为 PDB
            io.set_structure(structure)
            io.save(output_path)

            print(f"转换成功: {cif_file} -> {pdb_file}")
            success_count += 1

        except Exception as e:
            error_msg = f"转换失败: {cif_file} - {str(e)}"
            print(error_msg)
            failed_files.append((cif_file, str(e)))
            traceback.print_exc()  # 打印详细错误信息用于调试

    # 输出总结报告
    print("\n" + "=" * 50)
    print(f"转换完成! 成功: {success_count}/{len(cif_files)}")

    if failed_files:
        print("\n失败文件列表:")
        for file, error in failed_files:
            print(f"- {file}: {error}")

    print(f"所有 PDB 文件已保存至: {os.path.abspath(output_folder)}")



def simplify_chain_ids(structure):
    """
    简化结构中的链ID（如将"A-1"简化为"A"）

    参数:
        structure: Bio.PDB 结构对象

    返回:
        修改后的结构对象
    """
    # 用于记录原始链ID到简化链ID的映射
    chain_id_map = {}

    # 遍历所有模型
    for model in structure:
        # 用于跟踪已使用的简化链ID，避免冲突
        used_chain_ids = set()

        # 第一次遍历：收集所有需要简化的链ID
        for chain in model:
            original_id = chain.id
            # 匹配"A-1"、"B-1"等格式
            match = re.match(r'^([A-Za-z])-\d+$', original_id)
            if match:
                simplified_id = match.group(1)
                # 记录映射关系
                if original_id not in chain_id_map:
                    chain_id_map[original_id] = simplified_id
            else:
                # 不符合模式的链保持原样
                if original_id not in chain_id_map:
                    chain_id_map[original_id] = original_id

        # 第二次遍历：实际修改链ID
        for chain in model:
            original_id = chain.id
            simplified_id = chain_id_map[original_id]

            # 检查是否有冲突
            if simplified_id in used_chain_ids:
                # 如果冲突，使用原始ID
                new_id = original_id
            else:
                new_id = simplified_id
                used_chain_ids.add(simplified_id)

            # 修改链ID
            chain.id = new_id

    return structure

def convert_cif_to_pdb_with_simplified_chains(cif_file, pdb_file):
    """
    将CIF文件转换为PDB文件，并简化链ID

    参数:
        cif_file: 输入的CIF文件路径
        pdb_file: 输出的PDB文件路径
    """
    try:
        # 创建MMCIF解析器
        parser = MMCIFParser(QUIET=True)

        # 解析结构
        structure = parser.get_structure(os.path.basename(cif_file), cif_file)

        # 简化链ID
        structure = simplify_chain_ids(structure)

        # 保存为PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)

        print(f"转换成功: {os.path.basename(cif_file)} -> {os.path.basename(pdb_file)}")
        return True

    except Exception as e:
        print(f"处理文件时出错: {cif_file}")
        print(f"错误信息: {str(e)}")
        traceback.print_exc()
        return False

def remove_dash_in_chain_name_in_folder(input_folder, output_folder=None):
    """
    处理文件夹中的所有CIF文件

    参数:
        input_folder: 输入文件夹路径
        output_folder: 输出文件夹路径（默认在输入文件夹下创建'pdb_output'文件夹）
    """
    # 设置默认输出文件夹
    if output_folder is None:
        #output_folder = os.path.join(input_folder, 'pdb_output')
        output_folder = input_folder

    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)

    # 获取所有CIF文件
    cif_files = [f for f in os.listdir(input_folder)
                 if f.lower().endswith(('.cif', '.cif.gz'))]

    if not cif_files:
        print(f"在 {input_folder} 中没有找到CIF文件")
        return

    print(f"找到 {len(cif_files)} 个CIF文件，开始转换并简化链ID...")
    print("=" * 50)

    success_count = 0
    failed_files = []

    for cif_file in cif_files:
        input_path = os.path.join(input_folder, cif_file)

        # 生成输出文件名（保留原始文件名，扩展名改为.pdb）
        base_name = os.path.splitext(cif_file)[0]
        if base_name.lower().endswith('.cif'):
            base_name = os.path.splitext(base_name)[0]
        pdb_file = base_name + '.pdb'
        output_path = os.path.join(output_folder, pdb_file)

        # 处理文件
        if convert_cif_to_pdb_with_simplified_chains(input_path, output_path):
            success_count += 1
        else:
            failed_files.append(cif_file)

    # 输出总结报告
    print("\n" + "=" * 50)
    print(f"处理完成! 成功: {success_count}/{len(cif_files)}")

    if failed_files:
        print("\n失败文件列表:")
        for file in failed_files:
            print(f"- {file}")

    print(f"所有PDB文件已保存至: {os.path.abspath(output_folder)}")


def there_is_pdb_in_this_folder(input_folder):
    thereis_pdb = 0
    for each_file in os.listdir(input_folder):
        if each_file.endswith('.pdb'):
            thereis_pdb = 1

    return thereis_pdb
#################################################################
three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'UNK': 'X'
    }

one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
        'U': 'SEC', 'O': 'PYL', 'X': 'UNK'
    }
###############################################

def get_the_sequence_from_single_pdb(input_pdb = r'D:\my_experi\structure_prediction\20250515_Nanobody_CMV\CMV_CP_core',
                                     chain_id = 'A',silent = 0):
    '''
    这个脚本是将一个PDF文件，中的一个chain，转化成为序列
    需要同时指定PDF文件的路径和链名称。
    :param input_pdb:
    :param chain_id:
    :return:
    '''

    if not input_pdb.endswith('.pdb'):
        print('输入文件似乎并不是一个pdb格式的文件。')
        return 1

    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb)

    sequence_list = []
    sequence = ''
    residue_numbers = []

    for chain in structure.get_chains():
        if chain.get_id() == chain_id:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # 只考虑标准氨基酸
                    try:
                        aa = three_to_one[residue.get_resname()]
                        sequence_list.append(aa)
                        sequence += aa
                        residue_numbers.append(residue.get_id()[1])
                    except KeyError:
                        pass  # 跳过非标准氨基酸
            break
    if silent == 0:
        print('sequence_list',sequence_list)
        print('sequence', sequence)
        print('residue_numbers',residue_numbers)
    #return ''.join(sequence_list), residue_numbers
    return sequence,residue_numbers

def get_aa_sequence_from_pdb_folder(input_pdb_folder = '',
                                    chain_id = 'all',silent = 0):
    '''
    这个程序是将一个文件夹下面所有的pdb文件，全部提取氨基酸序列，并输出成为fasta文件和csv文件。
    :param input_pdb_folder:
    :param chain_id: 指定要提取序列的chain_id，默认为all，即提取所有的序列。
    :param silent:
    :return:
    '''
    pass


def get_sequence_from_pdb(pdb_path = r'D:\my_experi\structure_prediction\20250515_Nanobody_CMV\CMV_CP_core',
                          chain_required = ''):
    '''
    这个程序是用pymol程序包写的，其实不是那么的友好。
    #这个脚本是将一个目录下所有的PDF文件，转化成为fasta格式的序列信息，并输出成为文件
    # 定义三字母到单字母的映射
    '''

    cif_to_pdb(pdb_path)
    os.chdir(pdb_path)
    #确定一个输出文件
    out_file = open('out_fasta.fa',mode='w')
    for each_file in os.listdir(pdb_path):
        if each_file.endswith('pdb'):
            print(each_file)
            #说明这是一个PDF文件。
            cmd.reinitialize()
            cmd.load(each_file)
            object_name = cmd.get_object_list()[0]
            # 自动获取加载的对象名称
            #print('object_name',object_name)
            #cmd.save('abc.pdb')

            chains = pymol.cmd.get_chains(object_name)
            # 获取所有链名
            #all_fasta_str = []

            # 遍历每条链并输出其氨基酸序列
            for chain in chains:
                if chain_required in chain or chain in chain_required:
                # 初始化存储氨基酸序列的列表
                    all_fasta_seq = ''
                    out_file.write('>'+each_file.rsplit('.')[0]+'_'+object_name+'_'+chain+'\n')
                    fasta_list = []


                    # 使用 iterate 获取每个残基的信息（基于 Cα 原子）
                    cmd.iterate(f"{object_name} and chain {chain} and name CA",
                                'fasta_list.append((resi, resn))', space=locals())

                    # 构建 FASTA 格式的字符串，使用单字母码表示氨基酸
                    for resi, resn in fasta_list:
                        one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
                        #all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
                        all_fasta_seq += one_letter
                    out_file.write(all_fasta_seq+'\n')

    out_file.close()
    return all_fasta_seq

#另一个从pdb文件中读取protein序列的脚本。但是这里，只适用于只包含nanobody和antigen的pdb
#这个pdb文件只包含两条链
#这个程序会输出一个fasta文件，并适合于用于igGM分析的输入文件
def get_sequence_from_nanobody_antigen_pdb(pdb_path = r"D:\my_experi\structure_prediction\20250521_Nanobody_ToBRFV_CPv2\fold_tobrfv_cpv2_nb29\fold_tobrfv_cpv2_nb29_model_0.pdb",
                                           change_chain_name = 1, output_fasta = False, return_antigen_nanobody_sequence = True,
                                           antigen_name = 'T',nanobody_name = 'H'):
    '''
    这个脚本是将一个pdb文件
    从nanobody_antigen_pdb文件中获取序列信息，并保存为fasta格式
    change_chain_name 是否对chain的名称进行改名
    将nanobody的chain name改为H ，antigen的chain name改为A
    :param pdb_path:
    :return:
    '''
    if '\\' in pdb_path:
        pdb_path_folder= pdb_path.rsplit('\\',maxsplit=1)[0]
        pdb_file_name = pdb_path.rsplit('\\',maxsplit=1)[1]
    else:
        pdb_file_name = pdb_path
        pdb_path_folder = '.'
    os.chdir(pdb_path_folder)
    if output_fasta:
        out_file = open(pdb_file_name.rsplit('.',maxsplit=1)[0]+'.fasta', mode='w')

    #print(pdb_path_folder)
    #print(pdb_file_name)
    #print(pdb_path)

    cmd.reinitialize()
    cmd.load(pdb_path)
    object_name = cmd.get_object_list()[0]
    # 自动获取加载的对象名称
    # print('object_name',object_name)
    # cmd.save('abc.pdb')

    chains = pymol.cmd.get_chains(object_name)
    # 获取所有链名
    # all_fasta_str = []

    # 遍历每条链并输出其氨基酸序列
    for chain in chains:
        #print(chain,chains)
        # 初始化存储氨基酸序列的列表
        all_fasta_seq = ''
        # 使用 iterate 获取每个残基的信息（基于 Cα 原子）
        fasta_list = []
        pymol.cmd.iterate(f"{object_name} and chain {chain} and name CA",
                    'fasta_list.append((resi, resn))', space=locals())

        for resi, resn in fasta_list:
            one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
            # all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
            all_fasta_seq += one_letter
            # 构建 FASTA 格式的字符串，使用单字母码表示氨基酸

        #接下来判断这一条链是否是nanobody or antigen
        #print('all_fasta_seq',all_fasta_seq)
        my_nanobody= My_nanobody(all_fasta_seq,stringent=0) #用相对不严格的参数。
        if my_nanobody.is_nanobody() == True:
            nanobody_sequence = all_fasta_seq
            if output_fasta:
                out_file.write('>'+nanobody_name + '\n')
                out_file.write(all_fasta_seq + '\n')
            if change_chain_name == 1:
                pymol.cmd.alter("chain "+chain,f"chain = '{nanobody_name}'")

        if my_nanobody.is_nanobody() == False:
            antigen_sequence = all_fasta_seq
            #那么这一条链则被视为antigen
            if output_fasta:
                out_file.write('>' +antigen_name + '\n')
                out_file.write(all_fasta_seq + '\n')
            if change_chain_name == 1:
                pymol.cmd.alter("chain "+chain,f"chain = '{antigen_name}'")

    if change_chain_name == 1:
        os.remove(pdb_path)
        pymol.cmd.save(format='pdb',filename=pdb_path)
        #print('pdb_path',pdb_path)
    if output_fasta:
        out_file.close()

    if return_antigen_nanobody_sequence:
        #print('antigen_sequence',antigen_sequence)
        #print('nanobody_sequence',nanobody_sequence)
        return antigen_sequence, nanobody_sequence
    pass

def get_sequence_from_antigen_binder_pdb(pdb_path = r"D:\my_experi\structure_prediction\20250521_Nanobody_ToBRFV_CPv2\fold_tobrfv_cpv2_nb29\fold_tobrfv_cpv2_nb29_model_0.pdb",
                                           change_chain_name = 1, output_fasta = False, return_antigen_binder_sequence = True,
                                           antigen_name = 'T',binder_name = 'B',
                                         antigen_collection_fasta = r'/public-supool/home/gaolab/haochengzhu/Effector_pdb/pathogen_effectors.fasta'):
    '''
    这个脚本是将一个pdb文件
    从antigen_binder_pdb文件中获取序列信息，并保存为fasta格式
    change_chain_name 是否对chain的名称进行改名
    将binder的chain name改为B ，antigen的chain name改为T
    '''
    if '\\' in pdb_path:
        pdb_path_folder= pdb_path.rsplit('\\',maxsplit=1)[0]
        pdb_file_name = pdb_path.rsplit('\\',maxsplit=1)[1]
    else:
        pdb_file_name = pdb_path
        pdb_path_folder = '.'
    os.chdir(pdb_path_folder)
    print('pdb_path',pdb_path)



    #print(pdb_path_folder)
    #print(pdb_file_name)
    #print(pdb_path)

    pymol.cmd.reinitialize()
    pymol.cmd.load(pdb_path)
    object_name = pymol.cmd.get_object_list()[0]
    # 自动获取加载的对象名称
    # print('object_name',object_name)
    # pymol.cmd.save('abc.pdb')

    chains = pymol.cmd.get_chains(object_name)
    # 获取所有链名
    # all_fasta_str = []

    # 遍历每条链并输出其氨基酸序列
    chain_sequences = []
    for chain in chains:
        #print(chain,chains)
        # 初始化存储氨基酸序列的列表
        all_fasta_seq = ''
        # 使用 iterate 获取每个残基的信息（基于 Cα 原子）
        fasta_list = []
        pymol.cmd.iterate(f"{object_name} and chain {chain} and name CA",
                    'fasta_list.append((resi, resn))', space=locals())

        for resi, resn in fasta_list:
            one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
            # all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
            all_fasta_seq += one_letter
            # 构建 FASTA 格式的字符串，使用单字母码表示氨基酸

        #接下来判断这一条链是否是nanobody or binder
        #判断的依据是，那个短，哪个就是binder
        chain_sequences.append(all_fasta_seq)

    #print('chain_sequences',chain_sequences)
    #读取chain_sequences
    find_antigen_in_collection = 0
    for each_fasta in SeqIO.parse(antigen_collection_fasta, "fasta"):
        if str(each_fasta.seq) in chain_sequences[0] or chain_sequences[0] in str(each_fasta.seq):
            binder_sequence = chain_sequences[1]
            antigen_sequence = chain_sequences[0]
            find_antigen_in_collection = 1
            print(f'find_antigen_in_fasta_collection,:{str(each_fasta.id)}:{str(each_fasta.seq)}')
            break
        if str(each_fasta.seq) in chain_sequences[1] or chain_sequences[1] in str(each_fasta.seq):
            binder_sequence = chain_sequences[0]
            antigen_sequence = chain_sequences[1]
            find_antigen_in_collection = 1
            print(f'find_antigen_in_fasta_collection,:{str(each_fasta.id)}:{str(each_fasta.seq)}')
            break

    if find_antigen_in_collection == 0:
        if len(chain_sequences[1]) > len(chain_sequences[0]):
            binder_sequence = chain_sequences[0]
            antigen_sequence = chain_sequences[1]

        elif len(chain_sequences[1]) == len(chain_sequences[0]):
            binder_sequence = chain_sequences[0]
            antigen_sequence = chain_sequences[1]

        elif len(chain_sequences[1]) < len(chain_sequences[0]):
            binder_sequence = chain_sequences[1]
            antigen_sequence = chain_sequences[0]
    print('binder_sequence',binder_sequence)
    print('antigen_sequence', antigen_sequence)

    if output_fasta:
        out_file = open(pdb_file_name.rsplit('.', maxsplit=1)[0] + '.fasta', mode='w')
        out_file.write('>' + antigen_name + '\n')
        out_file.write(antigen_sequence + '\n')
        out_file.write('>' + binder_name + '\n')
        out_file.write(binder_sequence + '\n')
        out_file.close()
    if change_chain_name == 1:
        used = 0
        for chain in chains:
            #print('chain',chain)

            all_fasta_seq = ''
            fasta_list = []
            pymol.cmd.iterate(f"{object_name} and chain {chain} and name CA",
                        'fasta_list.append((resi, resn))', space=locals())
            for resi, resn in fasta_list:
                one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
                # all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
                all_fasta_seq += one_letter

            #print('all_fasta_seq',all_fasta_seq)
            if all_fasta_seq == antigen_sequence and used == 0:
                antigen_chain = chain
                #pymol.cmd.alter("chain " + chain, f"chain = '{antigen_name}'")
                used = 1

            elif all_fasta_seq == binder_sequence:
                binder_chain = chain
                #pymol.cmd.alter("chain " + chain, f"chain = '{binder_name}'")

    if change_chain_name == 1:
        pymol.cmd.alter("chain " + antigen_chain, f"chain = '{antigen_name}'")
        pymol.cmd.alter("chain " + binder_chain, f"chain = '{binder_name}'")
        os.remove(pdb_path)
        pymol.cmd.save(format='pdb',filename=pdb_path)
        #print('pdb_path',pdb_path)


    if return_antigen_binder_sequence:
        #print('antigen_sequence',antigen_sequence)
        #print('nanobody_sequence',nanobody_sequence)
        return antigen_sequence, binder_sequence
    pass

def get_sequence_from_binder_dimer_pdb(pdb_path = r"D:\my_experi\structure_prediction\20250521_Nanobody_ToBRFV_CPv2\fold_tobrfv_cpv2_nb29\fold_tobrfv_cpv2_nb29_model_0.pdb",
                                           change_chain_name = 1, output_fasta = False, return_antigen_binder_sequence = True,
                                           binder_name1 = 'X',binder_name2 = 'Y'):
    '''
    这个脚本是将一个pdb文件
    从antigen_binder_pdb文件中获取序列信息，并保存为fasta格式
    change_chain_name 是否对chain的名称进行改名
    将binder的chain name改为B ，antigen的chain name改为T
    '''
    if '\\' in pdb_path:
        pdb_path_folder= pdb_path.rsplit('\\',maxsplit=1)[0]
        pdb_file_name = pdb_path.rsplit('\\',maxsplit=1)[1]
    else:
        pdb_file_name = pdb_path
        pdb_path_folder = '.'
    os.chdir(pdb_path_folder)

    #print(pdb_path_folder)
    #print(pdb_file_name)
    #print(pdb_path)

    pymol.cmd.reinitialize()
    pymol.cmd.load(pdb_path)
    object_name = pymol.cmd.get_object_list()[0]
    # 自动获取加载的对象名称
    # print('object_name',object_name)
    # pymol.cmd.save('abc.pdb')

    chains = pymol.cmd.get_chains(object_name)
    # 获取所有链名
    # all_fasta_str = []

    # 遍历每条链并输出其氨基酸序列
    chain_sequences = []
    for chain in chains:
        #print(chain,chains)
        # 初始化存储氨基酸序列的列表
        all_fasta_seq = ''
        # 使用 iterate 获取每个残基的信息（基于 Cα 原子）
        fasta_list = []
        pymol.cmd.iterate(f"{object_name} and chain {chain} and name CA",
                    'fasta_list.append((resi, resn))', space=locals())

        for resi, resn in fasta_list:
            one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
            # all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
            all_fasta_seq += one_letter
            # 构建 FASTA 格式的字符串，使用单字母码表示氨基酸

        #接下来判断这一条链是否是nanobody or binder
        #判断的依据是，那个短，哪个就是binder
        chain_sequences.append(all_fasta_seq)


    if len(chain_sequences[1]) == len(chain_sequences[0]):
        binder_sequence1 = chain_sequences[0]
        binder_sequence2 = chain_sequences[1]

    if output_fasta:
        out_file = open(pdb_file_name.rsplit('.', maxsplit=1)[0] + '.fasta', mode='w')
        out_file.write('>' + binder_name1 + '\n')
        out_file.write(binder_sequence1 + '\n')
        out_file.write('>' + binder_name2 + '\n')
        out_file.write(binder_sequence2 + '\n')
        out_file.close()
    if change_chain_name == 1:
        used = 0
        for chain in chains:

            all_fasta_seq = ''
            fasta_list = []
            pymol.cmd.iterate(f"{object_name} and chain {chain} and name CA",
                        'fasta_list.append((resi, resn))', space=locals())
            for resi, resn in fasta_list:
                one_letter = three_to_one.get(resn, 'X')  # 如果没有对应的单字母缩写，使用 'X'
                # all_fasta_str.append(f"{chain} {resi} {resn} {one_letter}")
                all_fasta_seq += one_letter

            if all_fasta_seq == binder_sequence1 and used == 0:
                pymol.cmd.alter("chain " + chain, f"chain = '{binder_name1}'")
                used = 1
            if all_fasta_seq == binder_sequence2:
                pymol.cmd.alter("chain " + chain, f"chain = '{binder_name2}'")

    if change_chain_name == 1:
        os.remove(pdb_path)
        pymol.cmd.save(format='pdb',filename=pdb_path)
        #print('pdb_path',pdb_path)


    if return_antigen_binder_sequence:
        #print('antigen_sequence',antigen_sequence)
        #print('nanobody_sequence',nanobody_sequence)
        return binder_sequence1, binder_sequence2
    pass



def standardize_100w_nanobody(nanobody_100w_input = r"D:\my_script\Nanobody\nano_collection\INDI_100w_nanobody.fasta",stringent = 0,output_dir = ''):
    '''
    :param nanobody_100w_input: 输入文件
    :param stringent:
    :param output_dir:
    :return:
    '''
    if output_dir == '':
        output_dir = nanobody_100w_input.rsplit(sep = '\\',maxsplit=1)[0]
    os.chdir(output_dir)
    out_100w_fasta = open('100w_nanobody_standardized.fasta',mode = 'w')
    i = 0
    for each_fasta in SeqIO.parse(nanobody_100w_input, "fasta"):
        #print(record.id)
        #print(record)
        nanobody_origin_seq = str(each_fasta.seq)
        my_nanobody = My_nanobody(sequence=nanobody_origin_seq,stringent=stringent,silent=1)
        if my_nanobody.is_nanobody():
            #print(my_nanobody.standized_nb())
            out_100w_fasta.write('>'+each_fasta.id+'\n')
            out_100w_fasta.write(my_nanobody.standized_nb() + '\n')

        i+=1
        if i%5000 == 0:
            print('已完成'+str(i)+'个nanobody的收集')

    out_100w_fasta.close()



        #然后将标准化号的nanobody放置在一个新的fasta文件里




    pass
def get_sequence_from_fasta(antigen_name = 'BSCTV_CP_core',antigen_coll_fa = r"D:\my_experi\structure_prediction\pathogen_effectors.fasta"):
    #这个脚本是从一个目的fasta文件中获取相对应名称的序列。
    for each_fasta in SeqIO.parse(antigen_coll_fa, "fasta"):
        #print(record.id)
        #print(record)
        if each_fasta.id == antigen_name:
            print(str(each_fasta.seq))
            return str(each_fasta.seq)

def fasta_2_dataframe(input_fasta_file = r"D:\my_experi\structure_prediction\pathogen_effectors.fasta",
                      required_id_list:list = ['all']):
    '''
    这个脚本是将输入的fasta文件转变成为datafame格式 并输出
    :param fasta_file:
    required_id_list : 所需要的所有fasta 默认是['all'],意味着所有都要
    如果指定了，则只挑选所需要fasta序列进行dataframe 化
    :return: datafame
    '''
    out_df = pd.DataFrame()
    if required_id_list == ['all']:
        for each_fasta in SeqIO.parse(input_fasta_file, "fasta"):
            df2append = pd.DataFrame()
            df2append['id'] = str(each_fasta.id).split()
            df2append['sequence'] = str(each_fasta.seq).split()
            out_df = out_df._append(df2append)
    else:
        for each_fasta in SeqIO.parse(input_fasta_file, "fasta"):
            df2append = pd.DataFrame()
            if str(each_fasta.id) in required_id_list:

                df2append['id'] = str(each_fasta.id).split()
                df2append['sequence'] = str(each_fasta.seq).split()
                out_df = out_df._append(df2append)

    out_df.index = np.arange(len(out_df))
    #print('out_df',out_df)
    return out_df


def data_prepare_4_deepnano(antigen_name = 'BSCTV_CP_core',antigen_name_list = ['ToLCNDV_CP_core','SSCTV_CP_core'],
                            nanobody_collection_fasta = r"D:\my_script\Nanobody\nano_collection\nanobody_collection1.fasta",
                            base_path = r'D:\my_script\Nanobody\DeepNano_qinghua\My_cases'):
    '''
    #这个脚本的目的是用来，生成deepnano脚本的输入序列的
    :param antigen_name: antigen_name，需要在pathogen_effectors.fasta这个文件中
    :param antigen_name_list: list of antigen_name，需要在pathogen_effectors.fasta这个文件中
    antigen_name与antigen_name_list不兼容。当都有赋值的时候，选取antigen_name_list
    :param nanobody_collection_fasta:
    :param base_path:
    :return:
    '''

    '''
    #这个脚本的目的是用来，生成deepnano脚本的输入序列的
    :param antigen_name: antigen_name，需要在pathogen_effectors.fasta这个文件中
    :param nanobody_collection_fasta: 待虚拟筛选的nanobody collection 需要是fasta格式
    :return:以antigen_name命名的一个文件夹
    '''

    antigen_coll_fa = r"D:\my_experi\structure_prediction\pathogen_effectors.fasta"



    if antigen_name_list == ['']:
        #明确输出文件的文件名称。
        output_fasta_name = antigen_name+'_nanobody_collection.fasta'
        output_tsv_name = antigen_name+'_nanobody_collection.tsv'

        output_path = 'deepnano_input_' + antigen_name + '_' + get_current_date_hczhu()
    else:
        antigen_name = antigen_name_list[0]
        output_fasta_name = antigen_name+'_etc_nanobody_collection.fasta'
        output_tsv_name = antigen_name+'_etc_nanobody_collection.tsv'
        output_path = 'deepnano_input_' + antigen_name +'_etc_' + get_current_date_hczhu()

    #调整一下输出路径
    os.chdir(base_path)
    if output_path not in os.listdir():
        os.mkdir(output_path)
        #不存在这个路径，需要在这里创制一下。
    os.chdir(output_path)

    #创制输入和输出问题。
    if 'input' not in os.listdir():
        os.mkdir('input')
    if 'output' not in os.listdir():
        os.mkdir('output')
    os.chdir('input')


    #然后将这一条antigen 序列放进nanobody_collection_fasta之中
    nanobody_collection_fasta_file = open(nanobody_collection_fasta).readlines()
    output_fasta_write = open(output_fasta_name, 'w')
    output_tsv_write = open(output_tsv_name, 'w')
    if antigen_name_list == ['']:
        antigen_name_list = antigen_name.split()
        #将antigen_name也转变成为一个列表
        # 首先获取antigen_name的对应序列
    for each_antigen in antigen_name_list:
        antigen_seq = get_sequence_from_fasta(antigen_name=each_antigen, antigen_coll_fa=antigen_coll_fa)
        output_fasta_write.write('>'+each_antigen+'\n')
        output_fasta_write.write(antigen_seq + '\n')
    for each_line in nanobody_collection_fasta_file:
        #把逗号替换为‘_’
        each_line = each_line.replace(',','_')
        output_fasta_write.write(each_line[:-1]+'\n')
        if each_line.startswith('>'):
            #这是名称
            nanobody_name = each_line[1:-1]
            for each_antigen in antigen_name_list:
                output_tsv_write.write(nanobody_name+'\t'+each_antigen+'\t0\n')

    output_fasta_write.close()
    output_tsv_write.close()

def data_analysis_from_deepnano_predict(deepnanp_output_file = r"D:\my_script\Nanobody\DeepNano_qinghua\input\deepnano_input_TYLCV_CP_etc\TYLCV_CP_etc_predictions2.csv",
                                        output_dir = '',
                                        fasta_file = '',analysis_porportion:float = 0.001):
    '''
    这个脚本是deepnano pipline产出的predictions.csv类型的文件进行分析并作图。
    :param deepnanp_output:
    :param output_dir:
    :return:
    '''
    if output_dir == '':
        output_dir = deepnanp_output_file.rsplit(sep = '\\',maxsplit=1)[0]
    os.chdir(output_dir)

    #
    #需要去input 文件夹下面寻找
    os.chdir('..\\input')
    if fasta_file == '':
        #在当前目录下寻找’collection.fasta‘文件
        find_the_fasta_file = 0
        for each_file in os.listdir('.'):
            if each_file.endswith('collection.fasta'):
                find_the_fasta_file = 1
                print(each_file)

                fasta_file = each_file
                #fasta_df = fasta_2_dataframe(fasta_file)
                #这一步比较耗时间
                #print(fasta_df)
        if find_the_fasta_file == 0:
            print('未找到fasta_file，请检查目录下文件或者额外指定。')

    #然后读取deepnanp_out
    deepnano_out_df = pd.read_csv(deepnanp_output_file)
    deepnano_out_df.index = np.arange(len(deepnano_out_df))
    #重命名一下列的名字。
    deepnano_out_df.columns = ['Nanobody_ID','Antigen_ID','Prediction']
    #根据Prediction的分值进行排序
    deepnano_out_df = deepnano_out_df.sort_values(by = 'Prediction',ascending= False)
    deepnano_out_df['rank'] = np.arange(len(deepnano_out_df))
    deepnano_out_df['relative_rank'] = deepnano_out_df['rank']/len(deepnano_out_df)

    #只需要分析前1%即可
    deepnano_out_df_need_seq = deepnano_out_df.loc[deepnano_out_df['relative_rank'] <= analysis_porportion,:]
    list_of_seq_needed = list(deepnano_out_df_need_seq['Nanobody_ID'])
    list_of_seq_needed += list(deepnano_out_df_need_seq['Antigen_ID'])
    list_of_seq_needed = list(set(list_of_seq_needed))
    #去重
    #print(list_of_seq_needed)
    fasta_df = fasta_2_dataframe(fasta_file,required_id_list=list_of_seq_needed)
    print(fasta_df)
    #只获取所需要的fasta的序列

    #然后将Antigen_ID与Nanobody_ID进行合并 方便后面复制
    for each_line in deepnano_out_df.index:
        deepnano_out_df.loc[each_line,'Antigen_Nanobody_ID'] = deepnano_out_df.loc[each_line,'Antigen_ID'] +'_' + deepnano_out_df.loc[each_line,'Nanobody_ID']


    #然后根据输入的fasta文件 补全antigen和nanobody的序列
    #利用index来补全 这样快一些
    fasta_df.index = fasta_df['id']
    deepnano_out_df.index = deepnano_out_df['Antigen_ID']
    deepnano_out_df['antigen_sequence'] = fasta_df['sequence']
    print('antigen_sequence acquired, done')
    deepnano_out_df.index = deepnano_out_df['Nanobody_ID']
    deepnano_out_df['nanobody_sequence'] = fasta_df['sequence']
    print('nanobody_sequence acquired, done')


    #然后再把index全都换回来
    fasta_df.index = np.arange(len(fasta_df))
    deepnano_out_df.index = np.arange(len(deepnano_out_df))



    #再次调整路径
    os.chdir(output_dir)
    #输出
    output_df_name = deepnanp_output_file.split('\\')[-1].rsplit('.',maxsplit=1)[0]+'_analysis.csv'
    deepnano_out_df.to_csv(output_df_name)


    paint = 1
    if paint == 1:
        #然后作图
        '''
        plot1 = \
        (
            plotnine.ggplot(data = deepnano_out_df, mapping = aes(x='rank', y='Prediction')) +
            plotnine.geom_point() +
            plotnine.theme_bw()

        )
        plotnine.ggsave(plot1,format='svg')
        '''
        color = deepnano_out_df['Prediction']
        plot2 = \
        (
            plotnine.ggplot(data = deepnano_out_df, mapping = aes(x='Prediction')) +
            plotnine.geom_histogram(fill = 'grey') +
            plotnine.theme_bw()

        )
        plotnine.ggsave(plot2, filename='hist1.png', format='png',width= 7, height= 5)
        plotnine.ggsave(plot2,filename= 'hist1.svg',format='svg',width= 7, height= 5)








    #首先进行排序


def fasta_merge(input_folder = r'D:\my_script\Nanobody\IgGM\My_cases\tobrfv_cpv2_20250614\output',key = 'default',
                output_fasta_name = 'merged_fasta.fasta',
                output_df_name = 'merged_fasta.csv'):
    #return output_df
    #这个脚本是将目录下所有的fasta文件进行汇总，汇总到一个fasta文件里面。。
    #key 只汇总特定ID的fasta序列。
    os.chdir(input_folder)
    output_file = open(output_fasta_name,'w')
    output_df = pd.DataFrame()
    for each_file in os.listdir():
        if each_file.endswith('fasta'):
            #这是一个目的fast文件。
            for each_fasta in SeqIO.parse(each_file, "fasta"):
                fasta_id = str(each_fasta.id)
                if key == fasta_id or key == 'default':
                    #此时是默认值，对所有世界都进行分析。
                    output_file.write('>'+str(each_fasta.id)+'\n')
                    output_file.write(str(each_fasta.seq)+'\n')
                    df2append = pd.DataFrame()
                    df2append['id'] = str(each_fasta.id).split()
                    df2append['sequence'] = str(each_fasta.seq).split()
                    output_df = output_df._append(df2append)
    if key != 'default':
        output_df_name = output_df_name.rsplit('.',maxsplit=1)[0]+'_'+key+'.csv'



    output_df.index = np.arange(len(output_df))
    output_file.close()
    output_df.to_csv(output_df_name)
    return output_df


def extra_alignment(input_folder:str = r'D:\my_experi\structure_prediction\structure_predicted_tylcv\TYLCV_CP_NbTo5_standard',
                    align_reference_pdb = 'default',
                    align_reference_chain = 'default2',
                    remove_reference_chain_after_alignment = 'True',
                    only_align = False):
    '''
    :param input_folder: 要进行extra_alignment的文件夹。
    :param align_reference_pdb: 作为alignment reference的PDb文件。
    :param align_chain: #要进行alignment的参考chain
    :return:
    '''
    #首先将这个目录下的所有cif文件转化成为pdb格式
    os.chdir(input_folder)
    cif_to_pdb(input_folder)
    pdf_file_name = 'TYLCV_CP_NbTo5_standard_seed_51989_sample_0.pdb'
    #首先读取mode下所有的pdb文件
    pymol.cmd.reinitialize()
    for each_file in os.listdir():
        if each_file.endswith('.pdb') and not each_file.endswith('only.pdb'):
            pymol.cmd.load(each_file)
    #首先选择一个align_reference_pdb
    #如果指定了就用指定的。如果没有指定，则从文件夹中所有的PDF文件中随机挑选一个。
    if align_reference_pdb == 'default':
        #此时则从文件夹中所有的PDb文件中，随机挑选一个作为alignment reference。
        align_reference_pdb_name = ''
        #print('input_folder: ', input_folder)
        for each_file in os.listdir(input_folder):
            #print('each_file: ',each_file)

            if each_file.endswith('.pdb'):
                align_reference_pdb = each_file
            if each_file[0:-4].endswith('_0.pdb'):
                align_reference_pdb = each_file
                break
                #如果发现这个文件的结尾是_0，这就先选用这个文件作为alignment reference。
        print(f'未指定alignment reference PDB文件,自动选用{align_reference_pdb}作为alignment reference。')
    else:
        align_reference_pdb = align_reference_pdb

    #如果此时没有找到合适的pdb
    if align_reference_pdb == 'default':
        print('该文件夹下没有pdb文件')
        return '该文件夹下没有pdb文件'

    align_reference_pdb_name = align_reference_pdb.split('.pdb')[0]
    #创制一个effector object 用于align
    #需要文件中存在一个，末尾编号为0的文件

    #接下来需要确定用于alignment的chain,默认情况下是A chain
    if align_reference_chain == 'default2':
        chains = pymol.cmd.get_chains(align_reference_pdb_name)
        if 'T' in chains:
            align_reference_chain = 'T'
        elif 'A' in chains:
            align_reference_chain = 'A'
        else:
            align_reference_chain = 'X'




    alignment = pymol.cmd.extra_fit(selection='n. ca')
    print('alignment',alignment)
    pymol.cmd.save(f'{input_folder}/merge_bulk_alignment.pse')
    #pymol.cmd.create('effector','c. A and *_0')
    pymol.cmd.create('effector', rf'c. {align_reference_chain} and {align_reference_pdb_name}')
    #pymol.cmd.save('a.pse')

    #alignment
    #cmd.extra_fit(selection = 'n. ca', reference='effector', method='super')
    pymol.cmd.extra_fit(selection='n. ca', reference='effector')
    #居中 effector
    pymol.cmd.orient('effector')
    if remove_reference_chain_after_alignment:
        pymol.cmd.remove(f'c. {align_reference_chain} and not effector') #删除掉effector object以外的chain A
    pymol.cmd.util.color_objs('all') #按照object进行染色。
    pymol.cmd.bg_colour('white') #将背景调为白色。
    pymol.cmd.zoom()
    #pymol.cmd.show('surface')
    pymol.cmd.set('transparency',0.75)
    os.makedirs('./png',exist_ok=True)
    pymol.cmd.png('./png/merge_png_objects.png',dpi= 150, ray = 0)
    pymol.cmd.spectrum("b",selection=("all"),quiet=0) #根据b-factor进行染色
    pymol.cmd.png('./png/merge_png_b_factor.png', dpi= 150, ray = 0)

    pymol.cmd.save(f'{input_folder}/merge.pse')
    #确定输出的文件夹。


    #然后添加interface 以及输出图片
    for each_file in os.listdir(input_folder):
        pymol.cmd.reinitialize()
        if not each_file.endswith('.pdb') or each_file.endswith('only.pdb'):
            continue

        #if each_file.endswith('.pdb') and not each_file.endswith('only.pdb'):
        pymol.cmd.reinitialize()
        pymol.cmd.load(each_file)
        pdb_name = each_file.rsplit('.pdb',maxsplit=1)[0]
        out_put_name =  pdb_name + '.pse'
        out_put_name2 = pdb_name + '_nb_only.pdb'
        out_put_name3 = pdb_name + '_dock.pse'
        out_put_name4 = pdb_name + '_bfactor.pse'
        out_put_figure1 = pdb_name + '.png'
        out_put_figure2 = pdb_name + '_dock.png'
        pymol.cmd.set_name(pdb_name,'fold0') #更改对象的名称
        #加氢
        pymol.cmd.h_add("fold0")
        pymol.cmd.sort("fold0 extend 1")
        #pymol.cmd.run(interfaceResidues())
        #pymol.cmd.run("D:\my_script\Pymol\InterfaceResidues_my.py")




        if only_align == False:
            try:
                pymol.cmd.run(ppi(pdb_name='fold0',ChA="T",ChB="B"))
            except AttributeError:
                pass
            try:

                pymol.cmd.run(interfaceResidues(cmpx = '*', cA='c. T', cB='c. B'))
                #interfaceResidues(cmpx='*', cA='c. A', cB='c. B', cutoff=1.0, selName="interface")
            except AttributeError:
                pass

            os.makedirs('png',exist_ok=True)
            os.makedirs('b_factor', exist_ok=True)

            pymol.cmd.bg_colour('white')
            pymol.cmd.color('orange', 'inter_1')
            pymol.cmd.color('blue', 'inter_2')
            pymol.cmd.util.cnc("inter_1", _self=pymol.cmd)
            pymol.cmd.util.cnc("inter_2", _self=pymol.cmd)

            #调整颜色
            pymol.cmd.show('surface')
            pymol.cmd.set('transparency', 0.75)
            pymol.cmd.show('cartoon')
            pymol.cmd.set('cartoon_transparency', 0)
            pymol.cmd.zoom()
            pymol.cmd.save(out_put_name)
            pymol.cmd.png('./png/'+out_put_figure1,dpi= 160)


            #然后提取出effector 画互作图
            pymol.cmd.reinitialize()
            pymol.cmd.load(each_file)
            pymol.cmd.set_name(pdb_name, 'fold0')  # 更改对象的名称
            # 加氢
            pymol.cmd.h_add("fold0")
            pymol.cmd.sort("fold0 extend 1")
            #
            pymol.cmd.extract('Target','fold0 and c. T')
            pymol.cmd.show_as('cartoon', 'fold0')
            pymol.cmd.show_as('surface','Target')
            pymol.cmd.show('cartoon', 'Target')
            pymol.cmd.color('palegreen', 'fold0')
            pymol.cmd.color('lightpink','Target')
            pymol.cmd.set('cartoon_transparency', 0)
            pymol.cmd.set('transparency', 0.54)
            pymol.cmd.bg_colour('white')  # 将背景调为白色。
            pymol.cmd.save(out_put_name3)
            pymol.cmd.png('./png/' + out_put_figure2, dpi=160)

            #然后在保存一份 b factor
            pymol.cmd.show_as('cartoon', 'Target')
            pymol.cmd.set('cartoon_transparency', 0)
            pymol.cmd.bg_colour('white')
            pymol.cmd.spectrum('b')
            pymol.cmd.save('./b_factor/'+out_put_name4)




            #删掉effector 只保留nanobody
            pymol.cmd.reinitialize()
            pymol.cmd.load(each_file)
            pymol.cmd.remove('c. T')
            os.makedirs(f'{input_folder}/supp_pdb',exist_ok=True)
            os.chdir(f'{input_folder}/supp_pdb')
            pymol.cmd.save(out_put_name2, format='pdb')
            os.chdir(input_folder)


def simple_alignment(align_pdb_path1 = '',
                     align_pdb_path2 = '',
                     output_pse_name = ''):
    #这个程序就是简单地将两个pdb文件进行alignment。
    pdb_name1 = align_pdb_path1.split('/')[-1].split('.pdb')[0]
    pdb_name2 = align_pdb_path2.split('/')[-1].split('.pdb')[0]
    pymol.cmd.reinitialize()

    pymol.cmd.load(align_pdb_path1)
    pymol.cmd.set_name(pdb_name1,'Binder')
    pymol.cmd.load(align_pdb_path2)

    pymol.cmd.extra_fit(selection='n. ca')
    pymol.cmd.show_as('cartoon', 'all')
    pymol.cmd.set('cartoon_transparency', 0)
    pymol.cmd.bg_colour('black')
    pymol.cmd.spectrum('b')
    pymol.cmd.save(output_pse_name)



def data_prepare_and_analyze_4_IgGM(input_pdb_or_cif = r"D:\my_script\Nanobody\IgGM\My_cases\tobrfv_cpv2_20250614\input\tobrfv_cpv2_nb35_model_0.pdb",
                                    mut_cdr = 'none',output_fasta_name = 'default',):
    '''
    这个程序的目的是为IgGM准备输入文件
    pdb文件只能包含两条链 chain,一条是nanobody，另一条是antigen(effector)
    输出一个fasta文件
    只需要将一个 antigen-nanobody的pdb文件放在input文件夹下面，就可以进行下面的内容
    :param input_pdb_or_cif:
    :param cif_or_pdb:  输入文件是cif格式还是pdb格式。如果是cif，是要首先将格式转化成为pdb。
    :param mut_cdr: 突变哪一个CDR?  none mut_cdr1 mut_cdr2 mut_cdr3 mut_cdr123
    把IgGM输出的文件，放在output文件夹下面
    :return:
    '''
    base_path = input_pdb_or_cif.rsplit('\\',maxsplit=1)[0]
    os.chdir(base_path)
    input_pdb_or_cif_name = input_pdb_or_cif.rsplit('\\',maxsplit=1)[1]
    input_pdb_name = input_pdb_or_cif_name.rsplit('.',maxsplit=1)[0] + '.pdb'
    if output_fasta_name == 'default':
        output_fasta_name = input_pdb_or_cif_name.rsplit('.',maxsplit=1)[0] + '.fasta'
        output_origin_fasta_name = input_pdb_or_cif_name.rsplit('.', maxsplit=1)[0] + '_original.fasta'
    #不管输入文件是什么，先将这个文件夹下所有可能的cif转化成为pdb
    cif_to_pdb(input_folder=base_path)
    # 然后从pdb文件中获取antigen_sequence nanobody_sequence
    #并对pdb文件进行重新命名chain

    antigen_sequence, nanobody_sequence = get_sequence_from_nanobody_antigen_pdb(pdb_path=input_pdb_name,
                                                                                 return_antigen_nanobody_sequence=True,
                                                                                 output_fasta=False,change_chain_name=True)
    my_nanobody = My_nanobody(sequence=nanobody_sequence,stringent=1,silent = 1)
    std_sequence = my_nanobody.standized_nb()
    if mut_cdr == 'none':
        mutant_nb_sequence = my_nanobody.sequence
    if mut_cdr == 'mut_cdr1':
        mutant_nb_sequence = my_nanobody.mut_cdr1()
    if mut_cdr == 'mut_cdr2':
        mutant_nb_sequence = my_nanobody.mut_cdr2()
    if mut_cdr == 'mut_cdr3':
        mutant_nb_sequence = my_nanobody.mut_cdr3()
    if mut_cdr == 'mut_cdr123':
        mutant_nb_sequence = my_nanobody.mut_cdr123()
    #然后根据突变序列，写入fasta文件
    os.chdir(base_path)
    output_fasta = open(output_fasta_name,'w')
    output_fasta.write('>H\n')
    output_fasta.write(mutant_nb_sequence+'\n')
    output_fasta.write('>A\n')
    output_fasta.write(antigen_sequence+'\n')
    output_fasta.close()

    output_origin_fasta = open(output_origin_fasta_name,'w')
    output_origin_fasta.write('>H\n')
    output_origin_fasta.write(nanobody_sequence+'\n')
    output_origin_fasta.write('>A\n')
    output_origin_fasta.write(antigen_sequence+'\n')
    output_origin_fasta.close()
    #创建output 文件夹
    os.chdir('..')
    if 'output' not in os.listdir():
        os.mkdir('output')
    #接下来开始考虑是否要对结果进行分析。
    #首先判断output文件夹当中是否含有内容
    os.chdir('output')

    if len(os.listdir()) > 1:
        #代表此时该文件夹当中存在信息,那么就可以开始进行信息处理。
        #首先将所有的fast文件进行汇总成一个文件。
        fasta_merge(input_folder=os.getcwd())
        output_df = fasta_merge(input_folder=os.getcwd(),key='H')
        #然后将所有的PDF文件，merge成为一个文件。
        #然后去into文件夹中获取原始序列
        os.chdir(base_path)
        for each_file in os.listdir():
            if each_file.endswith('original.fasta'):
                original_fasta_df = fasta_2_dataframe(each_file)
                output_df = output_df._append(original_fasta_df)

        os.chdir(base_path+'\\..\\output')
        output_df.to_csv('merged_fasta_H_all.csv')

        extra_alignment(input_folder=os.getcwd())
        pass

def pdb_align(input_folder = r'D:\my_script\Nanobody\RFantibody\rf2_output_example',
              nanobody_standardized = True,
              is_nanobody_antigen = True,
              is_output_from_rf2 = True):
    #这个脚本的目的是统计一个文件夹下所有的pdb文件，汇总序列
    os.chdir(input_folder)
    pdb_file_list = []
    output_df = pd.DataFrame()
    for each_file in os.listdir():
        if each_file.endswith('.pdb') and '..' not in each_file:
            pdb_file_list.append(each_file)
    print('pdb_file_list',pdb_file_list)
    for each_pdb_file in pdb_file_list:
        if is_nanobody_antigen:
            print('each_pdb_file', each_pdb_file)
            antigen_sequence, nanobody_sequence = get_sequence_from_nanobody_antigen_pdb(pdb_path = each_pdb_file,
                                                                                         change_chain_name = False, output_fasta = False, return_antigen_nanobody_sequence = True)
            df2append = pd.DataFrame()
            df2append['pdb_file'] = each_pdb_file.split()
            df2append['antigen_sequence'] = antigen_sequence.split()
            df2append['nanobody_sequence'] = nanobody_sequence.split()

            if nanobody_standardized:
                my_nanobody = My_nanobody(nanobody_sequence,stringent= 0)
                nanobody_std_sequence = my_nanobody.standized_nb()
                df2append['nanobody_sequence'] = nanobody_std_sequence.split()



            #如果这个pdb文件是rf2的输出结果。那么还可以从pdb文件中读取出模型的结果。
        if is_output_from_rf2:
            pdb_file = open(each_pdb_file,mode='r').readlines()
            for eachline in pdb_file:
                if eachline.startswith('SCORE'):
                    key = eachline.split(' ')[1].split(':')[0]
                    value = eachline.split(' ')[-1]
                    df2append[key] = value.split()


        output_df = output_df._append(df2append)
    output_df.index = np.arange(len(output_df))
    output_df.to_csv('pdb_summary.csv')
    return output_df
    pass
def rf2antibody_filtering(input_folder = r'D:\my_script\Nanobody\RFantibody\rf2_output_example',percentage_chose = 0.1,
                          output_folder = ''):
    '''
    :param input_folder:
    :param percentage_chose: 选择前百分比进行下一轮
    :return:
    '''

    #这个脚本是根据rf2antibody计算这个结果，进行初步的筛选，并将符合预期的PDF文件放到一个新的文件夹里面。
    os.chdir(input_folder)
    out_folder_name = output_folder.split('/')[-1]
    os.chdir('..')
    if out_folder_name not in os.listdir():
        os.mkdir(out_folder_name)
    os.chdir(input_folder)
    #首先读取pdb_summary.csv文件
    pdb_summary_df = pd.read_csv('pdb_summary.csv',index_col=0)
    #根据target_aligned_antibody_rmsd 进行排序
    pdb_summary_df['interaction_pae_plus_antibody_rmsd'] = pdb_summary_df['interaction_pae'] + pdb_summary_df['target_aligned_antibody_rmsd'] + pdb_summary_df['target_aligned_cdr_rmsd']
    pdb_summary_df = pdb_summary_df.sort_values(by = 'interaction_pae_plus_antibody_rmsd',ascending=True)
    pdb_summary_df['rank'] = np.arange(len(pdb_summary_df))
    pdb_summary_df['rel_rank'] = pdb_summary_df['rank']/len(pdb_summary_df)

    pdb_summary_filter_df = pdb_summary_df.loc[pdb_summary_df['rel_rank'] <= percentage_chose,:]
    pdb_summary_filter_df.to_csv('pdb_summary_filter.csv')


    for each_line in pdb_summary_filter_df.index:
        os.chdir(input_folder)
        file_name = pdb_summary_filter_df.loc[each_line,'pdb_file']
        if sys.platform.startswith('linux'):
            destination_file = f'../{out_folder_name}/'+file_name
        elif sys.platform.startswith('win'):
            destination_file = f'..\\{out_folder_name}\\'+file_name
        shutil.copy(file_name, destination_file)


    pass


def generating_json_as_af3_input(sequences_list = ['MSYTIATPSQFVFLSSAWADPIELINLCTNSLGNQFQTQQARTTVQRQFSEVWKPVPQVTVRFPDSGFKVYRYNAVLDPLVTALLGAFDTRNRIIEVENQANPTTAETLDATRRVDDATVAIRSAINNLVVELVKGTGLYNQSTFESASGLQWSSAPAS',
                                                   'QVQLVESGGGLVQAGGSLRLSCAASGFDFSKAWMGWFRQAPGKEREFVAAISPDGKESYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGFADGKGGGEDYWGQGTQVTVS'],
                                 input_name = 'BBTV',types = 'antigen_nanobody',
                                 modelSeeds = 1,
                                 output_folder = ''):
    """
    创建AlphaFold 3兼容的JSON输入文件
    参数:
    sequences -- 蛋白质序列列表 (单体) 或字典列表 (复合物)
    output_file -- 输出JSON文件路径
    type：['antigen_nanobody','monomer']
    sequences_list的格式：
    antigen_nanobody type:
    [antigen,nanobody]

    monomer：
    [antigen]
    """
    os.chdir(output_folder)
    # 处理单体输入（单条序列）
    output_file = input_name + '.json'
    if types == 'antigen_nanobody':
        data = {"name":input_name,
                'sequences':[{"protein":{"id":['T'],"sequence":sequences_list[0]}},
                             {"protein":{"id":['H'],"sequence":sequences_list[1]}}],
                "modelSeeds": [modelSeeds],
                "dialect": "alphafold3",
                "version": 1
                }


    if types == 'monomer':
        data = {"name":input_name,
                'sequences':[{"protein":{"id":['H'],"sequence":sequences_list[0]}}],
                "modelSeeds": [modelSeeds],
                "dialect": "alphafold3",
                "version": 1
                }


    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据
    #####################
    # 写入JSON文件
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"成功创建AlphaFold3输入文件: {output_file}")

def generating_json_for_af3_server(sequences_list = ['MSYTIATPSQFVFLSSAWADPIELINLCTNSLGNQFQTQQARTTVQRQFSEVWKPVPQVTVRFPDSGFKVYRYNAVLDPLVTALLGAFDTRNRIIEVENQANPTTAETLDATRRVDDATVAIRSAINNLVVELVKGTGLYNQSTFESASGLQWSSAPAS',
                                                   'QVQLVESGGGLVQAGGSLRLSCAASGFDFSKAWMGWFRQAPGKEREFVAAISPDGKESYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGFADGKGGGEDYWGQGTQVTVS'],
                                 input_name = 'BBTV',types = 'antigen_nanobody',
                                 modelSeeds = 1,
                                 output_folder = ''):
    """
    创建AlphaFold 3兼容的JSON输入文件  适配于af3 serveser
    sequences_list: antigen_nanobody
    """
    os.chdir(output_folder)
    # 处理单体输入（单条序列）
    output_file = input_name + '.json'
    if types == 'antigen_nanobody':
        data = [{
                 "name": input_name,
                 "modelSeeds": [str(modelSeeds)],
                 "sequences": [
                  {
                   "proteinChain": {
                    "sequence": sequences_list[0],
                    "count": 1,
                    "useStructureTemplate": True
                   }
                  },
                 {
                     "proteinChain": {
                         "sequence": sequences_list[1],
                         "count": 1,
                         "useStructureTemplate": True
                     }
                 }
                 ],
                 "dialect": "alphafoldserver",
                 "version": 1
                }
                ]


    if types == 'monomer':
        data = [{
                 "name": input_name,
                 "modelSeeds": [str(modelSeeds)],
                 "sequences": [
                  {
                   "proteinChain": {
                    "sequence": sequences_list[0],
                    "count": 1,
                    "useStructureTemplate": True
                   }
                  }
                 ],
                 "dialect": "alphafoldserver",
                 "version": 1
                }
                ]

    if types == 'homodimer':
        data = [{
                 "name": input_name,
                 "modelSeeds": [str(modelSeeds)],
                 "sequences": [
                  {
                   "proteinChain": {
                    "sequence": sequences_list[1],
                    "count": 2,
                    "useStructureTemplate": True
                   }
                  }
                 ],
                 "dialect": "alphafoldserver",
                 "version": 1
                }
                ]



    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据
    #####################
    # 写入JSON文件
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"成功创建AlphaFold3_server的输入文件: {output_file}")

def generating_json_for_helixfold(sequences_list = ['MSYTIATPSQFVFLSSAWADPIELINLCTNSLGNQFQTQQARTTVQRQFSEVWKPVPQVTVRFPDSGFKVYRYNAVLDPLVTALLGAFDTRNRIIEVENQANPTTAETLDATRRVDDATVAIRSAINNLVVELVKGTGLYNQSTFESASGLQWSSAPAS',
                                                   'QVQLVESGGGLVQAGGSLRLSCAASGFDFSKAWMGWFRQAPGKEREFVAAISPDGKESYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGFADGKGGGEDYWGQGTQVTVS'],
                                 input_name = 'BBTV',types = 'antigen_nanobody',
                                 modelSeeds = 1,
                                 output_folder = '',
                                  recycle = 10,ensemble = 1):
    """
    创建AlphaFold 3兼容的JSON输入文件  适配于af3 serveser
    sequences_list: antigen_nanobody
    """
    os.chdir(output_folder)
    # 处理单体输入（单条序列）
    output_file = input_name + '.json'
    if types == 'antigen_nanobody':
        data = {
                  "job_name": input_name,
                'recycle':recycle,
                'ensemble':ensemble,
                  "entities": [
                    {
                      "type": "protein",
                      "sequence": sequences_list[0],
                      "count": 1
                    },
                    {
                      "type": "protein",
                      "sequence": sequences_list[1],
                      "count": 1
                    }
                  ]
                }


    if types == 'monomer':
        data = {
                  "job_name": input_name,
            'recycle': recycle,
            'ensemble': ensemble,
                  "entities": [
                    {
                      "type": "protein",
                      "sequence": sequences_list[0],
                      "count": 1
                    },
                  ]
                }

    if types == 'homodimer':
        data = {
                  "job_name": input_name,
            'recycle': recycle,
            'ensemble': ensemble,
                  "entities": [
                    {
                      "type": "protein",
                      "sequence": sequences_list[1],
                      "count": 2
                    },
                  ]
                }

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据
    #####################
    # 写入JSON文件
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"成功创建helixfold的输入文件: {output_file}")
    return data

def generating_json_for_proteinx(sequences_list = ['MSYTIATPSQFVFLSSAWADPIELINLCTNSLGNQFQTQQARTTVQRQFSEVWKPVPQVTVRFPDSGFKVYRYNAVLDPLVTALLGAFDTRNRIIEVENQANPTTAETLDATRRVDDATVAIRSAINNLVVELVKGTGLYNQSTFESASGLQWSSAPAS',
                                                   'QVQLVESGGGLVQAGGSLRLSCAASGFDFSKAWMGWFRQAPGKEREFVAAISPDGKESYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGFADGKGGGEDYWGQGTQVTVS'],
                                 input_name = 'BBTV',types = 'antigen_nanobody',
                                 modelSeeds = 1,
                                 output_folder = ''):
    """
    创建AlphaFold 3兼容的JSON输入文件  适配于af3 serveser
    sequences_list: antigen_nanobody
    """
    os.chdir(output_folder)
    # 处理单体输入（单条序列）
    output_file = input_name + '.json'
    if types == 'antigen_nanobody':
        data = [
                {
                    "name": input_name,
                    "covalent_bonds": [],
                    "sequences": [
                        {
                            "proteinChain": {
                                "count": 1,
                                "sequence": sequences_list[0],
                                "modifications": [],
                                "msa": {

                                }
                            }
                        },
                        {
                            "proteinChain": {
                                "count": 1,
                                "sequence": sequences_list[1],
                                "modifications": [],
                                "msa": {

                                }
                            }
                        }
                    ]
                }
            ]


    if types == 'monomer':
        data = [
            {
                "name": input_name,
                "covalent_bonds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "count": 1,
                            "sequence": sequences_list[0],
                            "modifications": [],
                            "msa": {

                            }
                        }
                    }

                ]
            }
        ]

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据
    #####################
    # 写入JSON文件
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"成功创建proteinx的输入文件: {output_file}")
    return data

#利用prodigy评估两个蛋白之间的互作强度
def run_prodigy(pdb_file, chain1, chain2):
    """
    使用Prodigy预测两个链之间的相互作用强度
    返回包含预测结果的字典
    """
    # 构建命令
    cmd = [
        "prodigy",
        pdb_file,
        "--selection", chain1, chain2
        #"--json"  # 获取JSON格式输出
    ]

    try:
        # 运行Prodigy命令
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # 解析JSON输出

        output = str(result.stdout)
        predicted_binding_affinity= output.split(sep = 'binding affinity (kcal.mol-1):')[-1].split('\n')[0].strip()
        predicted_dissociation_constant = output.split(sep='constant (M) at 25.0˚C:')[-1].split('\n')[0].strip()



        print('prodigy',str(output))
        return predicted_binding_affinity,predicted_dissociation_constant


    except subprocess.CalledProcessError as e:
        print(f"Prodigy运行错误: {e.stderr}")
        print('no result out')
        return 0,0
        #sys.exit(1)
    except FileNotFoundError:
        print("错误: 未找到Prodigy。请确保已安装Prodigy并添加到PATH")
        print("安装方法: pip install prodigy")
        #sys.exit(1)


def run_foldx_AnalyseComplex(pdb_file, chain1, chain2):
    """
    使用foldx预测两个链之间的相互作用强度
    返回包含预测结果的字典
    """
    foldx_exe = r'/public-supool/home/gaolab/haochengzhu/Foldx4/foldx_20261231'
    rotabase = r'/public-supool/home/gaolab/haochengzhu/Foldx4/rotabase.txt'
    # 构建命令
    foldx_repairpdb_cmd = [
        foldx_exe,
        "--command=RepairPDB",
        rf"--rotabaseLocation={rotabase}",
        f'--pdb={pdb_file}'
    ]

    foldx_Optimize_cmd = [
        foldx_exe,
        "--command=Optimize",
        rf"--rotabaseLocation={rotabase}",
        f'--pdb={pdb_file}'
    ]

    cmd1 = [
        foldx_exe,
        "--command=AnalyseComplex",
        f'--pdb={pdb_file}',
        f"--analyseComplexChains={chain2},{chain1}",
        "--complexWithDNA=false"
    ]
    print('cmd',cmd1)

    cmd2 = [
        foldx_exe,
        "--command=AnalyseComplex",
        f'--pdb=Optimized_{pdb_file}',
        f"--analyseComplexChains={chain2},{chain1}",
        "--complexWithDNA=false"
    ]
    print('cmd',cmd2)



    # 运行foldx命令
    #result0 = subprocess.run(foldx_repairpdb_cmd, capture_output=True, text=True, check=True)
    #output1 = str(result0.stdout)
    #print('foldx_output', output1)

    result0 = subprocess.run(foldx_Optimize_cmd, capture_output=True, text=True, check=True)
    output0 = str(result0.stdout)
    print('foldx_Optimize_output\n', output0)

    result1 = subprocess.run(cmd1, capture_output=True, text=True, check=True)
    output1 = str(result1.stdout)
    print('foldx_output\n', output1)

    result2 = subprocess.run(cmd2, capture_output=True, text=True, check=True)
    output2 = str(result2.stdout)
    print('foldx_output_opt\n', output2)

    foldx_Interaction_Energy=output1.split('interaction between')[-1].split('Total          =')[-1].split('\n')[0].strip()
    foldx_Interaction_Energy_opt = output2.split('interaction between')[-1].split('Total          =')[-1].split('\n')[
        0].strip()


    print('foldx_Interaction_Energy: ',foldx_Interaction_Energy)
    print('foldx_Interaction_Energy_opt: ',foldx_Interaction_Energy_opt)

    #然后去读文件夹下面的Interface_Residues_Optimized文件
    for each_file in os.listdir('.'):
        if not each_file.startswith('Interface_Residues_Optimized'):
            continue
        inte_resi = open(each_file,'r').readlines()
        for i in range(len(inte_resi)):
            if inte_resi[i].startswith('interface residues between'):
                inter_info = inte_resi[i+1][:-1]
    binder_interaction_residue = ''
    target_interaction_residue = ''
    #print('inter_info',inter_info)
    print(inter_info.split(sep='\t'))
    for each_residue in inter_info.split(sep='\t'):
        if each_residue == '':
            continue
        if each_residue[1] == chain2:
            binder_interaction_residue += (each_residue[1:] + ' ')
        if each_residue[1] == chain1:
            target_interaction_residue += (each_residue[1:] + ' ')

    binder_interaction_residue = binder_interaction_residue[:-1]
    target_interaction_residue = target_interaction_residue[:-1]

    return foldx_Interaction_Energy,foldx_Interaction_Energy_opt,binder_interaction_residue,target_interaction_residue
        #predicted_binding_affinity= output.split(sep = 'binding affinity (kcal.mol-1):')[-1].split('\n')[0].strip()
        #predicted_dissociation_constant = output.split(sep='constant (M) at 25.0˚C:')[-1].split('\n')[0].strip()


####################################
def barcode_distance_calculator(barcode1='ATATAT', barcode2='TATATA'):
    distance = {'AA': 0, 'TT': 0, 'CC': 0, 'GG': 0,
                'AG': 1, 'GA': 1, 'TC': 1, 'CT': 1,
                'AT': 3, 'TA': 3, 'GC': 3, 'CG': 3,
                'AC': 3, 'CA': 3, 'GT': 3, 'TG': 3,
                'RY': 3, 'YR': 3, 'RR': 0, 'YY': 0,
                'RA': 0, 'RT': 3, 'RC': 3, 'RG': 0,
                'YA': 3, 'YT': 0, 'YC': 0, 'YG': 3,
                'AR': 0, 'TR': 3, 'CR': 3, 'GR': 0,
                'AY': 3, 'TY': 0, 'CY': 0, 'GY': 3, }

    base = ['A', 'T', 'C', 'G', 'R', 'Y']

    dis = 0
    if len(barcode1) == len(barcode2):
        for i in range(len(barcode1)):
            dis += distance[barcode1[i] + barcode2[i]]

    return dis


###########################################
def CAI_calculator(input_seq='ATCGCACTGCTGATGATGATGATGATGATGATG',
                   organism='Nicotiana benthamiana',
                   Codon_Usage_Tables_folder=r'/public-supool/home/gaolab/haochengzhu/Codon_Usage_Tables'):
    # 实现密码子适应指数CAI (codon_adaptation_index)的计算
    # 首先要调取并读取该物种的密码子利用表
    # os.chdir(Codon_Usage_Tables_folder)
    for each_file in os.listdir(Codon_Usage_Tables_folder):
        # print(each_file)
        if organism in each_file:
            lxml_file = each_file
            tree = etree.parse(Codon_Usage_Tables_folder + '/' + lxml_file)
            root = tree.getroot()
            for element in root.iter():
                if 'frequencies' in element.attrib.keys():
                    # print(element.tag, element.attrib, element.text)
                    codon_frequencies = element.attrib['frequencies']
                    pass
            break
    # 然后需要根据codon_frequencies，获得一个dict
    codon_frequency_rscu_dict = {}
    for each_codon in codon_frequencies.split(sep=';'):
        codon_frequency_rscu_dict[each_codon.split(':')[0]] = float(each_codon.split(':')[1])
    # print('organism =',organism)
    # print(codon_frequency_rscu_dict)

    # 接下来，可以计算CAI了
    raw_cai = 1
    if len(input_seq) % 3 == 0:
        # 说明输入序列的长度是三的倍数，可以正常进行分析。
        pass
    else:
        print('输入序列的长度并非是三的倍数对序列末尾进行了修正')
        input_seq = input_seq[0:-len(input_seq) % 3]
    for i in range(0, len(input_seq), 3):
        codon = input_seq[i:i + 3]
        if codon in codon_frequency_rscu_dict.keys():
            RSCU = codon_frequency_rscu_dict[codon]
            # 计算Relative Synonymous Codon Usage (RSCU)
            raw_cai = raw_cai * RSCU
        else:
            print('出现非正常的密码子:', codon)
    cai = raw_cai ** (3 / len(input_seq))
    return cai
    # print('CAI',cai,raw_cai)


############################
def anchor_bar_generator0(original_seq='ATGCTGCAAAGG',
                         iteration=240,
                         GC_min=0.38, GC_max=0.59,
                         bsa1_discard=True,
                         organism='Nicotiana benthamiana',
                         sort_by='cai',
                         anchor_bar_distance=4,
                         original_seq_type='nucl'):
    # 该还函数的目的是能够产生出anchor barcode
    # 首先获取密码子表
    # original_seq_type, prot or nucl
    standard_table = CodonTable.unambiguous_dna_by_id[1].forward_table
    aa2seq_dict = {}
    for each_codon in standard_table.keys():
        aa = standard_table[each_codon]
        if aa in aa2seq_dict.keys():
            aa2seq_dict[aa].append(each_codon)

        else:
            aa2seq_dict[aa] = [each_codon]
    # 添加上终止密码子
    aa2seq_dict['*'] = ['TAG', 'TGA', 'TAA']
    if original_seq_type == 'nucl':
        aa_seq = str(Seq(original_seq).translate())
    else:
        aa_seq = original_seq
    dna_seq_df = pd.DataFrame()
    # 创建一个用于承接随机产生的dna序列的df
    for i in range(iteration):
        df2append = pd.DataFrame()
        seq0 = ''
        for each_aa in aa_seq:
            # print(each_aa)
            seq0 += random.choice(aa2seq_dict[each_aa])
        # print('seq0',seq0)
        gc_content = (seq0.count('G') + seq0.count('C')) / len(seq0)
        if 'GGTCTC' in seq0 or 'GAGACC' in seq0 or 'CGTCTC' in seq0 or 'GAGACG' in seq0:
            # BsaI and BsmBI
            bsai = True
        elif 'GTCTC' in seq0 or 'GGTCT' in seq0 or 'AGACC' in seq0 or 'GAGAC' in seq0:
            # BsaI (5bp)
            bsai = True
        elif seq0[0:4] in ['TCTC', 'GACC', 'GGTC', 'GAGA'] or seq0[-4:] in ['TCTC', 'GACC', 'GGTC', 'GAGA']:
            # BsaI (4bp)
            bsai = True
        else:
            bsai = False
        df2append['sequece'] = seq0.split()
        df2append['bsai'] = bsai
        df2append['gc_content'] = gc_content
        df2append['cai'] = CAI_calculator(input_seq=seq0, organism=organism)
        dna_seq_df = dna_seq_df._append(df2append)
    dna_seq_df.index = np.arange(len(dna_seq_df))
    #如果最高的GC含量还是不能达到下线的话
    max_gc = np.max(dna_seq_df['gc_content'])
    if max_gc <= GC_min:
        dna_seq_df = dna_seq_df.loc[dna_seq_df['gc_content'] >= max_gc-0.015, :]
    else:
        dna_seq_df = dna_seq_df.loc[dna_seq_df['gc_content'] >= GC_min, :]

    dna_seq_df = dna_seq_df.loc[dna_seq_df['gc_content'] <= GC_max, :]
    if bsa1_discard:
        dna_seq_df = dna_seq_df.loc[dna_seq_df['bsai'] == False, :]

    # 在进行一次编号,根据参数(cai or gc_content)进行排序
    dna_seq_df = dna_seq_df.sort_values(by=sort_by, ascending=False)
    # 隐藏参数，当gc最大值低于0.5时，以gc值来排序

    if np.max(dna_seq_df['gc_content']) < 0.45:
        dna_seq_df = dna_seq_df.sort_values(by='gc_content', ascending=False)

    dna_seq_df.index = np.arange(len(dna_seq_df))
    # 接下来从这些序列中 对比anchor_bar_distance，选出最合适的anchor_bar
    for each_line in dna_seq_df.index:
        seq = dna_seq_df.loc[each_line, 'sequece']
        if barcode_distance_calculator(seq, original_seq) >= anchor_bar_distance:
            print('anchor_barcode: ', seq)
            return seq
            break


def anchor_bar_generator(original_seq='ATGCTGCAAAGG',
                         iteration=350,
                         GC_min=0.40, GC_max=0.59,
                         bsa1_discard=True,
                         organism='Nicotiana benthamiana',
                         sort_by='cai',
                         anchor_bar_distance=4,
                         original_seq_type='nucl'):
    # 该还函数的目的是能够产生出anchor barcode
    # 首先获取密码子表
    # original_seq_type, prot or nucl
    standard_table = CodonTable.unambiguous_dna_by_id[1].forward_table
    aa2seq_dict = {}
    for each_codon in standard_table.keys():
        aa = standard_table[each_codon]
        if aa in aa2seq_dict.keys():
            aa2seq_dict[aa].append(each_codon)

        else:
            aa2seq_dict[aa] = [each_codon]
    # 添加上终止密码子
    aa2seq_dict['*'] = ['TAG', 'TGA', 'TAA']
    if original_seq_type == 'nucl':
        aa_seq = str(Seq(original_seq).translate())
    else:
        aa_seq = original_seq
    dna_seq_df = pd.DataFrame()
    # 创建一个用于承接随机产生的dna序列的df

    # 只要没算出来，就一直算
    get_the_sequence = 0
    while get_the_sequence == 0:
        for i in range(iteration):
            df2append = pd.DataFrame()
            seq0 = ''
            for each_aa in aa_seq:
                # print(each_aa)
                seq0 += random.choice(aa2seq_dict[each_aa])
            # print('seq0',seq0)
            gc_content = (seq0.count('G') + seq0.count('C')) / len(seq0)
            if 'GGTCTC' in seq0 or 'GAGACC' in seq0 or 'CGTCTC' in seq0 or 'GAGACG' in seq0:
                # BsaI and BsmBI
                bsai = True
            elif 'GTCTC' in seq0 or 'GGTCT' in seq0 or 'AGACC' in seq0 or 'GAGAC' in seq0:
                # BsaI (5bp)
                bsai = True
            elif seq0[0:4] in ['TCTC', 'GACC', 'GGTC', 'GAGA'] or seq0[-4:] in ['TCTC', 'GACC', 'GGTC', 'GAGA']:
                # BsaI (4bp)
                bsai = True
            else:
                bsai = False
            df2append['sequece'] = seq0.split()
            df2append['bsai'] = bsai
            df2append['gc_content'] = gc_content
            df2append['cai'] = CAI_calculator(input_seq=seq0, organism=organism)
            dna_seq_df = dna_seq_df._append(df2append)
        dna_seq_df.index = np.arange(len(dna_seq_df))

        if bsa1_discard:
            dna_seq_df = dna_seq_df.loc[dna_seq_df['bsai'] == False, :]

        gc_max = np.max(dna_seq_df['gc_content'])
        while gc_max <= GC_min:
            GC_min = gc_max - 0.012

        dna_seq_df = dna_seq_df.loc[dna_seq_df['gc_content'] >= GC_min, :]
        dna_seq_df = dna_seq_df.loc[dna_seq_df['gc_content'] <= GC_max, :]

        '''
        # 在进行一次编号,根据参数(cai or gc_content)进行排序
        dna_seq_df = dna_seq_df.sort_values(by=sort_by, ascending=False)
        # 隐藏参数，当gc最大值低于0.5时，以gc值来排序

        if np.max(dna_seq_df['gc_content']) < 0.5:
            dna_seq_df = dna_seq_df.sort_values(by='gc_content', ascending=False)
        '''
        dna_seq_df.index = np.arange(len(dna_seq_df))
        # 接下来从这些序列中 对比anchor_bar_distance，选出最合适的anchor_bar
        for each_line in dna_seq_df.index:
            seq = dna_seq_df.loc[each_line, 'sequece']
            if barcode_distance_calculator(seq, original_seq) >= anchor_bar_distance:
                print('anchor_barcode: ', seq)

                get_the_sequence = 1  # while 循环停止
                return seq

        print('anchor_bar_generator: no sequence generated, repeat again')


def finding_dna_repeats(dna_sequence, min_length=9):
    """
    在DNA序列中查找长度大于min_length的重复序列

    参数:
    dna_sequence (str): 输入的DNA序列
    min_length (int): 最小重复长度阈值，默认为10bp

    返回:
    list: 包含重复序列信息的列表，每个元素为(重复序列, 起始位置1, 起始位置2)
    """
    dna_sequence = dna_sequence.upper()
    if not all(base in 'ATCG' for base in dna_sequence):
        print("错误: 序列包含非DNA碱基(A, T, C, G)")
        return
    if len(dna_sequence) < 22:  # 至少需要2*10+2个碱基才能有10bp的重复
        print("序列太短，无法包含长于10bp的重复序列")
        return

    repeats = []
    n = len(dna_sequence)

    # 检查所有可能的子串长度，从min_length到序列长度的一半
    for length in range(min_length , n // 2 + 1):
        # 检查所有可能的起始位置
        for i in range(0, n - length * 2 + 1):
            substring = dna_sequence[i:i + length]
            substring_rc = str(Seq(substring).reverse_complement())

            # 检查这个子串是否在序列的其他位置出现
            for j in range(i + length, n - length + 1):
                if dna_sequence[j:j + length] == substring or dna_sequence[j:j + length] == substring_rc:
                    repeats.append((substring, i, j))

    if len(repeats) > 0:
        print(f"发现 {len(repeats)} 个长于{min_length}bp的重复序列:")
        for repeat, pos1, pos2 in repeats:
            print(f"重复序列: {repeat}, 位置: {pos1}-{pos1 + len(repeat) - 1} 和 {pos2}-{pos2 + len(repeat) - 1}")
        return 1
    else:
        print(f"未发现长于{min_length}bp的重复序列")
        return 0


def gc_content_upper(seq0 = 'TCTTTCTGGGTCGACCTATTcAAACAGCCGGGTTTCGActtcTCTGTGTACTGGAACTTCTGGGTTGAAACTATG',
                     random_threshold = 0.44,
                     gc_content_upper = 0.60,
                     gc_content_lower = 0.535,):
    #这个脚本的目的，是提高氨基酸编码序列的GC含量
    #random 多少的概率 进行TC或AG变换

    seq0 = seq0.upper()
    origin_aa = str(Seq(seq0).translate())
    now_is_good = 0
    #random.seed = 1
    j = 0
    k = 0
    while now_is_good == 0:
        j += 1

        for i in range(len(seq0)):

            if seq0[i] == 'A' and random.uniform(0,1)> random_threshold:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'G'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1

            if seq0[i] == 'A' and random.uniform(0,1)> random_threshold:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'C'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1


            elif seq0[i] == 'T' and random.uniform(0,1)> random_threshold:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'C'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1

            elif seq0[i] == 'T' and random.uniform(0,1)> random_threshold:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'G'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1

            elif seq0[i] == 'G' and random.uniform(0,1) < random_threshold/1.5:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'A'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1
            elif seq0[i] == 'G' and random.uniform(0,1) < random_threshold/1.5:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'T'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1


            elif seq0[i] == 'C' and random.uniform(0,1)< random_threshold/1.5:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'T'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1

            elif seq0[i] == 'C' and random.uniform(0,1)< random_threshold/1.5:
                #print(random.uniform(0,1))
                seq1 = seq0[:i] + 'A'+ seq0[i+1:]
                origin_aa1 = str(Seq(seq1).translate())
                if origin_aa1 == origin_aa:
                    seq0 = seq1


        gc_content = (seq0.count('G') + seq0.count('C')) / len(seq0)
        print('seq0: ', seq0,'gc_content',gc_content,'random_threshold',random_threshold,'j:',j,'k:',k)
        if random_threshold < 0.002 :
            random_threshold = 0.46
            gc_content_lower -= 0.008

        if random_threshold > 0.98:
            random_threshold = 0.46
            gc_content_upper += 0.008
            #重置

        if gc_content > gc_content_upper:
            random_threshold += 0.0001
            now_is_good = 0
        elif gc_content < gc_content_lower:
            random_threshold -= 0.0001
            now_is_good = 0
        elif gc_content <= gc_content_upper and gc_content >= gc_content_lower:
            if 'GGTCTC' in seq0 or 'GAGACC' in seq0:
                now_is_good = 0
            elif seq0.startswith('TCTC') or seq0.startswith('GTCTC') or seq0.startswith('AGACC') or seq0.startswith('GACC'):
                now_is_good = 0
            elif seq0.endswith('GAGA') or seq0.endswith('GAGAC') or seq0.endswith('GGTC') or seq0.endswith('GGTCT'):
                now_is_good = 0
            else:
                finding_repeat = 1
                k += 1
                point1, point2 ,point3,point4,point5 = 65,200,400,1200,2000
                if k < point1:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 8)
                if k >= point1 and k < point2:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 9)
                if k >= point2 and k < point3:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 10)
                if k >= point3 and k < point4:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 11)
                if k >= point4 and k < point5:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 12)
                if k >= point5:
                    finding_repeat = finding_dna_repeats(seq0,min_length = 13)

                if finding_repeat == 0:
                    now_is_good = 1

    print('gc_content_upper:',seq0)
    return seq0


###################################################################################################
class CA_Select(Select):
    """选择器用于提取Cα原子"""

    def accept_atom(self, atom):
        return atom.get_name() == "CA"

'''
def load_structure(pdb_file):
    """加载PDB结构并提取Cα原子"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_atoms.append(residue["CA"])
    return structure, ca_atoms
'''

def extract_ca_atoms(pdb_file = None, structure = None, chain_id='H', start_res=0, end_res=None):
    """提取PDB文件中指定区域的Cα原子"""
    parser = PDBParser(QUIET=True)
    if structure == None:
        structure = parser.get_structure("structure", pdb_file)


    ca_atoms = []

    for model in structure:
        for chain in model:
            # 检查链ID是否匹配
            if chain_id and chain.id != chain_id:
                continue

            for residue in chain:
                # 获取残基编号
                res_id = residue.id[1]
                # 检查残基范围
                if start_res and res_id < start_res:
                    continue
                if (end_res and res_id > end_res):
                    continue

                # 添加Cα原子
                if residue.has_id("CA"):
                    ca_atoms.append(residue["CA"])

    return structure, ca_atoms

def calculation_nanobody_rsmd(pdb_file1, pdb_file2):
    """计算两个PDB文件之间的RMSD"""
    # 加载结构并提取Cα原子
    struture1, atoms1 = extract_ca_atoms(pdb_file = pdb_file1)
    struture2, atoms2 = extract_ca_atoms(pdb_file = pdb_file2)
    #print('struture2',struture2)


    # 检查原子数量是否匹配
    if len(atoms1) != len(atoms2):
        print(f"警告: 结构原子数不匹配 - {pdb_file1}: {len(atoms1)}, {pdb_file2}: {len(atoms2)}")
        print("将使用共同长度的子集进行比较")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]
    # 设置叠加器
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)  # 设置参考原子和移动原子

    # 计算并返回RMSD
    rmsd_all = sup.rms
    print(f"计算完成: {pdb_file1} 和 {pdb_file2} 之间的RMSD")
    print(f"主链Cα原子RMSD: {rmsd_all:.3f} Å")
    return rmsd_all


def calculate_partial_rmsd(partial_ca1, partial_ca2):
    """计算两个部分Cα原子列表之间的RMSD"""
    # 检查原子数量
    if len(partial_ca1) == 0 or len(partial_ca2) == 0:
        raise ValueError("未找到指定区域的Cα原子！请检查参数")

    # 确保原子数量匹配
    if len(partial_ca1) != len(partial_ca2):
        print(f"注意: 部分区域的Cα原子数量不同 - {len(partial_ca1)} vs {len(partial_ca2)}")
        print("将使用共同长度的子集进行比较")
        min_len = min(len(partial_ca1), len(partial_ca2))
        partial_ca1 = partial_ca1[:min_len]
        partial_ca2 = partial_ca2[:min_len]

    # 计算RMSD
    sum_sq = 0.0
    for atom1, atom2 in zip(partial_ca1, partial_ca2):
        diff = atom1.coord - atom2.coord
        sum_sq += sum(diff * diff)

    rsmd = (sum_sq / len(partial_ca1)) ** 0.5
    return rsmd

def calculation_nanobody_cdrs_rsmd(pdb_file1, pdb_file2):
    """计算两个nanobody PDB文件之间的nanobody"""
    #首先获取nanobody的序列
    nanobody1_seq,nanobody1_residue_numbers = get_the_sequence_from_single_pdb(input_pdb=pdb_file1,
                                     chain_id='H', silent=1)

    nanobody2_seq,nanobody2_residue_numbers = get_the_sequence_from_single_pdb(input_pdb=pdb_file2,
                                     chain_id='H', silent=1)

    # 加载结构并提取Cα原子
    structure1, atoms1 = extract_ca_atoms(pdb_file1)
    structure2, atoms2 = extract_ca_atoms(pdb_file2)

    # 检查原子数量是否匹配
    if len(atoms1) != len(atoms2):
        print(f"警告: 结构原子数不匹配 - {pdb_file1}: {len(atoms1)}, {pdb_file2}: {len(atoms2)}")
        print("将使用共同长度的子集进行比较")
        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]
    # 设置叠加器
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)  # 设置参考原子和移动原子
    # 将变换应用到整个结构2
    sup.apply(structure2.get_atoms())
    structure2.get_atoms

    #然后开始计算CDR区域的RSMD
    my_nanobody = My_nanobody(sequence=nanobody1_seq)
    cdr1_begin, cdr1_end = my_nanobody.obtain_cdr_border(cdr = 'cdr1')
    cdr2_begin, cdr2_end = my_nanobody.obtain_cdr_border(cdr = 'cdr2')
    cdr3_begin, cdr3_end = my_nanobody.obtain_cdr_border(cdr='cdr3')

    print('regions',cdr1_begin,cdr1_end,cdr2_begin,cdr2_end,cdr3_begin,cdr3_end)
    _, atoms1_cdr1 = extract_ca_atoms(structure = structure1,start_res=cdr1_begin, end_res=cdr1_end)
    _, atoms2_cdr1 = extract_ca_atoms(structure = structure2, start_res=cdr1_begin, end_res=cdr1_end)
    _, atoms1_cdr2 = extract_ca_atoms(structure = structure1, start_res=cdr2_begin, end_res=cdr2_end)
    _, atoms2_cdr2 = extract_ca_atoms(structure = structure2, start_res=cdr2_begin, end_res=cdr2_end)
    _, atoms1_cdr3 = extract_ca_atoms(structure = structure1, start_res=cdr3_begin, end_res=cdr3_end)
    _, atoms2_cdr3 = extract_ca_atoms(structure = structure2, start_res=cdr3_begin, end_res=cdr3_end)



    rsmd_cdr1 = calculate_partial_rmsd(atoms1_cdr1,atoms2_cdr1)
    rsmd_cdr2 = calculate_partial_rmsd(atoms1_cdr2, atoms2_cdr2)
    rsmd_cdr3 = calculate_partial_rmsd(atoms1_cdr3, atoms2_cdr3)

    print('rsmd_cdr1', rsmd_cdr1)
    print('rsmd_cdr2', rsmd_cdr2)
    print('rsmd_cdr3', rsmd_cdr3)


    return rsmd_cdr1,rsmd_cdr2,rsmd_cdr3,my_nanobody.len_cdr1,my_nanobody.len_cdr2,my_nanobody.len_cdr3


def concatenate_csv(concatenate_from_csv = '',concatenate_to_csv = '',
                    concatenate_key_from = '',concatenate_key_to = '',column_to_concatenate = [''],
                    new_csv = 'default'):
    concatenate_from_df = pd.read_csv(concatenate_from_csv,index_col=0)
    concatenate_to_df = pd.read_csv(concatenate_to_csv, index_col=0)

    concatenate_from_df.index = concatenate_from_df[concatenate_key_from]
    concatenate_to_df.index = concatenate_to_df[concatenate_key_to]

    for each_column in column_to_concatenate:
        concatenate_to_df[each_column] = concatenate_from_df[each_column]

    concatenate_to_df.index = np.arange((len(concatenate_to_df)))
    if new_csv == 'default':
        concatenate_to_df.to_csv(concatenate_to_csv)
    else:
        concatenate_to_df.to_csv(new_csv)


AA_MAP = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

'''
def changing_aa_in_pdb(input_pdb = '',
                       chain_to_change = 'H',
                       change_to_sequence = '',
                       output_pdb_path = ''):
'''



def analyze_secondary_structure(pdb_file, chain_id='B'):
    #返回结构中的α螺旋和βsheet的数量
    # 创建PDB解析器
    parser = PDBParser(QUIET=True)

    try:
        # 解析PDB文件
        structure = parser.get_structure("protein", pdb_file)
    except FileNotFoundError:
        print(f"错误: 文件 '{pdb_file}' 未找到")
        return
    except Exception as e:
        print(f"解析PDB文件时出错: {str(e)}")
        return

    # 获取目标链
    target_chain = None
    for model in structure:
        if chain_id in model:
            target_chain = model[chain_id]
            break

    if target_chain is None:
        print(f"警告: 链 '{chain_id}' 在PDB文件中不存在")
        return

    # 提取二级结构信息
    helices = []
    sheets = []
    current_helix = []
    current_sheet = []

    for residue in target_chain:
        # 跳过非氨基酸残基（如水分子、离子等）
        if not is_aa(residue):
            continue

        # 获取二级结构类型（来自PDB的HELIX/SHEET记录）
        ss_type = residue.xtra.get("SSE", " ")

        # 处理α螺旋（H、G、I）
        if ss_type in ['H', 'G', 'I']:
            if not current_helix and current_sheet:
                sheets.append(current_sheet)
                current_sheet = []
            current_helix.append(residue)
        # 处理β折叠（E、B）
        elif ss_type in ['E', 'B']:
            if not current_sheet and current_helix:
                helices.append(current_helix)
                current_helix = []
            current_sheet.append(residue)
        # 处理二级结构结束
        else:
            if current_helix:
                helices.append(current_helix)
                current_helix = []
            if current_sheet:
                sheets.append(current_sheet)
                current_sheet = []

    # 添加最后的片段
    if current_helix:
        helices.append(current_helix)
    if current_sheet:
        sheets.append(current_sheet)

    # 筛选真正的α螺旋（至少4个残基）
    true_helices = [h for h in helices if len(h) >= 4]
    # 筛选真正的β折叠（至少3个残基）
    true_sheets = [s for s in sheets if len(s) >= 3]

    # 输出结果
    print(f"PDB文件: {pdb_file}")
    print(f"链 {chain_id} 分析结果:")
    print(f"• α螺旋数量: {len(true_helices)}")
    print(f"• β折叠数量: {len(true_sheets)}")
    print(f"• 总残基数: {sum(len(h) for h in true_helices) + sum(len(s) for s in true_sheets)}")


#20251225填写
def calculate_nt_ct_distance(pdb_file, chain_id):
    """
    计算PDB文件中指定链的N端和C端碳原子之间的距离

    参数:
    pdb_file: PDB文件路径
    chain_id: 链标识符（如'A', 'B'等）

    返回:
    distance: N端和C端之间的距离（单位：Å）
    """

    # 创建解析器
    parser = Bio.PDB.PDBParser(QUIET=True)

    try:
        # 解析PDB文件
        structure = parser.get_structure('protein', pdb_file)

        # 获取指定链
        target_chain = None
        for model in structure:
            if chain_id in model:
                target_chain = model[chain_id]
                break

        if target_chain is None:
            raise ValueError(f"链 {chain_id} 在PDB文件中不存在")

        # 收集所有残基
        residues = list(target_chain.get_residues())

        if len(residues) < 2:
            raise ValueError("链中残基数不足，无法计算N端和C端距离")

        # 获取N端残基的N原子（第一个残基）
        n_term_residue = residues[0]
        n_atom = None

        # 尝试获取N原子（蛋白质N端通常有N原子）
        #if 'N' in n_term_residue:
        #    n_atom = n_term_residue['N']
        # 如果没有N原子，尝试获取其他代表性原子
        if 'CA' in n_term_residue:
            n_atom = n_term_residue['CA']
        elif len(n_term_residue) > 0:
            # 获取第一个原子
            n_atom = list(n_term_residue.get_atoms())[0]

        if n_atom is None:
            raise ValueError("无法找到N端残基的合适原子")

        # 获取C端残基的C原子（最后一个残基）
        c_term_residue = residues[-1]
        c_atom = None

        # 尝试获取C原子（蛋白质C端通常有C原子）

        # 如果没有C原子，尝试获取其他代表性原子
        if 'CA' in c_term_residue:
            c_atom = c_term_residue['CA']
        elif 'O' in c_term_residue:
            c_atom = c_term_residue['O']
        elif len(c_term_residue) > 0:
            # 获取最后一个原子
            c_atom = list(c_term_residue.get_atoms())[-1]

        if c_atom is None:
            raise ValueError("无法找到C端残基的合适原子")

        # 计算距离
        distance = np.linalg.norm(n_atom.coord - c_atom.coord)

        # 获取残基信息
        n_residue_info = f"{n_term_residue.resname}{n_term_residue.id[1]}"
        c_residue_info = f"{c_term_residue.resname}{c_term_residue.id[1]}"

        # 原子信息
        n_atom_info = n_atom.name
        c_atom_info = c_atom.name


        nt_ct_distance_info = {
            'distance': round(distance, 3),
            'n_term_residue': n_residue_info,
            'c_term_residue': c_residue_info,
            'n_atom': n_atom_info,
            'c_atom': c_atom_info,
            'chain_id': chain_id
        }
        return nt_ct_distance_info
        #distance,n_term_residue,c_term_residue,n_atom,c_atom,chain_id

        #nt_ct_distance_info = calculate_nt_ct_distance(pdb_file, chain_id)
    except Exception as e:
        return f"错误: {str(e)}"





if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python analyze_secondary.py <pdb文件>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    analyze_secondary_structure(pdb_file)



if __name__ == '__main__':
    #cal_interaction_by_pyrosetta(input_pdb=r'/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/20250721BSCTV_C3v4/helixfold3api_1/good_pdb/20250721BSCTV_C3v4_33_antibmpnn226_rank1.pdb')

    calculation_nanobody_rsmd(r'/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/20250719ToLCNDV_CP_corev8_shortcdr3/deepnano_out7/20250719ToLCNDV_CP_corev8_shortcdr3_37_dldesign_0_HT_antifold4_rank5_antibmpnn675_tfoldab.pdb',
                              r'/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/20250719ToLCNDV_CP_corev8_shortcdr3/deepnano_out7/20250719ToLCNDV_CP_corev8_shortcdr3_37_dldesign_0_HT_antifold4_rank5_antibmpnn865_tfoldab.pdb')

    calculation_nanobody_cdrs_rsmd(r'/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/20250719ToLCNDV_CP_corev8_shortcdr3/deepnano_out7/20250719ToLCNDV_CP_corev8_shortcdr3_26_dldesign_1_HT_antifold6_rank4_antibmpnn202_tfoldab.pdb',
                              r'/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/20250719ToLCNDV_CP_corev8_shortcdr3/deepnano_out7/20250719ToLCNDV_CP_corev8_shortcdr3_26_dldesign_1_HT_antifold6_rank4_antibmpnn202_tfoldab.pdb')

if __name__ == '__main__123':
    #a = My_nanobody()
    b = My_nanobody(sequence= r'MQVQLVESGGGLVQPGGSLRLSCAASGRTFSEYNMGWFRQAPGQGLEAVAAIRSSGTTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAMSRVDTDSPAFYDYWGQGTLVTVS')
    c = My_nanobody(
        sequence=r'QVQLQESGGGSVQVGGSLRVACAASGDTFSGYLAAWFRQAPGKGREGVAAINSKRHTTSYADSVKGRFTISKDNADNIMYLEMNSLKPEDTAIYYCAAADAIGLAEYWSTPTLSAARYKYWGQGTQVTVSS')
    d = My_nanobody(
        sequence=r'QVQLQESGGGSVQVGGSLRVACAASGDTFSGYLAAWFRQAPGKGREGVAAINSKRHTTSYADSVKGRFTISKDNADNIMYLEMNSLKPEDTAIYYCAAADAIGLAEYWSTPTLSAARYKYWGQGTQVTVSS',
    stringent=0)
    print(d.obtain_the_fixed_residues())




if __name__ == '__main__123':
    #a = My_nanobody()
    b = My_nanobody(sequence= r'MQVQLVESGGGLVQPGGSLRLSCAASGRTFSEYNMGWFRQAPGQGLEAVAAIRSSGTTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAMSRVDTDSPAFYDYWGQGTLVTVS')
    c = My_nanobody(
        sequence=r'QVQLQESGGGSVQVGGSLRVACAASGDTFSGYLAAWFRQAPGKGREGVAAINSKRHTTSYADSVKGRFTISKDNADNIMYLEMNSLKPEDTAIYYCAAADAIGLAEYWSTPTLSAARYKYWGQGTQVTVSS')
    d = My_nanobody(
        sequence=r'QVQLQESGGGSVQVGGSLRVACAASGDTFSGYLAAWFRQAPGKGREGVAAINSKRHTTSYADSVKGRFTISKDNADNIMYLEMNSLKPEDTAIYYCAAADAIGLAEYWSTPTLSAARYKYWGQGTQVTVSS',
    stringent=0)
    print(d.standized_nb())

    #get_sequence_from_fasta()
    #get_sequence_from_pdb()
    #get_sequence_from_nanobody_antigen_pdb()
    #standardize_100w_nanobody()

#IgGM
if __name__ == '__main__123':
    #IgGM分析进行输入文件准备。
    #get_sequence_from_nanobody_antigen_pdb()
    #input_pdb_or_cif = r"D:\my_script\Nanobody\IgGM\My_cases\mScarlet_20250614\input\fold_mscarlet_lam_8_model_0.cif"

    data_prepare_and_analyze_4_IgGM(input_pdb_or_cif = r"D:\my_script\Nanobody\IgGM\My_cases\mScarlet_20250616\input\fold_mscarlet_lam_8_model_0.cif",
                                    mut_cdr = 'mut_cdr123')



if __name__ == '__main__123':
    data_prepare_4_deepnano(antigen_name_list = ['ICMV_CP_core'],
                            nanobody_collection_fasta = r"D:\my_script\Nanobody\nano_collection\100w_nanobody_standardized.fasta")
if __name__ == '__main__123':

    data_analysis_from_deepnano_predict(deepnanp_output_file= r"D:\my_script\Nanobody\DeepNano_qinghua\My_cases\deepnano_input_ICMV_CP_core_etc_20250614\output\ICMV_CP_core_predictions1.csv")

#RFantibody 的结果处理
if __name__ == '__main__123':

    cif_to_pdb(r'D:\my_experi\structure_prediction\structure_effectors\Nanoviridae\fold_bbtv_cp')

    input_folder = r'D:\my_script\Nanobody\RFantibody\My_cases\BBTV_20250616b\rf2_outputs2'
    pdb_align(input_folder = input_folder,is_output_from_rf2=True)
    rf2antibody_filtering(input_folder = input_folder)

#GeoFlow的结果处理
if __name__ == '__main__123':

    input_folder = r'D:\my_script\Nanobody\RFantibody\My_cases\BBTV_20250616b\rf2_outputs2'
    pdb_align(input_folder = input_folder,is_output_from_rf2=False,is_nanobody_antigen=True,nanobody_standardized=True)


#cif_to_pdb
if __name__ == '__main__123':
    cif_to_pdb(r"D:\my_experi\structure_prediction\20250615_Nanobody_mScarlet\fold_mscarlet_lam_8")

# Acquire the sequences of 3 cdrs.
if __name__ == '__main__123':
    sequence = 'AQVQLVESGGGLVQAGGSLRLSCAASGFDFSDYTMNWFRQAPGKEREFVAAIRSSGATAYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAMSRVDTDDIYFFDYWGQGTLVTVSK'
    my_nanobody = My_nanobody(sequence,stringent=0)

#  Generating the Json file for AF3
if __name__ == '__main__123':
    os.chdir(r'D:\my_script\Nanobody')
    generating_json_as_af3_input()