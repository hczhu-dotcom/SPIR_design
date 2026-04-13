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

current_path = os.getcwd()

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

def get_the_sequence_from_single_pdb(input_pdb = r'',
                                     chain_id = 'A',silent = 0):
    '''
    这个脚本是将一个PDF文件，中的一个chain，转化成为序列
    需要同时指定PDF文件的路径和链名称。
    :param input_pdb:
    :param chain_id:
    :return:
    '''

    if not input_pdb.endswith('.pdb'):
        print('it seems that this is not a pdb-format file.')
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
                                    chain_id = 'all',silent = 1,
                                    out_file_name = ''):
    '''
    这个程序是将一个文件夹下面所有的pdb文件，全部提取氨基酸序列，并输出成为fasta文件和csv文件。
    :param input_pdb_folder:
    :param chain_id: 指定要提取序列的chain_id，默认为all，即提取所有的序列。
    :param silent:
    :return:
    '''
    out_file = open(out_file_name, mode='w')
    for each_file in os.listdir(input_pdb_folder):
        if each_file.endswith('.pdb'):
            #this is a pdb file
            aa,residue_numbers = get_the_sequence_from_single_pdb(f"{input_pdb_folder}/{each_file}", chain_id = 'A',silent = silent)
            name = each_file.split('.pdb')[0]
            out_file.write('>' + name + '\n')
            out_file.write(aa + '\n')
    out_file.close()
    pass


def get_sequence_from_pdb(pdb_path = r'D:\my_experi\structure_prediction\20250515_Nanobody_CMV\CMV_CP_core',
                          out_fasta_file = 'pathogen.fasta',
                          chain_required = ''):
    '''
    这个程序是用pymol程序包写的，其实不是那么的友好。
    #这个脚本是将一个目录下所有的PDF文件，转化成为fasta格式的序列信息，并输出成为文件
    # 定义三字母到单字母的映射
    '''

    cif_to_pdb(pdb_path)
    os.chdir(pdb_path)
    #确定一个输出文件
    out_file = open(out_fasta_file,mode='w')
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
                                         antigen_collection_fasta = f"{current_path}/pathogen_protein.fasta"):
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

def get_sequence_from_binder_dimer_pdb(pdb_path = r'',
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
    if len(chain_sequences[1]) != len(chain_sequences[0]):
        print("the sequence of the two chains are not the same")


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
                   Codon_Usage_Tables_folder=f'{current_path}/Codon_Usage_Tables'):
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
        #df2append['cai'] = 1
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



if __name__ == "__main__123":
    if len(sys.argv) != 2:
        print("用法: python analyze_secondary.py <pdb文件>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    analyze_secondary_structure(pdb_file)


