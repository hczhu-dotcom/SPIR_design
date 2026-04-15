import os
import argparse
import sys
import shutil
import numpy as np
import pandas as pd
import json
import SPIR_analysis


# 这个脚本是根据一个csv中的antigen-nanobody信息，生成用于run af3的input序列
# 允许两种输入模式，一个是csv文件


def parse_args1():
    parser = argparse.ArgumentParser(description=' To generate sequences into the input jsons of alphafold3.')
    parser.add_argument('--case_name', type=str, required=False,
                        default = 'PVY_CPv1',
                        help='')
    parser.add_argument('--pathogen_protein_path', type=str, required=False,
                        default = "/public-supool/home/gaolab/haochengzhu/SPIR_design/Pathogen_protein/PVY_CP.pdb",
                        help='')
    parser.add_argument('--binder_lengths', type=str, required=False,
                        default = '[69, 117]',
                        help='')
    parser.add_argument('--number_of_final_designs', type=int, required=False,
                        default = 12,
                        help='')
    parser.add_argument('--target_hotspot_residues', type=str, required=False,
                        default = 'null',
                        help='')


    args = parser.parse_args()
    return args



def candidate_filtering(args):
    case_name = args.case_name
    pathogen_protein_path = args.pathogen_protein_path
    binder_lengths = args.binder_lengths
    number_of_final_designs = args.number_of_final_designs
    target_hotspot_residues = args.target_hotspot_residues

    binder_lengths_bottom = int(binder_lengths[1:].split(',')[0])
    binder_lengths_up = int(binder_lengths[:-1].split(',')[1])
    print('binder_lengths_bottom',binder_lengths_bottom)
    print('binder_lengths_up', binder_lengths_up)

    #pathogen_protein = pathogen_protein_path.split('//')[-1].split('.pdb')[0]
    pathogen_protein = os.path.split(pathogen_protein_path)[-1].split('.pdb')[0]
    print('pathogen_protein:',pathogen_protein)

    if pathogen_protein_path == '':
        pathogen_protein = case_name[8:].rsplit(sep = 'v', maxsplit=1)[0]
        print(f'未提供pathogen_protein信息,默认pathogen_protein为{pathogen_protein}')

    if target_hotspot_residues == 'null':
        target_hotspot = None
    else:
        #131-141
        target_hotpot_low = int(target_hotspot_residues.split('-')[0])
        target_hotpot_up = int(target_hotspot_residues.split('-')[-1])
        target_hotspot = target_hotspot_residues
        pass

    current_directory = os.getcwd()
    data = {
    'design_path': f"{current_directory}/../BindCraft/my_cases/{case_name}/",
    "binder_name": f"{pathogen_protein}_Bd",
    "starting_pdb": f"{pathogen_protein_path}",
    "chains": "A",
    "target_hotspot_residues": target_hotspot,
    "lengths": [binder_lengths_bottom, binder_lengths_up],
    "number_of_final_designs": number_of_final_designs
    }

    os.chdir(r"../BindCraft/settings_target/")
    with open(f"{case_name}.json", 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据

    print(f'Successfully finished the json file of {case_name}')


def generating_effector_fasta_collection(args):
    # Generating a fasta format file which contains all sequence of pathogen effector.
    pathogen_protein_path = args.pathogen_protein_path
    #print(pathogen_protein_path)
    pathogen_protein_path_folder,_ = os.path.split(pathogen_protein_path)
    #os.chdir(pathogen_protein_path_folder)
    #read all pdb file in this folder
    #checking the fasta collection
    #if 'pathogen_protein.fasta' in os.listdir('..'):
    #os.chdir('..')
    out_file = open(f'{pathogen_protein_path_folder}/../pathogen_protein.fasta',mode='w')

    for each_file in os.listdir(pathogen_protein_path_folder):
        if each_file.endswith('.pdb'):
            #this is a pdb file
            aa,residue_numbers = SPIR_analysis.get_the_sequence_from_single_pdb(f"{pathogen_protein_path_folder}/{each_file}", chain_id = 'A',silent = 1)
            name = each_file.split('.pdb')[0]
            out_file.write('>' + name + '\n')
            out_file.write(aa + '\n')
    out_file.close()











if __name__ == '__main__':
    args = parse_args1()
    candidate_filtering(args)
    generating_effector_fasta_collection(args)

