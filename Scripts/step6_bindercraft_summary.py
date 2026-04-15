import os
import argparse
import sys
import shutil
import numpy as np
import pandas as pd
import SPIR_analysis
#import beeswarm
import plotnine
#import PyRosetta_analysis



# 这个脚本是根据一个csv中的antigen-nanobody信息，生成用于run af3的input序列
# 允许两种输入模式，一个是csv文件


def parse_args2():
    parser = argparse.ArgumentParser(description=' To generate sequences into the input jsons of alphafold3.')
    parser.add_argument('--input_folder', type=str, required=True,
                        default = '',
                        help='')
    parser.add_argument('--output_csv', type=str, required=True,
                        default = '',
                        help='')

    args = parser.parse_args()
    return args



def candidate_filtering(args):
    input_folder = args.input_folder
    output_csv = args.output_csv

    '''
    # identiy the absolute path
    current_path = os.getcwd()
    if "SPIR_design" not in current_path:
        print('please check the current path')
    if "SPIR_design" in current_path:
        base_path = current_path.split('SPIR_design')[0]
    if "BindCraft" in current_path:
        base_path = current_path.split('BindCraft')[0]
    if "BindCraft" in input_folder and "BindCraft" in output_csv:
        input_folder = base_path+"/BindCraft/"+input_folder.split('BindCraft')[1]
        output_csv = base_path+"/BindCraft/"+output_csv.split('BindCraft')[1]
    '''


    if input_folder != '':
        #os.chdir(input_folder)
    #读取final_design_stats.csv

        output_df = pd.read_csv(f"{input_folder}/final_design_stats.csv", index_col=0)

    output_df.index = np.arange(len(output_df))

    #首先更改pdb文件中的链名
    for each_pdb in os.listdir(f'{input_folder}/Accepted'):
        if not each_pdb.endswith('.pdb'):
            continue
        pdb_path = f'{input_folder}/Accepted/{each_pdb}'
        case_name = each_pdb.split('.pdb')[0]
        antigen_sequence, binder_sequence = SPIR_analysis.get_sequence_from_antigen_binder_pdb(pdb_path = pdb_path)
        #获取信息的同时 change chain name, target_name = 'T',binder_name = 'B'

        #计算binder的N末端 C末端之间的距离
        nt_ct_distance_dict = SPIR_analysis.calculate_nt_ct_distance(pdb_file = pdb_path, chain_id = 'B')
        nt_ct_distance = nt_ct_distance_dict['distance']


        print('pdb文件是： ',each_pdb)
        print('case_name ',case_name)
        print('antigen_sequence: ',antigen_sequence)
        print('binder_sequence: ', binder_sequence)
        print('nt_ct_distance: ', nt_ct_distance)

        for each_line in output_df.index:
            if output_df.loc[each_line,'Design'] in case_name:
                output_df.loc[each_line, 'antigen_sequence'] = antigen_sequence
                output_df.loc[each_line, 'binder_sequence'] = binder_sequence
                output_df.loc[each_line, 'nt_ct_distance'] = nt_ct_distance


    #去重
    print(f'去重复前，binder个数为:{len(output_df)}个')
    output_df = output_df.drop_duplicates(subset=['binder_sequence'])
    print(f'去重复后，binder个数为:{len(output_df)}个')
    output_df = output_df.sort_values(by = 'Average_i_pTM',ascending=False)
    output_df.index = np.arange(len(output_df))
    output_df.to_csv(output_csv)







if __name__ == '__main__':
    args = parse_args2()
    candidate_filtering(args)

