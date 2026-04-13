import os
import argparse
import sys
import shutil
import pandas as pd
import numpy as np
import json
import SPIR_analysis
import time
from paddlehelix.task import helixfold3
#/ssd_diskB/gaolab/miniconda3/envs/phsdk/lib/python3.11/site-packages/paddlehelix/task/helixfold3.py


def parse_args2():
    parser = argparse.ArgumentParser(description=' To generate sequences into the input jsons of alphafold3.')
    parser.add_argument('--input_csv_file', '-icf', type=str, required=False,
                        default = '',
                        help='the input_csv_file containing antigen-nanobody information,\
                        There are at least three columns in this file: name, antigen, nanobody')
    parser.add_argument('--input_type', type=str, required=False,
                        default = 'csv',
                        help='csv, pdb or fasta')
    parser.add_argument('--output_folder', '-outf', type=str, required=False,
                        default = '/public_inspur_dss/Lab_share/GCX/haochengzhu/My_cases_copy/pvx_tgb1_20250616b/rf2_output_filtered_pdb',
                        help='the output_folder of collection of alphafold3 input json files')
    parser.add_argument('--modelSeeds', '-ms', type=int, required=False,
                        default = 1,
                        help='modelSeeds')
    parser.add_argument('--input_name_column', '-inc', type=str, required=False,
                        default = 'seq_new_name',
                        help='the name of column containing the input_name')
    parser.add_argument('--antigen_column', '-ac', type=str, required=False,
                        default = 'antigen_sequence',
                        help='the name of column containing the antigen sequence')
    parser.add_argument('--nanobody_column', '-nc', type=str, required=False,
                        default = 'nanobody_sequence_std',
                        help='the name of column containing the antigen sequence')
    parser.add_argument('--proportion2next_step', type=float, required=False,
                        default = 0.5,
                        help='Proportional of the nanobodies')
    parser.add_argument('--num_of_sequence2next_step', type=int, required=False,
                        default = 0,
                        help='num of the nanobodies')
    parser.add_argument('--iptm_cut_off', type=float, required=False,
                        default = 0,
                        help='iptm_cut_off')
    parser.add_argument('--only_filtering', type=int, required=False,
                        default = 1,
                        help='only_filtering or not')
    parser.add_argument('--summarize_results', type=int, required=False,
                        default = 1,
                        help='summarize the results or not')
    parser.add_argument('--num_of_structures_per_sequence', type=int, required=False,
                        default = 1,
                        help='num_of_structures_per_sequence')
    parser.add_argument('--max_job_num', type=int, required=False,
                        default = 3,
                        help='max_job_num')
    parser.add_argument('--execute_quiet', type=bool, required=False,
                        default = False,
                        help='execute_quiet')
    parser.add_argument('--binder_or_nanobody', type=str, required=False,
                        default = 'nanobody',
                        help='binder,nanobody,binder_dimer,binder_monomer')

    parser.add_argument('--helixfold_recycle', type=int, required=False,
                        default = 10,
                        help='helixfold_recycle')
    parser.add_argument('--helixfold_ensemble', type=int, required=False,
                        default = 1,
                        help='helixfold_ensemble')

    parser.add_argument('--combined_with_dimer', type=int, required=False,
                        default = 0,
                        help='combined_with_dimer')

    parser.add_argument('--hetero_folder', type=str, required=False,
                        default = '',
                        help='')

    parser.add_argument('--dimer_folder', type=str, required=False,
                        default = '',
                        help='')
    parser.add_argument('--monomer_folder', type=str, required=False,
                        default = '',
                        help='')
    parser.add_argument('--require_dna_seq', type=int, required=False,
                        default = 1,
                        help='')






    args = parser.parse_args()
    return args



def generating_json_as_helixfold3_input(sequences_list = ['',''],
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
    #直接会将所有的输入文件都写入json文件
    os.chdir(output_folder)
    # 处理单体输入（单条序列）
    output_file = input_name + '.json'
    if types == 'antigen_nanobody':
        data = {
            "job_name": input_name,
            "entities":[
            {
                "type": "protein",
                "sequence":sequences_list[0],
                "count": 1
            },
            {
                "type": "protein",
                "sequence": sequences_list[1],
                "count": 1
            }
        ]
        }


    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)  # 格式化写入json文件的数据
    #####################

    print(f"成功创建helixfold3输入文件: {output_file}")

def generating_helixfold3_input(input_csv_file = '',
                                input_name_column = '',
                                antigen_column = '',
                                nanobody_column = '',
                                output_folder = '',
                                device = '',
                                helixfold_recycle = 10,
                                helixfold_ensemble = 1,
                                binder_or_nanobody = ''):

    os.makedirs(f"{output_folder}/helixfold3_json_input",exist_ok=True)
    os.chdir(f"{output_folder}/helixfold3_json_input")
    input_csv = pd.read_csv(input_csv_file)
    input_csv = input_csv.dropna(subset=[antigen_column,nanobody_column,input_name_column])
    for each_line in input_csv.index:
        antigen_sequence = input_csv.loc[each_line,antigen_column]
        nanobody_sequence = input_csv.loc[each_line, nanobody_column]
        input_name = input_csv.loc[each_line, input_name_column]
        if binder_or_nanobody != 'binder_monomer':
            SPIR_analysis.generating_json_for_helixfold(sequences_list = [antigen_sequence, nanobody_sequence],
                                     input_name = input_name,types = 'antigen_nanobody',
                                     modelSeeds = 1,
                                     output_folder = f"{output_folder}/helixfold3_json_input",
                                      recycle = int(helixfold_recycle),
                                                           ensemble = int(helixfold_ensemble))
        if binder_or_nanobody == 'binder_monomer':
            SPIR_analysis.generating_json_for_helixfold(sequences_list = [nanobody_sequence],
                                     input_name = input_name,types = 'monomer',
                                     modelSeeds = 1,
                                     output_folder = f"{output_folder}/helixfold3_json_input",
                                      recycle = int(helixfold_recycle),
                                                           ensemble = int(helixfold_ensemble))

        '''
        generating_json_as_helixfold3_input(sequences_list=[antigen_sequence, nanobody_sequence],
                                            input_name=input_name, types='antigen_nanobody',
                                            modelSeeds=1,
                                            output_folder=f"{output_folder}/helixfold3_json_input")
        '''

#然后 写一下，对chai生成的序列进行summary
def helixfold3_api_summary(output_folder = '',proportion2next_step = 0.5,
                           num_of_sequence2next_step = 0,only_filtering = 0,
                           iptm_cut_off = 0,binder_or_nanobody = 'nanobody',require_dna_seq = 1):

    os.chdir(output_folder)
    summary_df = pd.DataFrame()
    print('开始进行helixfold3_api_summary')
    #print(output_folder)
    api_table_csv = pd.read_csv('table.csv')
    task_2_task_id = {}
    task_id_2_task = {}
    for each_index in api_table_csv.index:
        job_name_json = api_table_csv.loc[each_index,'data']
        if 'recycle' in job_name_json:
            task = job_name_json.split("'job_name': '")[1].split("', 'recycle':")[0]
        else:
            task = job_name_json.split("'job_name': '")[1].split("', 'entities':")[0]
        task_id = api_table_csv.loc[each_index,'task_id']
        task_2_task_id[task] = str(task_id)
        task_id_2_task[str(task_id)] = task
    print('task_2_task_id',task_2_task_id)
    print('task_id_2_task', task_id_2_task)

    #然后进入result文件夹
    os.chdir(f"{output_folder}/result")
    # 先把这个文件夹下面的所有文件夹都清空。
    if 'result_copy' in os.listdir(output_folder):
        shutil.rmtree(f"{output_folder}/result_copy")
    os.makedirs(f"{output_folder}/result_copy",exist_ok=True)
    #首先解压缩
    print('开始解压缩')
    print(os.listdir())
    for each_zip in os.listdir():
        if each_zip.endswith('.zip'):
            print('each_zip',each_zip)

            task_id = each_zip.split('_')[-2]
            if task_id in task_id_2_task.keys():
                task = task_id_2_task[task_id]
            else:
                task = each_zip.split('result_to_download_')[-1].rsplit(sep = '_',maxsplit=1)[0]
            new_zip_name = f"{task}.zip"

            shutil.copy(f"{output_folder}/result/{each_zip}",
                        f"{output_folder}/result_copy/{new_zip_name}")
            os.chdir(f"{output_folder}/result_copy/")
            os.system(f"unzip -o {output_folder}/result_copy/{new_zip_name}")
            print('unzip finished')
            #暂停4秒钟
            time.sleep(0.5)
            os.renames(f"{output_folder}/result_copy/{each_zip.split('.zip')[0]}",
                       f"{output_folder}/result_copy/{task}")
    #然后再改名




    #
    designed_dna_sequence_dict = {}
    for each_file in os.listdir(f"{output_folder}/result_copy/"):
        if not each_file.endswith('.zip'):
            os.chdir(f"{output_folder}/result_copy/{each_file}")#
            for each_file2 in os.listdir(f"{output_folder}/result_copy/{each_file}"):
                if '-rank' in each_file2:
                    print('each_file2',each_file2)
                    rank = each_file2.split('-rank')[-1]
                    SPIR_analysis.cif_to_pdb_v2(f"{output_folder}/result_copy/{each_file}/{each_file2}")
                    new_name = f"{each_file}_rank{rank}.pdb"
                    os.renames(f"{output_folder}/result_copy/{each_file}/{each_file2}/predicted_structure.pdb",
                               f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}")


                    shutil.copy(f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}", f"{output_folder}/result_copy/{each_file}/{new_name}")


                    if binder_or_nanobody == 'binder':
                        antigen_sequence1, binder_sequence1 = \
                            SPIR_analysis.get_sequence_from_antigen_binder_pdb(
                                f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}")
                        #prodigy
                        predicted_binding_affinity,predicted_dissociation_constant = SPIR_analysis.run_prodigy(
                            f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}", 'T', 'B')

                        ##计算binder的N末端 C末端之间的距离
                        nt_ct_distance_dict = SPIR_analysis.calculate_nt_ct_distance(pdb_file=f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}",
                                                                                         chain_id='B')
                        binder_nt_ct_distance = nt_ct_distance_dict['distance']

                        #foldx
                        os.chdir(f"{output_folder}/result_copy/{each_file}/{each_file2}")
                        foldx_Interaction_Energy, foldx_Interaction_Energy_opt, binder_interaction_residue, target_interaction_residue = 0, 0, 'A', 'A'
                        '''
                        foldx_Interaction_Energy,foldx_Interaction_Energy_opt,binder_interaction_residue,target_interaction_residue = \
                            SPIR_analysis.run_foldx_AnalyseComplex(
                            f"{new_name}", 'T', 'B')
                        '''



                    if binder_or_nanobody == 'binder_dimer':
                        binder_sequence1, binder_sequence2 = \
                            SPIR_analysis.get_sequence_from_binder_dimer_pdb(
                                f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}")
                        #prodigy
                        predicted_binding_affinity,predicted_dissociation_constant = SPIR_analysis.run_prodigy(
                            f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}", 'X', 'Y')


                        #foldx
                        os.chdir(f"{output_folder}/result_copy/{each_file}/{each_file2}")
                        foldx_Interaction_Energy, foldx_Interaction_Energy_opt, binder_interaction_residue, target_interaction_residue = 0, 0, 'A', 'A'
                        '''
                        foldx_Interaction_Energy,foldx_Interaction_Energy_opt,binder_interaction_residue,target_interaction_residue = \
                            SPIR_analysis.run_foldx_AnalyseComplex(
                            f"{new_name}", 'X', 'Y')
                        '''


                    if binder_or_nanobody == 'binder_monomer':
                        ##计算binder的N末端 C末端之间的距离
                        nt_ct_distance_dict = SPIR_analysis.calculate_nt_ct_distance(pdb_file=f"{output_folder}/result_copy/{each_file}/{each_file2}/{new_name}",
                                                                                         chain_id='A')
                        binder_nt_ct_distance = nt_ct_distance_dict['distance']


                    #然后进行汇总

                    with open(f"{output_folder}/result_copy/{each_file}/{each_file2}/all_results.json",
                              'r') as file:
                        json_data = json.load(file)
                        if binder_or_nanobody != 'binder_monomer':
                            iptm = json_data['iptm']
                        if binder_or_nanobody == 'binder_monomer':
                            iptm = json_data['iptm'] = 0

                        ptm = json_data['ptm']
                        ranking_confidence = json_data['ranking_confidence']
                        has_clash = json_data['has_clash']
                        mean_plddt = json_data['mean_plddt']

                    print(iptm, ptm, ranking_confidence, has_clash, mean_plddt)
                    df2append = pd.DataFrame()
                    df2append['results_name'] = f"{each_file}_rank{rank}".split()
                    df2append['iptm'] = iptm
                    df2append['ptm'] = ptm
                    df2append['ranking_confidence'] = ranking_confidence
                    df2append['has_clash'] = has_clash
                    df2append['mean_plddt'] = mean_plddt
                    df2append['iptm+ptm'] = iptm + ptm
                    df2append['case_name'] = each_file
                    if binder_or_nanobody != 'binder_dimer':
                        df2append['binder_nt_ct_distance'] = binder_nt_ct_distance

                    if binder_or_nanobody != 'binder_monomer':
                        df2append['Prodigy binding affinity (kcal.mol-1)'] = predicted_binding_affinity
                        df2append['Prodigy dissociation constant (M) at 25.0˚C'] = predicted_dissociation_constant
                        df2append['foldx_Interaction_Energy'] = foldx_Interaction_Energy
                        df2append['foldx_Interaction_Energy_opt'] = foldx_Interaction_Energy_opt
                        df2append['binder_interaction_residue'] = binder_interaction_residue
                        df2append['target_interaction_residue'] = target_interaction_residue




                    #print('prodigy_output',prodigy_output)



                    if binder_or_nanobody == 'binder':

                        df2append['binder_sequence'] = binder_sequence1
                        df2append['binder_length'] = len(binder_sequence1)
                        designed_sequence = binder_sequence1
                        df2append['antigen_sequence'] = antigen_sequence1

                    if binder_or_nanobody == 'binder_dimer':

                        df2append['binder_sequence1'] = binder_sequence1
                        df2append['binder_sequence2'] = binder_sequence2
                        df2append['binder_length'] = len(binder_sequence2)
                        designed_sequence = binder_sequence1
                        antigen_sequence1 = 'NANANAN'
                        df2append['antigen_sequence'] = antigen_sequence1.split()

                    if binder_or_nanobody != 'binder_monomer':
                        if require_dna_seq == 0:
                            designed_dna_sequence = 'AA'
                        elif designed_sequence in designed_dna_sequence_dict.keys():
                            designed_dna_sequence = designed_dna_sequence_dict[designed_sequence]
                            print(f'designed_dna_sequence already design: {designed_dna_sequence_dict[designed_sequence]}')
                        else:
                            designed_dna_sequence = SPIR_analysis.anchor_bar_generator \
                                (original_seq=designed_sequence,
                                 iteration=2000,
                                 GC_min=0.486, GC_max=0.605,
                                 bsa1_discard=True,
                                 organism='Oryza sativa',
                                 sort_by='cai',
                                 anchor_bar_distance=0,
                                 original_seq_type='prot')

                            designed_dna_sequence = SPIR_analysis.gc_content_upper(designed_dna_sequence)
                            designed_dna_sequence_dict[designed_sequence] = designed_dna_sequence

                        '''
                        antigen_dna_sequence = nanobody_analysis.anchor_bar_generator \
                            (original_seq=antigen_sequence1,
                             iteration=1000,
                             GC_min=0.48, GC_max=0.65,
                             bsa1_discard=True,
                             organism='Oryza sativa',
                             sort_by='cai',
                             anchor_bar_distance=0,
                             original_seq_type='prot')
                        '''

                        pikm_cc = 'MEAAAMAVTAATGALAPVLVKLAALLDDGECNLLEGSRSDAEFIRSELEAVHSLLTPNILGRMGDDDAACKDGLIAEVRELSYDLDDAVDDFLELNFEQRRSASPFGELKARVEERVSNRFSDWKLPAASLPPSSVHRRAGLPPPDAGLVGMHKRKEELIELLEQGSSDASRWRKRKPHVPLRG'
                        pikm_nbarc_lrr = 'KEITAMLAPVKSICEFHEVKTICILGLPGGGKTTIARVLYHALGTQFQCRVFASISPSSSPSPNLTETLADIFAQAQLGVTDTLSTPYGGSGTGRALQQHLIDNISAFLLNKKYLIVIDDIWHWEEWEVIRKSIPKNDLGGRIIMTTRLNSIAEKCHTDDNDVFVYEVGDLDNNDAWSLSWGIATKSGAGNRIGTGEDNSCYDIVNMCYGMPLALIWLSSALVGEIEELGGAEVKKCRDLRHIEDGILDIPSLQPLAESLCLGYNHLPLYLRTLLLYCSAYHWSNRIERGRLVRRWIAEGFVSEEKEAEGYFGELINRGWITQHGDNNSYNYYEIHPVMLAFLRCKSKEYNFLTCLGLGSDTSTSASSPRLIRRLSLQGGYPVDCLSSMSMDVSHTCSLVVLGDVARPKGIPFYMFKRLRVLDLEDNKDIQDSHLQGICEQLSLRVRYLGLKGTRIRKLPQEMRKLKHLEILYVGSTRISELPQEIGELKHLRILDVRNTDITELPLQIRELQHLHTLDVRNTPISELPPQVGKLQNLKIMCVRSTGVRELPKEIGELNHLQTLDVRNTRVRELPWQAGQISQSLRVLAGDSGDGVRLPEGVCEALINGIPGATRAKCREVLSIAIIDRFGPPLVGIFKVPGSHMRIPKMIKDHFRVLSCLDIRLCHKLEDDDQKFLAEMPNLQTLVLRFEALPRQPITINGTGFQMLESFRVDSRVPRIAFHEDAMPNLKLLEFKFYAGPASNDAIGITNLKSLQKVVFRCSPWYKSDAPGISATIDVVKKEAEEHPNRPITLLINAGYKEISTESHGSSENIAGSSGIDTEPAQAQHDNLPAVRDDYKGKGILLDGRCPTCGRATKIEEETQDRVADIEIQTETTS'
                        rga5_n = 'MSSSSLGAMDAPASFSLGAMGPLLRKLDSLLVAPEIRLPKPLKEGIELLKEDLEEIGVSLVEHSVVDSPTHKARFWMDEVRDLSYHIEDCIDTMFSMRSGGDDGKPRSERRHKVGRAKIDGFSKKPKPCTRMARIAELRALVREASERLERYQLGDVCGSSSPVVFTADGRARPLHHGVSANLVGVDEFKTKLNRWLSDEEGPHLKVAAIVGPAGIGKTALATELYRDHRWQFECRAFVRASRKPDMQRLLGGILSQVQRRQRSSDAYADSTVQSLIDNLREHLQDRRYLIIIDGLWETAVWNIANSAFPDVNSFSRILITADIEQVALECCGYKYDYIMRMEPLGSLDSKKVFFNKVFGSEDQCPPELKEVSNTILEKCGGLPLAIISIAGLLGSQPENPVLWDYVTKYLCSSLGTNPTLKDVVKETLNLSYNSLPHPFKTCLLYLGMYPDGHIMLKADLMKQWSAEGFVSANEAKDTEEIVDKYFDELVNRGILEPVEINKNGKVLSCTLHHAVHDLVMPKFNDDKFTMSVDYSQTITGPSTMVRRLSLHFSSTRYATKPAGIILSRVRSLAFFGLLNCMPCIGEFKLLRVLILEFWGSHGEQRSLNLIPVCRLFQLRYLKTSGDVVVQLPAQISGLQYLETLEIDARVSAVPFDLVHLPNLLHLQLQDETKLPDGIGCMRSLRTLQYFDLGNNSVDNLRGLGELTNLQDLHLSYSAPSSNEGLMINLNAITSSLSRLSNLKSLILSPGAISMVIFFDISSIISVVPVFLQRLELLPPICIFCRLPKSIGQLHKLCILKVSVRELLTTDIDNLTGLPSLTVLSLYAQTAPEGRFIFKDGTLPVLKYFKFGCGELCLAFMAGAMPNLQRLKLVFNIRKSEKYRHTLFGIEHLVSLQDIATRIGVDTSTGESDRRAAESAFKETVNKHPRCLRSSLQWVVSTEEESHPLEKQHHKREKGSSAGHGVLEKESVEDSEKNTDRVQTLLSPQLSNMESVVESAL'
                        rga5_c = 'LAGGKKGAYKKHPTYNLSPFDYVEYPPSAPIMQDINPCSTM'


                        df2append['designed_sequence'] = designed_sequence
                        df2append['designed_in_pikm'] = f'{pikm_cc}{designed_sequence}{pikm_nbarc_lrr}'
                        df2append['designed_in_rga5'] = f'{rga5_n}{designed_sequence}{rga5_c}'
                        df2append['designed_dna_sequence'] = designed_dna_sequence
                        #df2append['antigen_dna_sequence'] = antigen_dna_sequence



                    summary_df = summary_df._append(df2append)
                    summary_df = summary_df.sort_values(by='iptm+ptm', ascending=False)
                    summary_df.index = np.arange(len(summary_df))
                    summary_df.to_csv(f"{output_folder}/helixfold3_api_summary.csv")

            if ('merge.pse' not in os.listdir(f"{output_folder}/result_copy/{each_file}") or only_filtering == 1) \
                    and binder_or_nanobody != 'binder_monomer':
                SPIR_analysis.cif_to_pdb(input_folder=f"{output_folder}/result_copy/{each_file}")
                SPIR_analysis.extra_alignment(input_folder=f"{output_folder}/result_copy/{each_file}",only_align=True)

                #然后再重新将pdb文件当中的chain换成HT(nanobody)或者TB (binder)
                for each_pdb2 in os.listdir(f"{output_folder}/result_copy/{each_file}"):
                    if each_pdb2.endswith('.pdb'):
                        if binder_or_nanobody == 'nanobody':
                            antigen_sequence3, nanobody_sequence3 = \
                                SPIR_analysis.get_sequence_from_nanobody_antigen_pdb(
                                    f"{output_folder}/result_copy/{each_file}/{each_pdb2}")
                        if binder_or_nanobody == 'binder':
                            antigen_sequence3, binder_sequence3 = \
                                SPIR_analysis.get_sequence_from_antigen_binder_pdb(
                                    f"{output_folder}/result_copy/{each_file}/{each_pdb2}")
                        if binder_or_nanobody == 'binder_dimer':
                            _, _ = \
                                SPIR_analysis.get_sequence_from_binder_dimer_pdb(
                                    f"{output_folder}/result_copy/{each_file}/{each_pdb2}")

    if binder_or_nanobody == 'binder_monomer':
        #到这里就可以结束了。
        return 'good, we have finished'
    #然后整理文件
    nanobody_list = list(set(summary_df['designed_sequence']))
    summary_df_best = pd.DataFrame()
    for each_nb_sequence in nanobody_list:
        sub_df = summary_df.loc[summary_df['designed_sequence'] == each_nb_sequence, :].copy()
        sub_df = sub_df.sort_values(by='iptm+ptm', ascending=False)
        sub_df.index = np.arange(len(sub_df))
        summary_df_best = summary_df_best._append(sub_df.loc[0, :])

    summary_df_best.sort_values(by='ranking_confidence', ascending=False, inplace=True)
    summary_df_best.index = np.arange(len(summary_df_best))
    summary_df_best.to_csv(f"{output_folder}/helixfold3_summary_best.csv")
    summary_df.to_csv(f"{output_folder}/helixfold3_api_summary.csv")
    os.chdir(output_folder)

    # 接下来进行过滤
    if num_of_sequence2next_step == 0:
        #此时利用proportion2next_step分析
        if proportion2next_step == 0:
            print('num_of_sequence2next_step and proportion2next_step can not be both 0')
            num_of_sequence2next_step = 10
        else:
            num_of_sequence2next_step = round(len(summary_df_best)*proportion2next_step)

    summary_df['rank'] = np.arange(len(summary_df))+1
    summary_filtered_df = summary_df.loc[summary_df['rank']<= num_of_sequence2next_step,:]
    summary_filtered_df.to_csv(f"{output_folder}/helixfold3_summary_filtered.csv")

    if iptm_cut_off != 0:
        summary_filtered_df = summary_df.loc[summary_df['iptm'] >= iptm_cut_off, :]
        summary_filtered_df.to_csv(f"{output_folder}/helixfold3_summary_filtered.csv")

    #然后把好的部分挑出来，放到一个新的文件夹里面。
    summary_filtered_df = pd.read_csv(f"{output_folder}/helixfold3_summary_filtered.csv",index_col=0)
    if 'good_pdb' in os.listdir(output_folder):
        shutil.rmtree(f"{output_folder}/good_pdb")
    os.makedirs(f"{output_folder}/good_pdb", exist_ok=True)
    for each_index in summary_filtered_df.index:
        result_name = summary_filtered_df.loc[each_index,"results_name"]
        for each_file3 in os.listdir(f"{output_folder}/result_copy/"):
            if not each_file3.endswith('.zip'):
                for each_pdb in os.listdir(f"{output_folder}/result_copy/{each_file3}"):
                    if each_pdb == f"{result_name}.pdb":
                        shutil.copy(f"{output_folder}/result_copy/{each_file3}/{result_name}.pdb",
                                    f"{output_folder}/good_pdb/{result_name}.pdb")

    SPIR_analysis.extra_alignment(f"{output_folder}/good_pdb",only_align = True)
    #把best也复制一份
    os.makedirs(f"{output_folder}/best_good_pdb", exist_ok=True)
    for each_index in summary_df_best.index:
        result_name = summary_df_best.loc[each_index,"results_name"]
        for each_file3 in os.listdir(f"{output_folder}/result_copy/"):
            if not each_file3.endswith('.zip'):
                for each_pdb in os.listdir(f"{output_folder}/result_copy/{each_file3}"):
                    if each_pdb == f"{result_name}.pdb":
                        shutil.copy(f"{output_folder}/result_copy/{each_file3}/{result_name}.pdb",
                                    f"{output_folder}/best_good_pdb/{result_name}.pdb")






    #接下来，准备计算



#下一步，运行helixfold3
def run_helixfold3_api_in_folder(input_csv_file = '',
                             output_folder='', modelSeeds=0, only_filtering=0,summarize_results = 0,
                                 num_of_structures_per_sequence=1,max_job_num = 2,
                                 execute_quiet = False,num_of_sequence2next_step = 10,
                                 iptm_cut_off = 0,binder_or_nanobody = 'nanobody',require_dna_seq = 1):

    input_csv = pd.read_csv(input_csv_file,index_col=0)
    input_csv.index = np.arange(len(input_csv))
    #print(input_csv)
    os.makedirs(output_folder, exist_ok=True)

    if only_filtering == 1:
        input_csv = input_csv.loc[0:0,:]
        print(input_csv)

    #进入json文件夹
    os.chdir(f"{output_folder}/helixfold3_json_input")
    json_list = []
    #print(os.listdir())
    for each_json in os.listdir():
        if each_json.endswith('.json'):
            json_list.append(f"{output_folder}/helixfold3_json_input/{each_json}")
    #print('json_list: ',json_list)

    '''
    for each_line in input_csv.index:
        antigen_sequence = input_csv.loc[each_line,antigen_column]
        nanobody_sequence = input_csv.loc[each_line, nanobody_column]
        input_name = input_csv.loc[each_line, input_name_column]
    '''
    json_list_to_submit = []
    num_of_job = 0
    for each_json in json_list:
        #然后检查该任务是否已经存在
        input_name = each_json.split(sep='/')[-1].split('.json')[0]
        task_exist = 0
        for each_file0 in os.listdir(output_folder):
            if input_name in each_file0:
                print(f"{input_name}任务已存在，跳过该任务。")
                task_exist = 1
        if task_exist == 0 and num_of_job < max_job_num:
            json_list_to_submit.append(each_json)
            num_of_job += 1
            print(f'准备进行{input_name}的计算')
    #print('json_list_to_submit: ', json_list_to_submit)

    #如果这个文件夹存在，则全部删除之。
    if os.path.exists(f"{output_folder}/helixfold3_json_input4submit"):
        shutil.rmtree(f"{output_folder}/helixfold3_json_input4submit")
    os.makedirs(f"{output_folder}/helixfold3_json_input4submit",exist_ok=False)

    for each_json2 in json_list_to_submit:
        json_name = each_json2.split('/')[-1]
        shutil.copy(each_json2,f"{output_folder}/helixfold3_json_input4submit/{json_name}")

    #20250912
    #检查helixfold3_json_input4submit中的项目是否已经存在于.result文件夹中
    if 'result' in os.listdir(output_folder):
        for each_json3 in os.listdir(f"{output_folder}/helixfold3_json_input4submit"):
            json_name = each_json3.split('.json')[0]
            for each_result in os.listdir(f"{output_folder}/result"):
                if json_name+'_' in each_result:
                    print('results already existed, removed task '+each_json3)
                    os.remove(f"{output_folder}/helixfold3_json_input4submit/{each_json3}")

    go_on = 'no'
    for j in range(num_of_structures_per_sequence):
        if only_filtering == 0:
            #print(input_name, antigen_sequence, nanobody_sequence)

            helixfold3.execute(input_data=f"{output_folder}/helixfold3_json_input4submit/", output_dir=f"{output_folder}",quiet = False)
            #print("10秒后继续执行")
            #time.sleep(10)

            go_on = input('print_yes:')

    go_on = 'yes'
    if go_on == 'yes' or only_filtering == 1:
        if summarize_results == 1:

            helixfold3_api_summary(output_folder=output_folder, proportion2next_step=proportion2next_step,
                                   num_of_sequence2next_step=num_of_sequence2next_step, only_filtering=only_filtering,
                                   iptm_cut_off = iptm_cut_off,binder_or_nanobody = binder_or_nanobody,require_dna_seq = require_dna_seq)
    '''
    helixfold3_api_summary(output_folder=output_folder, proportion2next_step=proportion2next_step,
                       num_of_sequence2next_step=num_of_sequence2next_step, only_filtering=only_filtering,
                           iptm_cut_off = iptm_cut_off)
    '''

    #print(rf"congratulations!,已完成{num_of_job_finished/num_of_job*100}%的任务")

def combined_with_dimer_results(hetero_folder = '', dimer_folder = '',monomer_folder = '',bindcraft_summary_csv = 'bindcraft_summary.csv'):
    hetero_summary_df = pd.read_csv(f"{hetero_folder}/helixfold3_api_summary.csv",index_col=0)
    hetero_summary_df.index = hetero_summary_df['results_name']
    hetero_summary_df['hetero_dimer_iptm'] = hetero_summary_df['iptm']
    hetero_summary_df['hetero_dimer_ranking_confidence'] = hetero_summary_df['ranking_confidence']

    if os.path.exists(dimer_folder) and dimer_folder != '':
        dimer_summary_df = pd.read_csv(f"{dimer_folder}/helixfold3_api_summary.csv",index_col=0)
        dimer_summary_df.index = dimer_summary_df['results_name']
        hetero_summary_df['homo_dimer_iptm'] = dimer_summary_df['iptm']
        hetero_summary_df['homo_dimer_ranking_confidence'] = dimer_summary_df['ranking_confidence']
        hetero_summary_df['homo_dimer_foldx_Interaction_Energy_opt'] = dimer_summary_df['foldx_Interaction_Energy_opt']

    if os.path.exists(monomer_folder) and monomer_folder != '':
        monomer_summary_df = pd.read_csv(f"{monomer_folder}/helixfold3_api_summary.csv", index_col=0)
        monomer_summary_df.index = monomer_summary_df['results_name']
        hetero_summary_df['monomer_ptm'] = monomer_summary_df['ptm']
        hetero_summary_df['monomer_ranking_confidence'] = monomer_summary_df['ranking_confidence']
        hetero_summary_df['homo_dimer_iptm+monomer_ptm'] = hetero_summary_df['monomer_ptm'] + hetero_summary_df['hetero_dimer_iptm']
        hetero_summary_df['monomer_binder_nt_ct_distance'] = monomer_summary_df['binder_nt_ct_distance']


        #然后将monomer的结果和hetero dimer的结果进行merge
        merge_pse_path = f"{hetero_folder}/merge_pse"
        os.makedirs(merge_pse_path,exist_ok=True)

        for each_folder1 in os.listdir(f"{hetero_folder}/result_copy/"):
            if each_folder1.endswith('.zip'):
                continue
                #是压缩包，忽略过去
            case_name = each_folder1
            print('case_name: ',case_name)
            for each_pdb1 in os.listdir(f"{hetero_folder}/result_copy/{each_folder1}"):
                if not each_pdb1.endswith('rank1.pdb'):
                    continue
                align_pdb_path1 = f"{hetero_folder}/result_copy/{each_folder1}/{each_pdb1}"

                #接下里去monomer_folder里面找
                for each_folder2 in os.listdir(f"{monomer_folder}/result_copy/"):
                    if each_folder2 == case_name and not each_folder2.endswith('.zip'):
                        for each_pdb2 in os.listdir(f"{monomer_folder}/result_copy/{each_folder2}"):
                            if not each_pdb2.endswith('rank1.pdb'):
                                continue
                            align_pdb_path2 = f"{monomer_folder}/result_copy/{each_folder2}/{each_pdb2}"

                            SPIR_analysis.simple_alignment(align_pdb_path1=align_pdb_path2,
                                             align_pdb_path2=align_pdb_path1,
                                             output_pse_name=f"{hetero_folder}/merge_pse/{case_name}_merge.pse")





    hetero_summary_df.sort_values(by='ranking_confidence', ascending=False,inplace=True)
    hetero_summary_df.index = np.arange(len(hetero_summary_df))
    hetero_summary_df.to_csv(f"{hetero_folder}/helixfold3_with_dimer_summary.csv")

    #然后整理best
    best_df = pd.DataFrame()
    for each_case_name in list(set(hetero_summary_df['case_name'])):
        sub_df = hetero_summary_df.loc[hetero_summary_df['case_name'] == each_case_name,:]
        sub_df.index = np.arange(len(sub_df))
        #平均一下binder_nt_ct_distance
        sub_df['average_nt_ct_distance'] = np.average(sub_df['binder_nt_ct_distance'])
        best_df = best_df._append(sub_df.loc[0,:])

    best_df.sort_values(by='ranking_confidence', ascending=False,inplace=True)
    best_df.index = np.arange(len(best_df))
    best_df.to_csv(f"{hetero_folder}/helixfold3_bestone_with_dimer_summary.csv")

    #然后将文件与bindcraft生成的‘final_design_stats.csv’文件中的信息进行合并
    column_need_to_merge = ['MPNN_score','MPNN_seq_recovery',
                            'Average_pLDDT','Average_pTM',
                            'Average_i_pTM','Average_pAE',
                            'Average_i_pAE','Average_i_pLDDT',
                            'Average_ss_pLDDT','Average_Unrelaxed_Clashes',
                            'Average_Relaxed_Clashes','Average_Binder_Energy_Score',
                            'Average_Surface_Hydrophobicity','Average_ShapeComplementarity',
                            'Average_PackStat','Average_dG',
                            'Average_dSASA','Average_dG/dSASA',
                            'Average_Interface_SASA_%','Average_Interface_Hydrophobicity',
                            'Average_n_InterfaceResidues','Average_n_InterfaceHbonds',
                            'Average_InterfaceHbondsPercentage','Average_n_InterfaceUnsatHbonds',
                            'Average_InterfaceUnsatHbondsPercentage','Average_Interface_Helix%',
                            'Average_Interface_BetaSheet%','Average_Interface_Loop%',
                            'Average_Binder_Helix%','Average_Binder_BetaSheet%',
                            'Average_Binder_Loop%','Average_Hotspot_RMSD',
                            'Average_Target_RMSD','Average_Binder_pLDDT',
                            'Average_Binder_pTM','Average_Binder_pAE',
                            'Average_Binder_RMSD','Average_Binder_BetaSheet%',]

    if os.path.exists(bindcraft_summary_csv): #如果存在bindcraft_summary_csv，则执行下面的命令。
        bindcraft_summary_df = pd.read_csv(bindcraft_summary_csv,index_col=0)
        bindcraft_summary_df.index = bindcraft_summary_df['Design']

        for each_index in best_df.index:
            #print('each_index:',each_index)
            case_name = best_df.loc[each_index,'case_name']
            print('case_name:', case_name)

            if case_name in bindcraft_summary_df.index:
                for each_item in column_need_to_merge:
                    best_df.loc[each_index, each_item] = bindcraft_summary_df.loc[case_name,each_item]

        best_df.to_csv(f"{hetero_folder}/helixfold3_bestone_with_dimer_summary.csv")








if __name__ == '__main__':
    args = parse_args2()
    input_csv_file = args.input_csv_file
    input_type = args.input_type
    output_folder = args.output_folder

    modelSeeds = args.modelSeeds
    input_name_column = args.input_name_column
    antigen_column = args.antigen_column
    nanobody_column = args.nanobody_column

    proportion2next_step = args.proportion2next_step
    num_of_sequence2next_step = args.num_of_sequence2next_step
    only_filtering = args.only_filtering
    summarize_results = args.summarize_results
    num_of_structures_per_sequence = args.num_of_structures_per_sequence
    execute_quiet = args.execute_quiet
    iptm_cut_off = args.iptm_cut_off
    binder_or_nanobody = args.binder_or_nanobody
    helixfold_recycle = args.helixfold_recycle
    helixfold_ensemble = args.helixfold_ensemble
    max_job_num = args.max_job_num

    combined_with_dimer = args.combined_with_dimer
    hetero_folder = args.hetero_folder
    dimer_folder = args.dimer_folder
    monomer_folder = args.monomer_folder
    require_dna_seq = args.require_dna_seq

    #nanobody_analysis.extra_alignment(f"{output_folder}/best_good_pdb")
    #'''
    if combined_with_dimer == 0:
        if only_filtering == 0:
            generating_helixfold3_input(input_csv_file = input_csv_file,
                                        binder_or_nanobody = binder_or_nanobody,
                                  input_name_column = input_name_column,
                                  antigen_column = antigen_column,
                                  nanobody_column = nanobody_column,
                                  output_folder = output_folder,
                                        helixfold_recycle = helixfold_recycle,
                                        helixfold_ensemble = helixfold_ensemble)


        run_helixfold3_api_in_folder(input_csv_file=input_csv_file,
                                     output_folder=output_folder, modelSeeds=modelSeeds, only_filtering=only_filtering,
                                     summarize_results = summarize_results,require_dna_seq = require_dna_seq,
                                     num_of_structures_per_sequence=1, max_job_num=max_job_num,
                                     execute_quiet = execute_quiet,num_of_sequence2next_step = num_of_sequence2next_step,
                                     iptm_cut_off = iptm_cut_off,binder_or_nanobody = binder_or_nanobody)

    if combined_with_dimer ==1:
        combined_with_dimer_results(hetero_folder = hetero_folder, dimer_folder = dimer_folder,monomer_folder = monomer_folder)

    #'''
    if binder_or_nanobody == 'binder' and require_dna_seq == 1 and only_filtering == 1 and summarize_results == 1:
        hetero_folder = output_folder.rsplit('//', maxsplit=1)[0] + '/HF3_binder_effector_dimer '
        dimer_folder = output_folder.rsplit('//', maxsplit=1)[0] + '/HF3_binder_dimer'
        monomer_folder = output_folder.rsplit('//', maxsplit=1)[0] + '/HF3_binder_monomer '
        '''
        hetero_folder = output_folder[:-1]+'1'
        dimer_folder = output_folder[:-1]+'2'
        monomer_folder = output_folder[:-1]+'3'
        '''
        bindcraft_summary_csv = output_folder+'/../bindcraft_summary.csv'

        combined_with_dimer_results(hetero_folder=hetero_folder, dimer_folder=dimer_folder,
                                    monomer_folder=monomer_folder, bindcraft_summary_csv = bindcraft_summary_csv)





