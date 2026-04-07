# SPIR_design
Code for SPIR design

This repository is used for designing of Synthetic Plant Immune Receptor (SPIR)

This repository is based on BindCraft (https://github.com/martinpacesa/BindCraft) and Helixfold 3 (https://gitee.com/paddlehelix/paddlehelix), so please Install these two software first, and create two conda environments (named BindCraft and phsdk, respectively, if use the default names)

## Step1: Cloning this repository.
```bash
git clone https://github.com/hczhu-dotcom/SPIR_design
```

## Step2: Predecting the structure of plant pathogen protein 
Predecting the structure of plant pathogen protein using AlphaFold 3, by using its webserver: https://alphafoldserver.com/.

## Step3: Triming the pathogen protein
Triming the N and C unstructured region of the plant pathogen protein using PyMol, and converting the CIF (Crystallographic Information File) format into the PDB (Protein Data Bank) format. This process would result the pathogen_protein.pdb, such as PVY_CP.pdb.

## Step4: Generating the input json file for BindCraft, by running: 

```bash
conda activate BindCraft
python step4_bindercraft_input_json \
	--your_case_name #the_name_of_your_job, such as PVY_CP \
--pathogen_protein_path #path_of_your_pathogen_protein_structure, such as Pathogen_protein/PVY_CP.pdb \
	--binder_lengths the_max_and_min_of_binder # such as '[76, 130]' \
	--number_of_final_designs num_of_binder_design # default 130 \
	--target_hotspot null
```
Alternatively, you can also generate the BindCraft input json file manually following the guideline in https://github.com/martinpacesa/BindCraft.

## Step5: Runing BindCraft binder design pipeline, by running:

```bash
conda activate BindCraft
python3.10 -u ../BindCraft/bindcraft.py \
	--settings ../BindCraft/settings_target/your_case_name.json \
	--filters ../BindCraft/settings_filters/effector_filters3.json \
	--advanced ../BindCraft/settings_advanced/default_4stage_multimer.json 
```
This process would take several hours to generate ~100 final binder designs passing all filters as recommended.

## Step6: Summarizing the results of BindCraft

After finishing the binder design process, you can summarize the results of BindCraft by using the following scripts. This script would generate a summarize the sequences of binder which pass the filtration of BindCraft,

```bash
conda activate BindCraft
python \
step22_bindercraft_summary.py \
	--input_folder ../BindCraft/your_case_name/ \
	--output_csv ../BindCraft/your_case_name/bindcraft_summary.csv
```

## Step7: Validating the candidate binders using HelixFold3
The output binders designed by BindCraft should be further validated by HelixFold3. For each candidate binder, three structural models were predicted: the binder-pathogen protein heterodimer, the binder-binder homodimer, and the binder monomer, to evaluate their binding affinity, oligomerization potential, and structural stability, respectively. Binders were filtered out if the inter-chain predicted TM-score (ipTM) of the heterodimer was below 0.85, the ipTM of the homodimer was above 0.6, or the predicted TM-score (pTM) of the monomer was below 0.85. We suggest to use the api of HelixFold3.

To predict the structures of binder monomer, run the script:

```bash
conda activate phsdk
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_monomer \
	--input_name_column 'Design' \
	--antigen_column 'antigen_sequence' \
	--nanobody_column 'binder_sequence' \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--num_of_structures_per_sequence 1 \
	--iptm_cut_off 0.71 \
	--only_filtering 0 \
	--summarize_results 0 \
	--max_job_num 150 \
	--helixfold_recycle 10 \
	--helixfold_ensemble 1 \
	--binder_or_nanobody 'binder_monomer' \
	--execute_quiet False \
	--require_dna_seq 0
```

To predict the structures of binder dimer , run the script:

```bash
conda activate phsdk
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_dimer \
	--input_name_column 'Design' \
	--antigen_column 'antigen_sequence' \
	--nanobody_column 'binder_sequence' \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--num_of_structures_per_sequence 1 \
	--iptm_cut_off 0.71 \
	--only_filtering 0 \
	--summarize_results 0 \
	--max_job_num 150 \
	--helixfold_recycle 10 \
	--helixfold_ensemble 1 \
	--binder_or_nanobody 'binder_dimer' \
	--execute_quiet False \
	--require_dna_seq 0
```

To predict the structures of pathogen protein-binder dimer, run the script:

```bash
conda activate phsdk
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_effector_dimer \
	--input_name_column 'Design' \
	--antigen_column 'antigen_sequence' \
	--nanobody_column 'binder_sequence' \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--num_of_structures_per_sequence 1 \
	--iptm_cut_off 0.71 \
	--only_filtering 0 \
	--summarize_results 0 \
	--max_job_num 150 \
	--helixfold_recycle 10 \
	--helixfold_ensemble 1 \
	--binder_or_nanobody 'binder' \
	--execute_quiet False \
	--require_dna_seq 0
```

After finishing the HelixFold3 and downloading all the result, run the scripts to summary the results of HelixFold3.

```bash
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_monomer \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--iptm_cut_off 0.71 \
	--only_filtering 1 \
	--summarize_results 1 \
	--binder_or_nanobody 'binder_monomer' \
	--require_dna_seq 0 
```
```bash
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_dimer \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--iptm_cut_off 0.71 \
	--only_filtering 1 \
	--summarize_results 1 \
	--binder_or_nanobody 'binder_dimer' \
	--require_dna_seq 0 
```
```bash
python3 step14_run_helixfold3_api_in_folder_hczhu.py \
	--input_csv_file --output_csv ../BindCraft/your_case_name/bindcraft_summary.csv\
	--input_type csv \
	--output_folder ../BindCraft/your_case_name/HF3_binder_effector_dimer \
	--proportion2next_step 0 \
	--num_of_sequence2next_step 0 \
	--iptm_cut_off 0.71 \
	--only_filtering 1 \
	--summarize_results 1 \
	--binder_or_nanobody 'binder' \
	--require_dna_seq 1 
```

If “--require_dna_seq” is set of 1, the script would generate the DNA sequences of binder. The ipTM and pTM scores of each state of each binder are recorded in the output file of HelixFold3. Binders with ideal predicted performance could be synthesized for further experimental validation.


## Credits
This repository uses code from:
AlphaFold 3 server (https://alphafoldserver.com/),
BindCraft (https://github.com/martinpacesa/BindCraft),
Helixfold3 (https://github.com/PaddlePaddle/PaddleHelix/tree/dev/apps/protein_folding/helixfold3, https://gitee.com/paddlehelix/paddlehelix)






