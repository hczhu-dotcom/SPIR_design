# SPIR_design
Code for SPIR design

This github is used for designing of Synthetic Plant Immune Receptor (SPIR)

Step1: predecting the structure of plant pathogen protein using AlphaFold 3,by using its webserver https://alphafoldserver.com/.
Step2: triming the N and C unstructured region of the plant pathogen protein, resulting pathogen_protein.pdb, such as PVY_CP.pdb



```bash
conda activate hczhu_177
python /public-supool/home/gaolab/haochengzhu/My_script_hczhu_177/step22_bindercraft_summary.py \
	--input_folder /public-supool/home/gaolab/haochengzhu/BindCraft/my_cases/$case_name \
	--output_csv /public-supool/home/gaolab/haochengzhu/BindCraft/my_cases/$case_name/bindcraft_summary.csv
```

