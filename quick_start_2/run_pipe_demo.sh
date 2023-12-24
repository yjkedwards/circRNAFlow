#!/bin/bash

#define variables for timestamped log files
DT=`/bin/date --rfc-3339=ns|tr ":" "_"|tr " " "T"`;
NF_OUT="${DT}.nf.out" ;
NF_ERR="${DT}.nf.err" ; 

# run pipeline with specific version and indicated sample/test data
PIPE_FILE="circRNAFlow.DSL2.nf" ; 
set -x
nextflow -C  demo.sing.config   run ${PIPE_FILE} --cohort_comp_conf pipe_data/cohort_comp_conf.json \
	--comp_list pipe_data/comp_list.txt  --adapter_fasta sample_data/adapters.fasta \
	--fqgzglob 'sample_data/*.gz'   \
	--rrna_glob 'ref_data/hg38_rrna_prep/human_rrnas.f*'  \
	--gtf_file ref_data/Homo_sapiens.GRCh38.107.gtf  \
	--stardb_glob 'ref_data/star_db/*'   \
	--fa_ref_file ref_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa  \
	--repeat_file ref_data/human_repeats.gtf   \
	--circatlas_bed ref_data/circatlas_human_bed_v2.0.txt \
	--kegg_cp_organism_str hsa  \
	--kegg_db ref_data/symbol_to_kegg_gid.HUMAN.csv \
	--craft_input_glob_str "/mnt/sdc/five_tera/UAB_YE_CF/CRAFT_proj_dir_MRO_READONLY_NOCIRCRNA_DATA/input/*" \
	--craft_params "/mnt/sdc/five_tera/UAB_YE_CF/CRAFT_proj_dir_MRO_READONLY_NOCIRCRNA_DATA/params.txt" \
	--craft_ref_path_file "Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
	--deeptarget_mirna_fa "/mnt/sdc/five_tera/UAB_YE_CF/circRNAFlow/singularity/KipoiDeepTarget/example/mirna_fasta_file"  \
	-with-report report.html  \
	-with-trace trace.txt -with-dag flowchart.png -resume #1>${NF_OUT} 2>${NF_ERR}

