#!/usr/bin/env python3


import sys
import argparse
import csv
import json
import glob
from collections import defaultdict
import re
import os
import gc
from functools import cmp_to_key
from math import log2
from Bio import SeqIO
import pandas as pd

#copied/adapted FROM https://github.com/kipoi/models/blob/master/deepTarget/bio_utils.py
class bio_utils:
	# definitions for candiate lists
	SEED_START = 1
	SEED_SIZE = 6  # 6-mer seed
	CTS_SIZE = 40

	wc_pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}

	@staticmethod
	def getRCSeed(mirna_sequence):
		rev_sequence = mirna_sequence[::-1]
		seed = rev_sequence[-bio_utils.SEED_START - bio_utils.SEED_SIZE:-bio_utils.SEED_START]
		rc_seed = ""
		for i in seed:
			rc_seed += bio_utils.wc_pairs[i]
		return rc_seed

	@staticmethod
	def find_candidate_positions(mirna_sequence, mrna_sequence):
		positions = []
		rc_seed=bio_utils.getRCSeed(mirna_sequence)
		it = re.finditer(rc_seed, mrna_sequence)
		for match in it:
			positions.append(match.span()[1] + 1)
		return positions

def getRecs(fasta_file,verify_unique_names=True):
	rec_list=list()
	seq_list=list()
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		seq_list.append(str(seq_record.seq).upper())
		rec_list.append(str(seq_record.id))
	if verify_unique_names:
		rec_set=set(rec_list)
		if len(rec_set)==len(rec_list):
			#if the len of the data as a set is the same as the len of the data as a list then each item is unique!
			pass
		else:
			count_dict=defaultdict(int)
			for rec in rec_list:
				count_dict[rec]=count_dict[rec]+1
			dup_items=set([k for k in count_dict if count_dict[k]>1])
			raise Exception(f"In fasta file {fasta_file}, found the following {len(dup_items)} duplicate items : {dup_items}")
	return seq_list,rec_list

def setup_analysis(row,tmp_dir):
	mrna_fasta_t_to_u_name="mRNA.T_to_U.fa"
	mirna_name=row['mirna_id']
	analysis_id=row['summary_id']
	analysis_dict=dict()
	analysis_dict['mirna_fasta_file']=f"{tmp_dir}/analysis_{analysis_id}.fa"
	analysis_dict['mrna_fasta_file']=f"{tmp_dir}/{mrna_fasta_t_to_u_name}"
	analysis_dict['query_pair_file']=f"{tmp_dir}/query_pair_file_{analysis_id}"
	return analysis_dict
#==> example/mrna_fasta_file <==
#>NM_003629
#AGAGGAAGUGGGAAGA
#==> example/mrna_fasta_file.orig.txt <==
#>NM_003629
#AGAGGAAGUGGGAAG
#>NM_001135041
#==> example/query_pair_file <==
#hsa-miR-33b-5p	NM_005502
#hsa-miR-26a-5p	NM_005502

def writeAnalysisFiles(analysis_dict,mirna_name,mirna_seq,mrna_name,mrna_seq,write_mrna_file=False):

	#use the analysis dict to know the names of the files to write

	#write mirna query file
	with open(analysis_dict['mirna_fasta_file'],'w') as writer:
		writer.write(f">{mirna_name}\n{mirna_seq}")
	
	#write the query pair file
	with open(analysis_dict['query_pair_file'],'w') as writer:
		writer.write(f"{mirna_name}\t{mrna_name}")

	#write the mrna fasta if needed
	if write_mrna_file:
		with open(analysis_dict['mrna_fasta_file'],'w') as writer:
			writer.write(f">{mrna_name}\n{mrna_seq}")






if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='analyze a directory of circtest outputs',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("MRNA_FASTA",type=str,help="MRNA fasta ; should contain ONLY ONE record ; any T found are converted to U!")
	parser.add_argument("MIRNA_FASTA",type=str,help="path to mir fasta/database file")
	parser.add_argument("OUT_PLAN",type=str,help="path to save analysis plan/results summary table")
	parser.add_argument("TMP_DIR",type=str,default="tmp",help="path to temp dir")
	args = parser.parse_args()

	
	if(args):

		#verify TMP is a directory
		if not os.path.isdir(args.TMP_DIR):
			raise Exception(f"The temp directory {args.TMP_DIR} does not exist or is not a directory!")
		if not os.access(args.TMP_DIR, os.W_OK|os.R_OK):
			raise Exception(f"The temp directory {args.TMP_DIR} is neither readable nor writable!")
		   
		#first aquire the sequences and see what miRNAs have seeds in which mRNAs
		mirna_fasta,mirna_rec_list=getRecs(args.MIRNA_FASTA)
		mrna_fasta,mrna_rec_list=getRecs(args.MRNA_FASTA)
		if len(mrna_fasta)!=1:
			raise Exception(f"In file {args.MRNA_FASTA}, found {len(mrna_fasta)} sequences, but expected exactly one!")
		print(f"Loaded {len(mirna_fasta)} sequences from {args.MIRNA_FASTA}")
		print(f"Loaded {len(mrna_fasta)} sequence from {args.MRNA_FASTA}")
		if len(mirna_fasta)<=0:
			raise Exception(f"Error, too few mirna sequences found!")
		
		#verify no T in mRNA ; convert to U if any found
		if 'T' in mrna_fasta[0]:
			print(f"NOTE, converting T to U in the mRNA!")
			new_seq=[base if base!='T' else 'U' for base in mrna_fasta[0]]
			print(f"Old seq is\n{mrna_fasta[0]} and new seq is \n{new_seq}")
			mrna_fasta[0]="".join(new_seq)


		#from the mirnas, get the seeds
		mirna_seeds=[bio_utils.getRCSeed(mirna_fasta_seq) for mirna_fasta_seq in mirna_fasta]

		#for each mirna determine how many candidate sites are in the mRNA
		seed_positions=[bio_utils.find_candidate_positions(mirna_sequence,mrna_fasta[0]) for mirna_sequence in mirna_fasta]
		seed_positions_csvs=[",".join([str(x) for x in p]) for p in seed_positions]
		num_seed_positions=[len(positions) for positions in seed_positions]

		#for each mirna, if there are no candidate positions mark it to be skipped in analysis
		process_this_mirna=[p>0 for p in num_seed_positions]

		#make a short table indicating the analysis plan
		analysis_plan_df=pd.DataFrame.from_dict(
			{
				'summary_id':[x for x in range(len(mirna_rec_list))],
				'mirna_id':mirna_rec_list,
				'mirna_seq':mirna_fasta,
				'mirna_seed':mirna_seeds,
				'mirna_num_candidate_positions':num_seed_positions,
				'seed_positions':seed_positions_csvs,
				'analyze_mirna':process_this_mirna,
				'mrna_name':[mrna_rec_list[0] for x in range(len(mirna_rec_list))]
			}
		)
		print(analysis_plan_df.to_string(index=False))

		#filter out work to do
		analysis_plan_todo=analysis_plan_df[analysis_plan_df.apply(lambda row:row['analyze_mirna'],axis=1)]

		#set up analysis dicts to write input files
		analysis_dicts=analysis_plan_todo.apply(lambda row:setup_analysis(row,args.TMP_DIR),axis=1).tolist()

		#load the model and run each analysis!
		import kipoi
		model = kipoi.get_model('deepTarget')
		analysis_prediction_results=dict()
		for a in range(len(analysis_dicts)):

			#write the analysis files for deeptarget
			analysis_dict=analysis_dicts[a]
			main_id=int(analysis_dict['query_pair_file'].split("_")[-1])
			if a==0:
				write_mrna_file=True
			else:
				write_mrna_file=False
			writeAnalysisFiles(analysis_dict,mirna_rec_list[main_id],mirna_fasta[main_id],mrna_rec_list[0],mrna_fasta[0],write_mrna_file)

			num_pred_expected=num_seed_positions[main_id]
			preds = model.pipeline.predict(analysis_dict, batch_size=4)
			num_pred_received=len(preds)
			print(f"Raw predictions are {preds}")
			print(f"Number predictions expected is {num_pred_expected} ; number received is {num_pred_received}")
			for pred_pair in preds:
				print(f"A pred pair is {pred_pair}")
				p_0=pred_pair[0].astype(float)
				p_1=pred_pair[1].astype(float)
				if not main_id in analysis_prediction_results:
					analysis_prediction_results[main_id]=list()
				analysis_prediction_results[main_id].append([p_0,p_1])


		#build a column in the summary df to hold the analysis results and save all the results
		results_list_for_table=list()
		for analysis_id in range(len(mirna_rec_list)):
			if analysis_id in analysis_prediction_results:
				#get the data in JSON form
				json_results=json.dumps(analysis_prediction_results[analysis_id])
				results_list_for_table.append(json_results)
			else:
				results_list_for_table.append("[]")
		analysis_plan_df['predictions']=results_list_for_table
		print(f"Writing analysis plan/results summary table to {args.OUT_PLAN}")
		analysis_plan_df.to_csv(args.OUT_PLAN,sep="\t")
