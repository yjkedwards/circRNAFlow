#!/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
import csv
import json
import glob
from collections import defaultdict
import re
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool
import gc
from functools import cmp_to_key
from math import log2

def getSettingsFromFileName(fn):
	fn_no_tsv=fn.replace(".tsv","")
	#filtercount_10_filtersample_9_percentage_0.05_alpha_1.tsv
	settings=['filtercount','filtersample','percentage','alpha']
	setting_dict=dict()
	for s in settings:
		temp_str=s+"_([0-9\.]+)"
		setting_re=re.compile(temp_str)
		sr=re.search(setting_re,fn_no_tsv)
		if(sr):
			tg=float(sr.group(1))
			#print(f"for setting {s} for file name {fn}, the value is {tg}")
			setting_dict[s]=tg
		else:
			raise Exception(f"For filename {fn}, could not find setting(key-value) for {s} !")
	#print("For file "+fn+", got settings : "+json.dumps(setting_dict,indent=4))
	return setting_dict



def isLocusInCAPanda(cadf,locus_dict):
	#tm56d=dict({'Chr':"3","Start":121235043,"End":121235764,"Strand":"-"})
	#locus_dict=tm56d
	matching_rows=cadf[(cadf['Chro']=="chr"+locus_dict['Chr']) & (cadf['Start']==locus_dict['Start']) & (cadf['End']==locus_dict['End']) & (cadf['Stand']==locus_dict['Strand'])]
	#matching_rows=cadf[cadf.apply(lambda row:filterFunc(row),axis=1)]
	#print(f"For {locus_dict},(panda)  matching rows are {matching_rows.to_string()}")
	if(matching_rows.shape[0]>0):
		return matching_rows['Name'].tolist()[0]
	else:
		return None
		

def isLocusInCA(cadf,locus_dict):
	#bed stuff
	#Chro	Start	End	Stand	Name
	#chr2	60229684	60232367	+	mmu-March7_0001
	#chr3	158014010	158016381	-	mmu-Srsf11_0024
	#chr4	141090422	141098047	+	mmu-Spata21_0002
	#locus cols
	#"Chr"	"Start"	"End"	"Gene"	"JunctionType"	"Strand"
	def filterFunc(row):
		if(str(row['Chro'])=="chr"+str(locus_dict['Chr']) and
			int(row['Start'])==int(locus_dict['Start']) and
			int(row['End'])==int(locus_dict['End']) and
			row['Stand']==locus_dict['Strand']):
			return True
		else:
			return False
	#print(cadf.columns)
	matching_rows=cadf[cadf.apply(lambda row:filterFunc(row),axis=1)]
	print(f"For {locus_dict}, matching rows are {matching_rows.to_string()}")



def getFilesInDir(d):
	fglob=d+"/*.tsv"
	tsvs=glob.glob(fglob)
	return tsvs
	

def getAllDataDFFromFile(f,cohort_name):
	main_dict=readCSVAsDict(f)
	settings_dict=getSettingsFromFileName(f)
	if(not('Start' in main_dict)):
		raise Exception(f"Could not find key 'Start' in file {f}, is it empty? or mal-formatted?")
	num_rows=len(main_dict['Start'])
	for s in settings_dict:
		for r in range(num_rows):
			main_dict[s].append(settings_dict[s])
			if(s=='filtercount'):
				main_dict['cohort'].append(cohort_name)
	return main_dict
	


def  readCSVAsDict(f):
	mydict=None
	keys=None
	with open(f, mode='r') as infile:
		reader = csv.reader(infile,delimiter='\t')
		for row in reader:
			if(keys==None):
				keys=row
			else:
				if(mydict==None):
					mydict=defaultdict(list)
				for k in range(len(keys)):
					mydict[keys[k]].append(row[k])
			#print(json.dumps(row,indent=4))
	#print(json.dumps(mydict,indent=4))
	return mydict


def makeDictUnion(dict_list):
	merged_dict=defaultdict(list)
	for d in dict_list:
		for k in d:
			vals=d[k]
			#print(f"for a dict, a key is {k}, and its vals {vals}")
			for v in vals:
				merged_dict[k].append(v)
	return merged_dict

def getCohortName(d):
	bn=os.path.basename(d)
	if(len(bn)==0):
		#print("empty")
		bnp=os.path.split(d)[-2]
		bn=getCohortName(bnp)
		return bn
	return bn

def getQueryDicts(df):
	def makeDict(row):	
		#dict({'Chr':"3","Start":121235043,"End":121235764,"Strand":"-"})
		#matching_rows=cadf[(cadf['Chro']=="chr"+locus_dict['Chr']) & (cadf['Start']==locus_dict['Start']) & (cadf['End']==locus_dict['End']) & (cadf['Stand']==locus_dict['Strand'])]
		d=dict()
		d['Chr']=row['Chr']
		d['Start']=row['Start']
		d['End']=row['End']
		d['Strand']=row['Strand']
		return d
	dict_list=df.apply(lambda row:makeDict(row),axis=1).tolist()
	#print(json.dumps(dict_list,indent=4))
	return dict_list


def getUnique(merged_df,circAtlasDF):
	df_i=merged_df.copy()[['Chr','Start','End','Strand']]
	unique_rows=df_i.drop_duplicates()
	def getQD(row):
		keys=['Chr','Start','End','Strand']
		d=dict()
		for key in keys:
			d[key]=row[key]
			if(key in ['Start','End']):
				d[key]=int(d[key])
		return d
	dict_list=unique_rows.apply(lambda row:getQD(row),axis=1).tolist()
	#print(json.dumps(dict_list[0],indent=4))
	#print(json.dumps(dict_list,indent=4))
	tflist=[isLocusInCAPanda(circAtlasDF,locus_dict) for locus_dict in dict_list]
	main_lookup=dict()
	for d_i in range(len(dict_list)):
		csv_key=",".join([str(dict_list[d_i][k]) for k in dict_list[d_i]])
		csv_val=tflist[d_i]
		main_lookup[csv_key]=csv_val
	#print(json.dumps(main_lookup,indent=4))
	return main_lookup
		

def addCAFlagToDF(merged_df,CA_lookup_dict):
	def getFlag(row):
		keys=['Chr','Start','End','Strand']
		luk=",".join([str(row[k]) for k in keys])
		return CA_lookup_dict[luk]
	merged_df['CircAtlas']=merged_df.apply(lambda row:getFlag(row),axis=1)
	return merged_df

def makeListOfLociRows(merged_df):
	cols_to_uniq=['Chr',"Start","End","Strand"]
	merged_copy=merged_df.copy(deep=True)
	the_uniques=merged_copy[cols_to_uniq].drop_duplicates()
	def makeDict(row):
		d=dict()
		for c in cols_to_uniq:
			d[c]=row[c]
		d['Start']=int(d['Start'])
		d['End']=int(d['End'])
		return d
	list_of_dicts=the_uniques.apply(lambda row:makeDict(row),axis=1).tolist()
	#print(json.dumps(list_of_dicts,indent=4))
	return list_of_dicts
	

		
def makeSinglePlot(lookupkey,df,group_cols,out_plot_dir,adj_pv_thresh):
	#lookupkey=dict({'Chr':"3","Start":121235043,"End":121235764,"Strand":"-"})
	#lookupkey=dict({'Chr':"5","Start":"14045001","End":"14045372","Strand":"+"})
	#print(json.dumps(lookupkey,indent=4))
	kc=lookupkey['Chr']
	ks=str(lookupkey['Start'])
	ke=str(lookupkey['End'])
	kss=lookupkey['Strand']
	mrows=df[(df['Chr']==kc) & (df['Start']==ks) & (df['End']==ke) & (df['Strand']==kss)]
	ca_set=set(mrows['CircAtlas'].tolist())
	isInCircAtlas=any([str(ca).startswith("mmu-") for ca in ca_set])
	mrows_sorted=mrows.sort_values('sig_p', axis=0, ascending=True, inplace=False)
	mrows_sorted.reset_index(inplace=True)
	num_rows_under_thresh=sum(mrows_sorted.apply(lambda row:float(row['sig_p'])<=adj_pv_thresh,axis=1).tolist())
	#print(f"For key \n{json.dumps(lookupkey,indent=4)},\n num rows with sig_p<={adj_pv_thresh} is {num_rows_under_thresh} and its circatlas status is {isInCircAtlas}")
	if(not(isInCircAtlas) and num_rows_under_thresh==0):
		#if insignificant and not in circatlas, then skip
		pass
	else:
		if(isInCircAtlas):
			#if in circatlas, then plot at least one
			if(mrows_sorted.shape[0]>=1):
				head_val=max(num_rows_under_thresh,1)
				plotRows(mrows_sorted.head(head_val),group_cols,out_plot_dir)
		else:
			#if not in circatlas, plot everything sufficiently significant
			if(mrows_sorted.shape[0]>=1 and num_rows_under_thresh>=1):
				plotRows(mrows_sorted.head(num_rows_under_thresh),group_cols,out_plot_dir)

def getIndicesOfDataCols(df):
	all_cols=list(df.columns)
	ratio2idx=all_cols.index("group_2_ratio_mean")
	pvalidx=all_cols.index("p.val")
	data_idxs=[x for x in range(len(all_cols)) if ratio2idx<x and x<pvalidx]
	#raise Exception(f"For all cols \n {all_cols}, the data_idxs are {data_idxs} and they are {[all_cols[d] for d in data_idxs]}")
	return data_idxs

	


def plotRows(df,group_cols,out_plot_dir):
	#assumes already filtered by p-val(adj) etc.
	def getTitle(row):
		t=row['Chr']+":"+row['Start']+"-"+row['End']+"/"+row['Strand']+"/"+row['Gene']+" "
		#t=t+"P(adj)="+str(row['sig_p'])
		#t=t+"P(adj)="+str(float(np.round(float(row['sig_p']),4)))+" "
		#t=t+"("+str(float(np.round(float(row['group_1_ratio_mean']),4)))+" vs "+str(float(np.round(float(row['group_2_ratio_mean']),4)))+")"
		#filtercount	cohort	filtersample	percentage
		#row_names=['filtercount','filtersample','percentage']
		#dnames=['FC','FS','PCT']
		#farr=[dnames[i]+"="+str(row[row_names[i]]) for i in range(len(row_names))]
		#t=t+" ".join(farr)
		if('CircAtlas' in row):
			if(row['CircAtlas'] is not None):
				if(len(row['CircAtlas'])>1):
					t=t+" "+row['CircAtlas']

		return t

	#group by profile
	#To group by profile, first get the col idxs of the data cols
	data_idxs=getIndicesOfDataCols(df)
	all_names=list(df.columns)
	data_names=[all_names[d] for d in data_idxs]
	#Second, for each row, get the profile using those cols
	def getProfile(row):
		the_profile=[row[n] for n in data_names]
		return the_profile
	row_profiles=df.apply(lambda row:getProfile(row),axis=1).tolist()
	#Third for each row, assign its profile an ID by clustering (exact match)
	clusters=list()
	cluster_idxs=list()
	for profile in row_profiles:
		if(not(profile in clusters)):
			#if a never-seen profile add it to the list and create new id
			clusters.append(profile)
			cluster_idxs.append(len(clusters)-1)
		else:
			cluster_idxs.append(clusters.index(profile))
	#now clusters is a list of clusters
	#and cluster_idxs is a list of the ids for each row
	#Finally to each row add the cluster id
	df_copy=df.copy(deep=True)
	df_copy['ProfileClusterID']=cluster_idxs
	df_copy.reset_index(inplace=True,drop=True)
	unique_cluster_ids=set(cluster_idxs)

	

	pdf_name_keys=['Chr','Start','End','Gene']
	pdf_name=out_plot_dir+"/"+".".join([str(df[k].tolist()[0]) for k in pdf_name_keys])
	strand_str="plus"
	if(df.head(1)['Strand'].tolist()[0]=="-"):
		strand_str="minus"
	pdf_name=pdf_name+"."+strand_str+".pdf"
	if(len(unique_cluster_ids)>1):
		print(f"Note : multiple plots to be generated in {pdf_name}")
	with PdfPages(pdf_name) as pdf:
		# As many times as you like, create a figure fig and save it:
		for cluster_id in unique_cluster_ids:
			rows_this_is=df_copy[df_copy['ProfileClusterID']==cluster_id]
			for index, row in rows_this_is.head(1).iterrows():
				#print(row['colA'], row['colB'], row['colC'])
				all_cols=[]
				colors=[]
				for b_i in range(len(group_cols)):
					grp=group_cols[b_i]
					for g in grp:
						all_cols.append(g)
						if(b_i==0):
							colors.append('red')
						else:
							colors.append('blue')
				plt.figure().clear()
				fig = plt.figure()
				x=all_cols
				#print(f"all_cols is {all_cols}")
				y=[int(v) for v in row[all_cols].values.tolist()]
				plt.ylim([0,max(y)+5])
				#print(f"y is {y}")
				x_pos = [i for i, _ in enumerate(x)]
				plt.bar(x_pos, y, color=colors,zorder=10)
				plt.xlabel('Group')
				plt.ylabel('Count')
				plt.grid(zorder=0)
				plt.title(getTitle(row))
				for i, v in enumerate(y):
					num_digs=len(str(v))
					dig_off_val=0.1
					x_offset=(-1.0)*(dig_off_val*float(num_digs))-0.05
					plt.text(i+x_offset ,v+0.5, str(v), color='black',fontweight='bold',zorder=10)
					
				#plt.yticks([i for i in range(max(y)+1)])
				plt.xticks(x_pos,x,rotation=90)
				pdf.savefig(fig,bbox_inches='tight')
			
#		fig = plt.figure()
#		x = ['Nuclear', 'Hydro', 'Gas', 'Oil', 'Coal', 'Biofuel']
#		energy = [5, 6, 15, 22, 24, 8]
#		x_pos = [i for i, _ in enumerate(x)]
#		plt.bar(x_pos, energy, color='green')
#		plt.xlabel("Energy Source")
#		plt.ylabel("Energy Output (GJ)")
#		plt.title("Energy output from various fuel sources")
#		plt.xticks(x_pos, x)
#		pdf.savefig(fig)


def getGroupCols(df,count_in_first_group):
	the_cols=list(df.columns)
	#sig_p	group_1_ratio_mean	group_2_ratio_mean	WT_veh_ML_6381	WT_veh_FL_6209	WT_veh_ML_6425	WT_veh_FL_6456	WT_veh_ML_6205	WT_veh_FL_6529	CF_veh_ML_7030	CF_veh_FL_6294	CF_veh_ML_6455	CF_veh_FL_7483	CF_veh_ML_6528	CF_veh_FL_7576	p.val
	#the cols of the groups are between sig_p and p.val
	sig_p_idx=the_cols.index("group_2_ratio_mean")
	p_val_idx=the_cols.index("p.val")
	names_cols=[the_cols[i] for i in range(len(the_cols)) if sig_p_idx<i and i<p_val_idx]
	g1_cols=names_cols[:count_in_first_group]
	g2_cols=names_cols[count_in_first_group:]
	return [g1_cols,g2_cols]


def addLog2FCToDF(df,count_in_first_group,comp_obj,cohort_name):
	group_cols=getGroupCols(df,count_in_first_group)
	#raise Exception(f"eddie, group_cols are {group_cols}")
	df_copy=df.copy(deep=True)
	def safeFC(row):
		if('NA' in [row['group_1_ratio_mean'],row['group_2_ratio_mean']]):
			return float('inf')
		if(float(row['group_1_ratio_mean'])>0):
			return float(row['group_2_ratio_mean'])/float(row['group_1_ratio_mean'])
		else:
			return float('inf')
	computed_fc=df_copy.apply(lambda row:safeFC(row),axis=1).tolist()
	def safeLog2(fc):
		if(fc==float('inf')):
			return float('inf')
		if(fc!=0):
			return log2(abs(fc))
		else:
			return float('inf')
	log2fc=[safeLog2(fc) for fc in computed_fc]
	df_copy['FoldChange']=computed_fc
	df_copy['Log2FC']=log2fc
	first_group_cols=group_cols[0]
	second_group_cols=group_cols[1]	
	first_group_mean=df_copy.apply(lambda row:float(sum([float(row[fgc]) for fgc in first_group_cols]))/float(len(first_group_cols)),axis=1)
	secnd_group_mean=df_copy.apply(lambda row:float(sum([float(row[sgc]) for sgc in second_group_cols]))/float(len(second_group_cols)),axis=1)
	grp_ratios=[secnd_group_mean[i]/first_group_mean[i] for i in range(len(first_group_mean))]
	log2_grp_ratios=[safeLog2(gr) for gr in grp_ratios]
	first_grp_name=comp_obj['comp_names'][cohort_name][0]
	secnd_grp_name=comp_obj['comp_names'][cohort_name][1]
	df_copy['Mean_'+first_grp_name]=first_group_mean
	df_copy['Mean_'+secnd_grp_name]=secnd_group_mean
	df_copy['Ratio']=grp_ratios
	df_copy['Log2Ratio']=log2_grp_ratios
	return df_copy
	
	


def paraPlotting(t):
	the_dict=t[0]
	df=t[1]
	group_cols=t[2]
	out_pdf=t[3]
	pvc=t[4]
	makeSinglePlot(the_dict,df,group_cols,out_pdf,pvc)
	gc.collect()


def compareLociDict(d1,d2):
	#first sort on chrom
	map_chrom={	"1": 1,
	"2": 2,
	"3": 3,
	"4": 4,
	"5": 5,
	"6": 6,
	"7": 7,
	"8": 8,
	"9": 9,
	"10": 10,
	"11": 11,
	"12": 12,
	"13": 13,
	"14": 14,
	"15": 15,
	"16": 16,
	"17": 17,
	"18": 18,
	"19": 19,
	"20": 20,
	"21": 21,
	"22": 22,
	"23": 23,
	"X": 24,
	"Y": 25,
	"MT":26}



	if(d1["Chr"] in map_chrom and d2["Chr"] in map_chrom):
		#first compare on chrom
		d1_mapped=map_chrom[d1["Chr"]]
		d2_mapped=map_chrom[d2["Chr"]]
		if(d1_mapped<d2_mapped):
			return -1
		elif(d2_mapped<d1_mapped):
			return 1
	else:
		#for any chrom not on here return 1 to signify always last
		return 1
	#now compare on starts
	if(d1["Start"]<d2["Start"]):
		return -1
	else:
		return 1	



if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='analyze a directory of circtest outputs',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("CT_DIR",type=str,help="directory of a cohort containing outputs of circtesting")
	parser.add_argument("CA_BED",type=str,help="path to CircAtlas BED file")
	parser.add_argument("OUT_TSV",type=str,help="output summary TSV file")
	parser.add_argument("OUT_PDF_DIR",type=str,help="directory in which to write PDF files of plots")
	parser.add_argument("COHORT_COMP_JSON",type=str,help="configuration JSON for comparisons")
	parser.add_argument("-PVC",type=float,default=0.1,help="adjusted p-value threshold")
	args = parser.parse_args()
	print("Running plotter!")
	if(args):
		CT_DIR=args.CT_DIR
		CA_BED=args.CA_BED
		circAtlasDF=pd.read_csv(CA_BED, sep='\t' )
		cohort_name=getCohortName(CT_DIR)
		print(f"Loading with cohort name {cohort_name}")
		tsvs=getFilesInDir(CT_DIR)
		dict_list=list()
		count_in_first_group=None
		reader=open(args.COHORT_COMP_JSON)
		comp_obj=json.load(reader)
		count_in_first_group=comp_obj['comparison_counts'][cohort_name]
		print(f"Obtained count {count_in_first_group} for number in first group comparison !")
		print(f"Loading {len(tsvs)} TSVs...")
		if(len(tsvs)==1):
			if(os.stat(tsvs[0]).st_size==0):
				print("Length of a single input TSV is zero!  No plotting to do!")
				sys.exit(0)
		for tsv in tsvs:
			temp_dict=getAllDataDFFromFile(tsv,cohort_name)
			dict_list.append(temp_dict)
		print("Processing for data frame")
		merged_dict=makeDictUnion(dict_list)
		merged_df=pd.DataFrame.from_dict(merged_dict)
		print(f"Indexing CircAtlas for lookup...")
		CA_lookup_dict=getUnique(merged_df,circAtlasDF)
		print("Incorporating CircAtlas data...")
		ann_merge_df=addCAFlagToDF(merged_df,CA_lookup_dict)
		ann_merge_df_wfc=addLog2FCToDF(ann_merge_df,count_in_first_group,comp_obj,cohort_name)
		ann_merge_df_wfc.to_csv(args.OUT_TSV,index=False,sep="\t")
		print(f"Wrote to summary file {args.OUT_TSV}")
		group_cols=getGroupCols(ann_merge_df,count_in_first_group)
		dicts_for_plotting=makeListOfLociRows(ann_merge_df)
		dicts_for_plotting=sorted(dicts_for_plotting,key = cmp_to_key(compareLociDict))
		print(f"Found {len(dicts_for_plotting)} unique loci to consider plotting for...")
		#raise Exception(f"{json.dumps(dicts_for_plotting[:10],indent=4)}")
		print(f"To plot {len(dicts_for_plotting)} difference putative circRNA loci.")
		print(f"Using (adjusted) p-value threshold : {args.PVC}")
		chroms_plotted=set()
		for x in range(len(dicts_for_plotting)):
			#print(f"Running plotting for {dicts_for_plotting[x]}")
			this_chrom=dicts_for_plotting[x]['Chr']
			if(not(this_chrom in chroms_plotted)):
				print(f"Scanning/plotting chromosome {this_chrom}...")
				chroms_plotted.add(this_chrom)
			makeSinglePlot(dicts_for_plotting[x],merged_df,group_cols,args.OUT_PDF_DIR,args.PVC)
	else:
		parser.print_args()
