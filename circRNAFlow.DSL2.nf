

nextflow.enable.dsl=2


process runFQC {

  input:
  path FQGZ

shell:
'''
#simply run FASTQC on each .FASTQ.GZ file
fastqc !{FQGZ}
'''
}


process runFlexBar {



input:
	path '*'
	path 'adapter_fasta'
output:
    path 'flexbar_output/*.fastq.gz'

shell:
'''

cat !{adapter_fasta} > adapters.fasta
head --verbose adapters.fasta
mkdir -v output
ls -alht 
set -x
echo ${PWD}
OPBN=`find *.gz|head -1|sed -r "s/_R1//g"|sed -r "s/\\.gz//g"`;
flexbar -r `find *_R1_*.gz` -p `find *_R2_*.gz` -t output/${OPBN}  -z GZ -m 30 -u 0 -q TAIL -qt 28 -a adapters.fasta -qf sanger

#rename output files to have same name as input files ;
# also rename output directory so output is a bit more easily findable
mv -vi output flexbar_output
cd flexbar_output
for RN in `seq 1 2`; do
	IN_READ=`find ../*_R${RN}_* | xargs basename`;
	echo "IN_READ IS ${IN_READ}" ; 
	OUT_READ=`find *fastq_${RN}.*`; 
	echo "OUT_READ IS ${OUT_READ}" ; 
	#rename output file name to same as input
	mv -v ${OUT_READ} ${IN_READ} ; 	
done ;

'''

}


process mapAgainstRRNA {



input:
	file '*' 
	file '*' 
output:
	file 'adapter_removed_rRNA_filtered/*'




shell:
'''
ls -alht 
#RRNA DB NAME
RRNA_DB=`find *.rev.1.bt2|sed -r "s/\\.rev.1.bt2//g"`;
echo "Got RRNA_DB ${RRNA_DB}" ; 
# align reads against the rRNA reference in paired end fashion
set -x
bowtie2 -x ${RRNA_DB} -1 `find *_R1_*.fastq.gz` -2 `find *_R2_*.fastq.gz` --no-unal --threads 8 --un-conc-gz adapter_removed_rRNA_filtered.fastq.gz 1>/dev/null 2>/dev/null

#move output to a new dir and rename to input (maintain short filenames)
mkdir -v adapter_removed_rRNA_filtered
mv -vi adapter_removed_rRNA_filtered.fastq.* ./adapter_removed_rRNA_filtered
cd adapter_removed_rRNA_filtered
for RN in `seq 1 2`; do
	IN_READ=`find ../*_R${RN}_* | xargs basename`;
	echo "IN_READ IS ${IN_READ}" ; 
	OUT_READ=`find *.${RN}.gz`; 
	echo "OUT_READ IS ${OUT_READ}" ; 
	#rename output file name to same as input
	mv -v ${OUT_READ} ${IN_READ} ; 	
done ;




'''


}


//include STAR modules
include { star_align_mode as star_align_pair } from './modules/star.align.nf'
include { star_align_mode as star_align_first } from './modules/star.align.nf'
include { star_align_mode as star_align_second } from './modules/star.align.nf'



process DCC_step {

input:
	file '*'
	file 'gtf_file'
	file 'fa_file'
	file 'repeat_file'
output:
	tuple file('LinearCount'), file('CircRNACount'), file('CircCoordinates')
	
shell:
'''
#memory report
free -g
#setup sample sheet and BAM list
if [ -f samplesheet.txt ] ; then rm -v samplesheet.txt ; fi ;
if [ -f bam_files ] ; then rm -v bam_files ; fi ;
for INPUTF in `find input.*` ; do
	SNAME=`cat ${INPUTF}`;
	echo ${SNAME}.pair.chimeric.out.junction >> samplesheet.txt
	echo ${SNAME}.pair.Aligned.sortedByCoord.out.bam >> bam_files
	samtools index ${SNAME}.pair.Aligned.sortedByCoord.out.bam
done ;
#setup repeat file and GTF  and reffa
REPEAT_FILE="!{repeat_file}"
GTF_AN_FILE="!{gtf_file}"
REFSEQ="!{fa_file}";
samtools faidx ${REFSEQ} 
#don't forget to include fai file here
ls -alht ${REFSEQ}
ls -alht ${REFSEQ}.fai
ls -alht ${REPEAT_FILE} 
#mate mappings
if [ -f mate1.txt ] ; then rm -v mate1.txt ; fi ;
if [ -f mate2.txt ] ; then rm -v mate2.txt ; fi ;
for INPUTF in `find input.*` ; do
	echo "Found input ${INPUTF}" ; 
	SNAME=`cat ${INPUTF}`;
	echo "${SNAME}.first.chimeric.out.junction"  >> mate1.txt
	echo "${SNAME}.second.chimeric.out.junction"  >> mate2.txt
done ;
countthreshold=2
replicatethreshold=3
#countthreshold=2
#replicatethreshold=1
#  -Nr countthreshold replicatethreshold
#  countthreshold replicatethreshold [default: 2,5]
set -x
wc -l mate1.txt mate2.txt bam_files  samplesheet.txt
export CPUCOUNT=`cat /proc/cpuinfo |grep -Pic '^processor'`
echo "Computed CPUCOUNT as ${CPUCOUNT}" ;
DCC   --keep-temp     -T ${CPUCOUNT}   -D -G  -mt1 @mate1.txt -mt2 @mate2.txt -R ${REPEAT_FILE} -an ${GTF_AN_FILE} -Pi   -F -M -Nr ${countthreshold} ${replicatethreshold} -fg -B @bam_files  -A ${REFSEQ} @samplesheet.txt
set +x
'''
}

//include circtest modules
include { circtest as circtest } from './modules/circtest_process.nf'
include { circtest as circtest_plain } from './modules/circtest_process.nf'
include { circtest_plotting as circtest_plotting } from './modules/circtest_process.nf'
include { circtest_plotting as circtest_plain_plotting } from './modules/circtest_process.nf'

process run_cluster_profiler {

	input:
		file '*'
		file 'kegg_db'
	output:
		file 'cp_outputs.zip'

shell:
'''

#unzip the zip file
ZIP_FILE=`find *.zip`; 
unzip ${ZIP_FILE} ; 

#get the merged TSV
MERGED_TSV=`find */*.merged.tsv`; 
echo "Found merged TSV ${MERGED_TSV}" ; 

#extract the TSV data ; sig_DE data and filter.
# Also, get uniq in the process ; use simple rules to address "conflicts" in
# LOG2FC or sig_p
EXTRACT_ROWS_PY_SCRIPT="""
import pandas as pd
main_df=pd.read_csv('${MERGED_TSV}',sep='\t')
desired_cols=['Chr','Start','End','Gene','Strand','sig_p','Log2FC']

####################################
#filter by sig_p
max_sig_p=0.1
filt_df=main_df[main_df.apply(lambda row:row['sig_p']<max_sig_p,axis=1)]

#get the groups to collapse into
group_cols=['Chr','Start','End','Gene','Strand']
filt_groups=filt_df.groupby(group_cols).size().reset_index().rename(columns={0:'count'})

#perform the collapsing into new_df
import json
import numpy as np
new_df=None
new_df_all=None
def match_row(given_row,group_row):
    return all([group_row[gc]==given_row[gc] for gc in group_cols])
    
for index, group_row in filt_groups.iterrows():
    
    #first for the group, get all mathcing filtered
    #filter_for_group=filt_df[filt_df.apply(lambda)
    filt_rows_this_group=filt_df[filt_df.apply(lambda row:match_row(row,group_row),axis=1)]

    
    #to collapse it compute a new p-value ; use the mean
    new_sig_p_val=filt_rows_this_group['sig_p'].mean()
    
    #also compute a collapsed/new Log2FC ; convert from 
    # log2 to real, then average, then go back to log 2
    log2FCs=[float(x) for x in filt_rows_this_group['Log2FC'].tolist()]
    FCs=[2**l for l in log2FCs]
    avgFC=sum(FCs)/len(FCs)
    collapsedL2FC=np.log2(avgFC)

    #compute collaped Log2Ratio
    log2Ratios=[float(x) for x in filt_rows_this_group['Log2Ratio'].tolist()]
    ratios=[2**r for r in log2Ratios]
    avgRatio=sum(ratios)/float(len(ratios))
    collapsedLog2Ratios=np.log2(avgRatio)
    
    #create a row using the group and the collapsed values
    new_row_dict={
        grp_col_name:[group_row[grp_col_name]]
        for
        grp_col_name in group_cols
    }
    new_row_dict['sig_p']=[new_sig_p_val]
    new_row_dict['Log2Ratio']=[collapsedLog2Ratios]
    new_row_df=pd.DataFrame.from_dict(new_row_dict)
    
    #accumulate to new_df_all as well, not just the cols for cp_input, but all the cols
    temp_dict=dict()
    for col in filt_rows_this_group.columns:
        temp_dict[col]=filt_rows_this_group[col].tolist()
    temp_dict['sig_p']=[new_sig_p_val for x in range(len(log2Ratios))]
    temp_dict['Log2Ratio']=[collapsedLog2Ratios for x in range(len(log2Ratios))]
    temp_df=pd.DataFrame.from_dict(temp_dict)

    #merge the new row with the growing list of such rows
    if new_df is None:
        new_df=new_row_df.copy(deep=True)
        new_df_all=temp_df.copy(deep=True)
    else:
        new_df=pd.concat([new_df.copy(deep=True),new_row_df.copy(deep=True)])
        new_df_all=pd.concat([new_df_all.copy(deep=True),temp_df.copy(deep=True)])

#write to  cp_input.tsv (cluster_profiler input)
new_df.to_csv('cp_input.tsv', sep='\t',index=False)

#add the other output file
new_df_all.to_csv('full_annotated_cp_input.tsv',sep='\t',index=False)

"""
python3 -c "${EXTRACT_ROWS_PY_SCRIPT}" ; 

#run cluster profiler
KEGG_ORGANISM="!{params.kegg_cp_organism_str}"  KEGG_CP_FILE=kegg_db  CPSSEP=$'\\t'  Rscript /usr/local/bin/cp_script.R  cp_input.tsv Log2Ratio 0 Gene sig_p

#package the outputs
zip cp_outputs.zip  cp_input.tsv full_annotated_cp_input.tsv  cp_input.tsv.sig.*

'''

	}




workflow {

	//run FASTQC solo
	runFQC(Channel.fromPath(params.fqgzglob))

	//main pipeline:

	//1) after running FQC, pair up the FQs an then run them through flexbar for adapter removal.
	paired_fqgz=Channel.fromPath(params.fqgzglob).map { [  new File(""+it).getName().replace('.fastq.gz', '').replace('_R1_','_RX_').replace('_R2_','_RX_') , it ] }
		.groupTuple(by: 0)
		.map { it[1] }
	adapter_fasta=Channel.value(file(params.adapter_fasta))
	flexed=runFlexBar(paired_fqgz,adapter_fasta)

	//after running them through flexbar, use 
	// 2) BWA to map against an rRNA database for rRNA filtering
	rrRNACleaned=mapAgainstRRNA(flexed,
		channel.fromPath(params.rrna_glob).collect())

	//3) map the filtered data using STAR in 3 modes:
		star_ref=channel.fromPath(params.stardb_glob).collect()
		star_gtf=channel.fromPath(params.gtf_file).collect()
	// a) as a pair
	star_pair_tuple=star_align_pair(rrRNACleaned,star_ref,star_gtf,"pair")
	// b) just read 1
	star_first_tuple=star_align_first(rrRNACleaned,star_ref,star_gtf,"first")
	// c) just read 2
    star_second_tuple=star_align_second(rrRNACleaned,star_ref,star_gtf,"second")

	//4) combine the star results and merge with FQ files and input to DCC
	all_star_results=star_pair_tuple.mix(star_first_tuple)
		.mix(star_second_tuple)
		//this map associates a source R1FQ with its output
		.map { [new File(""+it[0]).text , it[1],it[2],it[3] ] }
	DCC_input=rrRNACleaned.map { [ (new File(""+it[0]).getName())+"" , it[0],it[1]] }
		.combine(all_star_results)
		//as STAR results are combined here, this filter makes sure that samples are kept with eachother (by string comparison on sample name)
		.filter { it[0].startsWith(it[3]) }
		.groupTuple()
		// -> pair name, R1 , R2, chimeric junction (3 per sample) , SJ.out(3 per sample) , Aligned.sortedByCoord.out.bam (3 per sample)
		.map{ [ it[0],it[1][0],it[2][0] , it[4][0],it[4][1],it[4][2] , it[5][0],it[5][1],it[5][2] , it[6][0],it[6][1],it[6][2]    ] }
		.collect()
	circtest_staging=DCC_step(
		DCC_input,
		star_gtf,
		channel.fromPath(params.fa_ref_file),
		channel.fromPath(params.repeat_file)
	)

	//5) prepare for circtest by combining DCC output with comparison configuration data
	comp_channel=Channel
    	.fromPath(params.comp_list)
    	.splitText()
	comp_circtest_channel_for_split=comp_channel.combine(circtest_staging)
    raw_circtest_results_by_cohort=circtest(comp_circtest_channel_for_split,channel.fromPath(params.cohort_comp_conf),"")
    raw_circtest_results_by_cohort_plain=circtest_plain(comp_circtest_channel_for_split,channel.fromPath(params.cohort_comp_conf),"plain")

    //6) run plotting
    plotting_results=circtest_plotting(raw_circtest_results_by_cohort,channel.fromPath(params.circatlas_bed),channel.fromPath(params.cohort_comp_conf))
    plotting_results_plain=circtest_plain_plotting(raw_circtest_results_by_cohort_plain,channel.fromPath(params.circatlas_bed),channel.fromPath(params.cohort_comp_conf))

    //7) run cluster profiler
    run_cluster_profiler(plotting_results,channel.fromPath(params.kegg_db))


}
