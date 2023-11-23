

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



process star_align_pair {



input:
	path '*' 
	path 'star_db/*'
	path 'gtf_file'
output:
	tuple file('source_pair.txt'),file('*.pair.chimeric.out.junction'),file('*.pair.SJ.out.tab'),file('*.pair.Aligned.sortedByCoord.out.bam') 

shell:
'''
ls -alht
find .
READ_1=`find *_R1_*gz`
READ_2=`find *_R2_*gz`
STAR --runThreadN 8 \
	--genomeDir ${PWD}/star_db  \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn ${READ_1} ${READ_2} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${PWD}/paired_output/ \
	--outReadsUnmapped Fastx \
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS \
	--outSJfilterOverhangMin 15 15 15 15 \
	--alignSJoverhangMin 15 \
	--sjdbGTFfile !{gtf_file} \
	--alignSJDBoverhangMin 15 \
	--outFilterMultimapNmax 20 \
	--outFilterScoreMin 1 \
	--outFilterMatchNmin 1 \
	--outFilterMismatchNmax 2 \
	--chimSegmentMin 15 \
	--chimScoreMin 15 \
	--chimScoreSeparation 10 \
	--chimJunctionOverhangMin 15 1>star.out 2>star.err
mv -vi paired_output/Chimeric.out.junction paired_output/pair.chimeric.out.junction && touch paired_output/pair.chimeric.out.junction
mv -vi paired_output/SJ.out.tab paired_output/pair.SJ.out.tab && touch  paired_output/pair.SJ.out.tab
mv -vi paired_output/Aligned.sortedByCoord.out.bam paired_output/pair.Aligned.sortedByCoord.out.bam  && touch paired_output/pair.Aligned.sortedByCoord.out.bam
mv -vi paired_output/pair.* .
/bin/echo -ne "${READ_1}" > source_pair.txt
for F in `find pair.*`; do
	mv -v ${F} ${READ_1}.${F} ;
done ;
chmod -vR 777 paired_output
sleep 5
'''

}




process star_align_first {


input:
	path '*' 
	path 'star_db/*'
	path 'gtf_file'
output:
	tuple file('source_pair.txt'),file('*.first.chimeric.out.junction'),file('*.first.SJ.out.tab'),file('*.first.Aligned.sortedByCoord.out.bam') 

shell:
'''
ls -alht
find .
READ_1=`find *_R1_*gz`
READ_2=`find *_R2_*gz`
STAR --runThreadN 8 \
	--genomeDir ${PWD}/star_db  \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn ${READ_1} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${PWD}/first_output/ \
	--outReadsUnmapped Fastx \
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS \
	--outSJfilterOverhangMin 15 15 15 15 \
	--alignSJoverhangMin 15 \
	--sjdbGTFfile !{gtf_file} \
	--alignSJDBoverhangMin 15 \
	--outFilterMultimapNmax 20 \
	--outFilterScoreMin 1 \
	--outFilterMatchNmin 1 \
	--outFilterMismatchNmax 2 \
	--chimSegmentMin 15 \
	--chimScoreMin 15 \
	--chimScoreSeparation 10 \
	--chimJunctionOverhangMin 15 1>star.out 2>star.err
mv -vi first_output/Chimeric.out.junction first_output/first.chimeric.out.junction && touch first_output/first.chimeric.out.junction
mv -vi first_output/SJ.out.tab first_output/first.SJ.out.tab && touch first_output/first.SJ.out.tab
mv -vi first_output/Aligned.sortedByCoord.out.bam first_output/first.Aligned.sortedByCoord.out.bam && touch first_output/first.Aligned.sortedByCoord.out.bam
mv -vi first_output/first.* .
/bin/echo -ne "${READ_1}" > source_pair.txt
for F in `find first.*`; do
	mv -v ${F} ${READ_1}.${F} ; 
done ;
chmod -vR 777 first_output
sleep 5
'''

}



process star_align_second {



input:
	path '*' 
	path 'star_db/*'
	path 'gtf_file'
output:
	tuple file('source_pair.txt'),file('*.second.chimeric.out.junction'),file('*.second.SJ.out.tab'),file('*.second.Aligned.sortedByCoord.out.bam') 



shell:
'''
ls -alht
find .
READ_1=`find *_R1_*gz`
READ_2=`find *_R2_*gz`
STAR --runThreadN 8 \
	--genomeDir ${PWD}/star_db  \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn ${READ_2} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${PWD}/second_output/ \
	--outReadsUnmapped Fastx \
	--outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS \
	--outSJfilterOverhangMin 15 15 15 15 \
	--alignSJoverhangMin 15 \
	--sjdbGTFfile !{gtf_file} \
	--alignSJDBoverhangMin 15 \
	--outFilterMultimapNmax 20 \
	--outFilterScoreMin 1 \
	--outFilterMatchNmin 1 \
	--outFilterMismatchNmax 2 \
	--chimSegmentMin 15 \
	--chimScoreMin 15 \
	--chimScoreSeparation 10 \
	--chimJunctionOverhangMin 15 1>star.out 2>star.err
mv -vi second_output/Chimeric.out.junction second_output/second.chimeric.out.junction && touch second_output/second.chimeric.out.junction
mv -vi second_output/SJ.out.tab second_output/second.SJ.out.tab && touch second_output/second.SJ.out.tab
mv -vi second_output/Aligned.sortedByCoord.out.bam second_output/second.Aligned.sortedByCoord.out.bam && touch second_output/second.Aligned.sortedByCoord.out.bam
mv -vi second_output/second.* .
/bin/echo -ne "${READ_1}" > source_pair.txt
for F in `find second.*`; do
	mv -v ${F} ${READ_1}.${F} ; 
done ;
chmod -vR 777 second_output
sleep 5
'''

}







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




process circtest {

input:
	file '*' 
	file 'cohort_comp_conf'
output:
	file '*.circtest_results.zip'

shell:
'''
ls -alht
ls -alht /usr/local/bin/
##############################
#edit the files to have shorter names in the column headers; use the passed-in JSON file
for DCCF in LinearCount CircRNACount ; do
	python3 /usr/local/bin/name_edit.py !{cohort_comp_conf}  ${DCCF} s.${DCCF} ;
done ;
##############################
#distribute the columns into desired cohorts using
# an inline python3 script
COHORT_TO_RUN=`cat input.1` ; 
echo "COHORT_TO_RUN IS ${COHORT_TO_RUN}" ; 
COHORT_DIST="""
import os
import csv
import copy
import shutil
import json
linear_file='s.LinearCount'
circ_file='s.CircRNACount'
coord_file='CircCoordinates'

# read cohort/comparison data from JSON file
cohort_comp_conf_json_file='!{cohort_comp_conf}'
jreader=open(cohort_comp_conf_json_file)
cohort_comp_conf=json.load(jreader)
jreader.close()
cohort_comps=cohort_comp_conf['cohort_comps']
comp_names=cohort_comp_conf['comp_names']


def writeToFile(s,f):
    writer=open(f,'w')
    writer.write(s.strip()+'\\n')
    writer.close()
for comp in cohort_comps:
    if(not(comp=='${COHORT_TO_RUN}')):
        continue
    print(f'Processing for comparision {comp}')
    common_keys=['Chr','Start','End']
    os.makedirs(comp,exist_ok=True)
    shutil.copyfile(coord_file,comp+'/'+coord_file)
    tsv_files=[linear_file,circ_file]
    out_files=[comp+'/'+t for t in tsv_files]
    writeToFile(comp_names[comp][0],comp+'/G1.txt')
    writeToFile(comp_names[comp][1],comp+'/G2.txt')
    for file_pair in zip(tsv_files,out_files):
        with open(file_pair[0], mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file,delimiter='\\t')
            wroteHeader=False
            writer=open(file_pair[1],'w')
            recs_written=0
            for row in csv_reader:
                    if(not(wroteHeader)):
                        header=copy.deepcopy(common_keys)
                        header.extend(cohort_comps[comp])
                        wroteHeader=True
                        writer.write('\\t'.join(header)+'\\n')
                    row_keys_to_write=copy.deepcopy(common_keys)
                    row_keys_to_write.extend(cohort_comps[comp])
                    vals_to_write=[row[rk] for rk in row_keys_to_write]
                    writer.write('\\t'.join(vals_to_write)+'\\n')
                    recs_written=recs_written+1
            writer.close()

    #now write the counts to  analysis_col.txt
    counts_dict=cohort_comp_conf['comparison_counts']
    the_counts_this_comp=counts_dict[comp]
    val_to_write=the_counts_this_comp+3
    writer=open(comp+'/analysis_col.txt','w')
    writer.write(str(val_to_write))
    writer.close()    

"""
python3 -c "${COHORT_DIST}" ; 


BASE_DIR="${PWD}" ; 
for CT_DIR in `find */s.*|grep -Pi 'inear'|xargs -I {} dirname {}`; do
	cd ${BASE_DIR}/${CT_DIR} ; 
	DATA_IDX=`cat analysis_col.txt`; 
	NAME_A=`cat G1.txt`;
	NAME_B=`cat G2.txt`;
	#start a job
	set -x
	circtest.cli.R s.CircRNACount s.LinearCount  CircCoordinates ${DATA_IDX} ${NAME_A} ${NAME_B}
	set +x
done ;
#zip up the cohort results
cd ${BASE_DIR} ;
DG1=`find */G1.txt|xargs dirname`; 
zip -r "${COHORT_TO_RUN}.circtest_results.zip" "${DG1}" ; 
'''

}



process circtest_plain {

input:
	file '*'
	file 'cohort_comp_conf'
output:
	file '*.circtest_results.zip'

shell:
'''
ls -alht
ls -alht /usr/local/bin/
##############################
#edit the files to have shorter names in the column headers; use the passed-in JSON file
for DCCF in LinearCount CircRNACount ; do
	python3 /usr/local/bin/name_edit.py !{cohort_comp_conf}  ${DCCF} s.${DCCF} ;
done ;
##############################
#distribute the columns into desired cohorts using
# an inline python3 script
COHORT_TO_RUN=`cat input.1` ; 
echo "COHORT_TO_RUN IS ${COHORT_TO_RUN}" ; 
COHORT_DIST="""
import os
import csv
import copy
import shutil
import json
linear_file='s.LinearCount'
circ_file='s.CircRNACount'
coord_file='CircCoordinates'

# read cohort/comparison data from JSON file
cohort_comp_conf_json_file='!{cohort_comp_conf}'
jreader=open(cohort_comp_conf_json_file)
cohort_comp_conf=json.load(jreader)
jreader.close()
cohort_comps=cohort_comp_conf['cohort_comps']
comp_names=cohort_comp_conf['comp_names']


def writeToFile(s,f):
    writer=open(f,'w')
    writer.write(s.strip()+'\\n')
    writer.close()
for comp in cohort_comps:
    if(not(comp=='${COHORT_TO_RUN}')):
        continue
    print(f'Processing for comparision {comp}')
    common_keys=['Chr','Start','End']
    os.makedirs(comp,exist_ok=True)
    shutil.copyfile(coord_file,comp+'/'+coord_file)
    tsv_files=[linear_file,circ_file]
    out_files=[comp+'/'+t for t in tsv_files]
    writeToFile(comp_names[comp][0],comp+'/G1.txt')
    writeToFile(comp_names[comp][1],comp+'/G2.txt')
    for file_pair in zip(tsv_files,out_files):
        with open(file_pair[0], mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file,delimiter='\\t')
            wroteHeader=False
            writer=open(file_pair[1],'w')
            recs_written=0
            for row in csv_reader:
                    if(not(wroteHeader)):
                        header=copy.deepcopy(common_keys)
                        header.extend(cohort_comps[comp])
                        wroteHeader=True
                        writer.write('\\t'.join(header)+'\\n')
                    row_keys_to_write=copy.deepcopy(common_keys)
                    row_keys_to_write.extend(cohort_comps[comp])
                    vals_to_write=[row[rk] for rk in row_keys_to_write]
                    writer.write('\\t'.join(vals_to_write)+'\\n')
                    recs_written=recs_written+1
            writer.close()

    #now write the counts to  analysis_col.txt
    counts_dict=cohort_comp_conf['comparison_counts']
    the_counts_this_comp=counts_dict[comp]
    val_to_write=the_counts_this_comp+3
    writer=open(comp+'/analysis_col.txt','w')
    writer.write(str(val_to_write))
    writer.close()    

"""
python3 -c "${COHORT_DIST}" ; 


BASE_DIR="${PWD}" ; 
for CT_DIR in `find */s.*|grep -Pi 'inear'|xargs -I {} dirname {}`; do
	cd ${BASE_DIR}/${CT_DIR} ; 
	DATA_IDX=`cat analysis_col.txt`; 
	NAME_A=`cat G1.txt`;
	NAME_B=`cat G2.txt`;
	#start a job
	set -x
	circtest.cli.PLAIN.R s.CircRNACount s.LinearCount  CircCoordinates ${DATA_IDX} ${NAME_A} ${NAME_B}
	set +x
done ;
#zip up the cohort results
cd ${BASE_DIR} ;
DG1=`find */G1.txt|xargs dirname`; 
zip -r "${COHORT_TO_RUN}.circtest_results.zip" "${DG1}" ; 
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
	star_pair_tuple=star_align_pair(rrRNACleaned,star_ref,star_gtf)
	// b) just read 1
	star_first_tuple=star_align_first(rrRNACleaned,star_ref,star_gtf)
	// c) just read 2
    star_second_tuple=star_align_second(rrRNACleaned,star_ref,star_gtf)

	//combine the star results and merge with FQ files and input to DCC
	all_star_results=star_pair_tuple.mix(star_first_tuple)
		.mix(star_second_tuple)
		.map { [new File(""+it[0]).text , it[1],it[2],it[3] ] }
	DCC_input=rrRNACleaned.map { [ (new File(""+it[0]).getName())+"" , it[0],it[1]] }
		.combine(all_star_results)
		.filter { it[0].startsWith(it[3]) }
		.groupTuple()
		// pair name, R1 , R2, chimeric junction (3 per sample) , SJ.out(3 per sample) , Aligned.sortedByCoord.out.bam (3 per sample)
		.map{ [ it[0],it[1][0],it[2][0] , it[4][0],it[4][1],it[4][2] , it[5][0],it[5][1],it[5][2] , it[6][0],it[6][1],it[6][2]    ] }
		.collect()
	circtest_staging=DCC_step(
		DCC_input,
		star_gtf,
		channel.fromPath(params.fa_ref_file),
		channel.fromPath(params.repeat_file)
	)

	//prepare for circtest
	comp_channel=Channel
    	.fromPath(params.comp_list)
    	.splitText()
	comp_circtest_channel_for_split=comp_channel.combine(circtest_staging)
	//result=comp_circtest_channel_for_split.multiMap { it ->
    //    comp_circtest_channel: it
    //    comp_circtest_channel_plain: it
    //}
    circtest(comp_circtest_channel_for_split,channel.fromPath(params.cohort_comp_conf))
	circtest_plain(comp_circtest_channel_for_split,channel.fromPath(params.cohort_comp_conf))
    


}
