

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
	star_align_pair(rrRNACleaned,star_ref,star_gtf)
	// b) just read 1
	star_align_first(rrRNACleaned,star_ref,star_gtf)
	// c) just read 2
    star_align_second(rrRNACleaned,star_ref,star_gtf)



}
