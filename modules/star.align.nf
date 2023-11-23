process star_align_mode {



input:
	path '*' 
	path 'star_db/*'
	path 'gtf_file'
	val mode
output:
	tuple file('source_pair.txt'),file('*.*.chimeric.out.junction'),file('*.*.SJ.out.tab'),file('*.*.Aligned.sortedByCoord.out.bam') 

shell:
'''
set -x
ls -alht
find .
READ_1=`find *_R1_*gz`
READ_2=`find *_R2_*gz`

if [ "!{mode}" == "pair" ] ; then
	IN_FASTQ="${READ_1} ${READ_2}" ; 
fi ;
if [ "!{mode}" == "first" ] ; then
	IN_FASTQ="${READ_1}" ; 
fi ;
if [ "!{mode}" == "second" ] ; then
	IN_FASTQ="${READ_2}" ; 
fi ;



STAR --runThreadN 8 \
	--genomeDir ${PWD}/star_db  \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn ${IN_FASTQ}  \
	--readFilesCommand zcat \
	--outFileNamePrefix ${PWD}/align_output/ \
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
mv -vi align_output/Chimeric.out.junction align_output/!{mode}.chimeric.out.junction && touch align_output/!{mode}.chimeric.out.junction
mv -vi align_output/SJ.out.tab align_output/!{mode}.SJ.out.tab && touch  align_output/!{mode}.SJ.out.tab
mv -vi align_output/Aligned.sortedByCoord.out.bam align_output/!{mode}.Aligned.sortedByCoord.out.bam  && touch align_output/!{mode}.Aligned.sortedByCoord.out.bam
mv -vi align_output/!{mode}.* .
/bin/echo -ne "${READ_1}" > source_pair.txt
for F in `find !{mode}.*`; do
	mv -v ${F} ${READ_1}.${F} ;
done ;
chmod -vR 777 align_output
sleep 5

'''

}
