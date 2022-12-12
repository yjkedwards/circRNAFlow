# CircRNAFlow

A pipeline for Circular RNA analysis

 * based on DCC, supported by STAR, bowtie, CircTest, ClusterProfiler, additional scripts, and of course Nextflow
 * requires .fastq.gz files to follow a naming convention

See the [Quick Start](quick_start/README.md) for a demonstration!

### Software Version Notes

| Software      | Language | Version/Commit | Parameter Notes |
| ------------- |  -------- | -------------- | --------------- |
| FastQC        | Java     | 0.11.9      | Default           |
| FlexBar        | C/C++     | 3.5.0      | -z GZ -m 30 -u 0 -q TAIL -qt 28 -a adapters.fasta         |
| Bowtie2        | C/C++     | bowtie2-2.3.0      | --no-unal --threads 8           |
| STAR        | C/C++    | 2.7.9a      | see paper appendix/suppl. materials           |
| DCC        | Python     | v0.5.0      | DCC --keep-temp -T $CPUCOUNT -D -G -mt1 @mate1.txt -mt2 @mate2.txt -R repeats.gtf -an annotation.gtf -Pi -F -M -Nr 2 3 -fg -B @bam_files -A ref.fasta @samplesheet.txt    |
| CircTest        | R     |  commit 0f5a86c    |  filtering run with various settings (see text); circ.test run with alpha=1.0            |

### Development Notes

The pipeline so har has been developed/tested with nextflow version 21.10.6.
