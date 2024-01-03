## CircRNAFlow Quick Start for DSL2

🚧This is not supported and is under construction! See ../quick_start/README.md until this DSL2 is more supported.🚧


This "quick start" for the DSL2 version of the CircRNAflow pipline provides a step-by-step method along with test data to execute a test run of the pipeline.

So far the pipeline is developed and tested with *SLURM* and *local* executors with containers (either *docker* or *singularity*).  This quick-start is intended for using SLURM, but by changing the profile, docker containers or the local executor could be used!

### Quick Start Steps For Using the SLURM executor and singularity containers

1.  Download to this "quick_start_2" directory the test data from Zenodo using these three *links*:

* link1
* link2
* link3

These are referenced *here* and *here* on [Zenodo](https://zenodo.org/ "Zenodo").

2.   Unpack the test data to reveal directories for: a) sample data, b) pipeline/configuration data, c) reference data, and d) singularity image files.

```
find *.tar.gz | xargs -I {} tar -xvzf {}
```

Unpackad files should result in the file/directory structure as seen below with the 4 directories (pipe_data, sample_data, ref_data, and sif_images):
```
pipe_data
├── cohort_comp_conf.json
├── comp_list.txt
└── params.txt
ref_data
├── circatlas_human_bed_v2.0.txt
├── craft_input
│   ├── AGO2_binding_sites.bed
│   ├── hg38.00.idx
│   ├── hg38.01.idx
│   ├── hg38.02.idx
│   ├── hg38.fa
│   ├── hg38.fa.gz
│   ├── hg38.nhr
│   ├── hg38.nin
│   ├── hg38.nsq
│   ├── hg38.shd
│   ├── Homo_sapiens.GRCh38.104.gtf
│   └── Homo_sapiens.GRCh38.dna.primary_assembly.fa
├── hg38_rrna_prep
│   ├── human_rrnas.fasta
│   ├── human_rrnas.fasta.1.bt2
│   ├── human_rrnas.fasta.2.bt2
│   ├── human_rrnas.fasta.3.bt2
│   ├── human_rrnas.fasta.4.bt2
│   ├── human_rrnas.fasta.rev.1.bt2
│   └── human_rrnas.fasta.rev.2.bt2
├── Homo_sapiens.GRCh38.107.gtf
├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
├── human_repeats.gtf
├── mature.hsa.mirbase.fa
├── star_db
│   ├── chrLength.txt
│   ├── chrNameLength.txt
│   ├── chrName.txt
│   ├── chrStart.txt
│   ├── Genome
│   ├── genomeParameters.txt
│   ├── Log.out
│   ├── SA
│   └── SAindex
└── symbol_to_kegg_gid.HUMAN.csv
sample_data
├── adapters.fasta
├── SRR8383065_R1_f.fastq.gz
├── SRR8383065_R2_f.fastq.gz
├── SRR8383066_R1_f.fastq.gz
├── SRR8383066_R2_f.fastq.gz
├── SRR8383068_R1_f.fastq.gz
├── SRR8383068_R2_f.fastq.gz
├── SRR8383069_R1_f.fastq.gz
├── SRR8383069_R2_f.fastq.gz
├── SRR8383071_R1_f.fastq.gz
├── SRR8383071_R2_f.fastq.gz
├── SRR8383072_R1_f.fastq.gz
└── SRR8383072_R2_f.fastq.gz
sif_images
├── bowtieuab.sif
├── circtestuab.sif
├── clusterprofileruab.sif
├── craftuab.sif
├── ct_plotteruab.sif
├── dccuab.sif
├── fastqcuab.sif
├── flexbaruab.sif
├── kipoideeptargetuab.sif
└── staruab.sif

```


3. Create symlinks to the DSL2 code:
```
ln -vs ../modules
ln -vs ../circRNAFlow.DSL2.nf
```
4. Create an SBATCH file (say "sbatch_me.sh") which will run the pipeline.  The example below is available to be edited and customized:
```

```
**NOTE**: within the file you create, *be sure* to customize the directories/paths in the file so that they exist on your system and so that they point to this clone of the repo.

5. Customize run_pipe_demo.sh and config file as necessary.

Set the PROFILES variable to be "singularity,slurm".  The paths for inputs should properly resolve to data on disk if steps 1 and 2 above were carried out.  Otherwise, if the data above were downloaded, but in different areas, update the paths as necessary.

The config file (demo.config) in this directory, is set up to use singularity images pulling them from dockerhub.  If that is not desired and .sif files are preferred, update the config file to use .sif images under the 

The *SLURM* profile of the config file has been customized to use the [Cheaha HPC](https://www.uab.edu/it/home/research-computing/cheaha "CHEAHA") center.  For example, the queue names have been set to use queues there.  Queue names may need adjusting (e.g. "medium" changed or "express" changed to valid names for your SLURM installation.).

6. Submit the job!

```
sbatch sbatch_me.sh
```

7. Per the config file, output will appear in "quick_start_output"

