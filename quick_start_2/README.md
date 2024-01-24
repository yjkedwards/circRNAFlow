## CircRNAFlow Quick Start for DSL2

⚠️🚧⚠️🚧This is not supported and is under construction! See ../quick_start/README.md until this DSL2 is more supported.🚧⚠️ 🚧⚠️



This "quick start" for the DSL2 version of the CircRNAflow pipline provides a step-by-step method along with test data to execute a test run of the pipeline.

So far the pipeline is developed and tested with *SLURM* and *local* executors with containers (either *docker* or *singularity*).  This quick-start is intended for using SLURM, but by changing the profile, docker containers or the local executor could be used!

Though the pipeline has been tested with singularity and docker containers and has been run successfully locally or on SLURM other environment might work, but we may not be able to support those.  Additionally, the pipeline has been successfully run with mouse and human data.  In theory other reference could work, but only mouse and human have been used.

## Quick Start Steps For Using the SLURM executor and singularity containers

#### 1.  Clone the repo and change to this directory

```
git clone https://github.com/yjkedwards/circRNAFlow.git
cd circRNAFlow/quick_start_2
```

#### 2.  Be sure that nextflow is installed and is in the path.  

Nextflow is available [here](https://www.nextflow.io/ "Nextflow").  The DSL2 version of the pipeline has been developed with version 23.10.0.

####  3.  Download test data.

Download tests data to this "quick_start_2" directory the test data from Zenodo using these three *links*:

* [lumacaftor_small_test_data.circrnaflow.tar.gz](https://zenodo.org/records/7339842/files/lumacaftor_small_test_data.circrnaflow.tar.gz) **27.1 GB**
* [lumacaftor_small_test_data_CRAFT_and_deeptarget_DSL2.tar.gz](https://zenodo.org/records/10449545/files/lumacaftor_small_test_data_CRAFT_and_deeptarget_DSL2.tar.gz) **6.8 GB**
* [lumacaftor_small_test_data_JUST_SIFs_DSL2.tar.gz](https://zenodo.org/records/10449545/files/lumacaftor_small_test_data_JUST_SIFs_DSL2.tar.gz) **12.3 GB**

These are referenced *here* and *here* on [Zenodo](https://zenodo.org/ "Zenodo").

#### 4.   Unpack the test data.

Unpack the test data to reveal directories for: a) sample data, b) pipeline/configuration data, c) reference data, and d) singularity image files.

```
find *.tar.gz | xargs -tI {} tar -xvzf {}
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
#### 5. Create symlinks.

Create symbolic links to the DSL2 code:
```
ln -vs ../modules
ln -vs ../circRNAFlow.DSL2.nf
```
#### 6. Create an SBATCH file.

Create an SBATCH file (say "sbatch_me.sh") which will run the pipeline.  An example is below is available to be  edited and customized.  Note that "##########" show some parts to take note of and possibly customize.  Some paths and modules will likely need to be customized.  The email address will need to be updated too.
```
#!/bin/bash
#SBATCH --job-name=circRNAFlow_demo
########## NOTE : change the output path here!
#SBATCH --output=/path/that/exists/where/a/log/can/be/written/job.out
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --mail-type=FAIL
########## change the email address here!
#SBATCH --mail-user=YOUR_EMAIL@SERVER
#SBATCH --mem=700GB
##########  change the path here where the script changes directory
cd /path/to/circRNAFlow/quick_start_2
########## Our HPC uses modules.  We load singularity to be able to use singularity containers.
module load Singularity
########## We load a Java module to be able to use singularity containers.  You may have a different 
##########    module of Java or an even more recent Java in your PATH already.
module load Java/11.0.2
########## I add nextflow to my path to be able to run it with the aforementioned Java!
PATH=${PATH}:/path/to/dir/of/Nextflow_23.10.0
########## Finally the demo script is run!
/bin/bash ./run_pipe_demo.sh

```
**NOTE**: within the file you create, *be sure* to customize the directories/paths in the file so that they exist on your system and so that they point to this clone of the repo.

#### 7. Customize run_pipe_demo.sh and the config file.

Customize run_pipe_demo.sh and the config file (demo.config) as necessary.

##### run_pipe_demo.sh

Set the PROFILES variable (in run_pipe_demo.sh) to be "singularity,slurm".  The paths in the nextflow run command for inputs should already properly resolve to data on disk if steps 1 and 2 above were carried out.  Otherwise, if the data above were downloaded, but in different areas, or saved with different names, then update the paths as necessary.

##### Config File (demo.config)

The config file (demo.config) in this directory, is set up to use singularity images pulling them from dockerhub.  If that is not desired and .sif files are preferred, update the config file to use .sif images under the "sif_images" directory.

The *SLURM* profile of the config file has been customized to use the [Cheaha HPC](https://www.uab.edu/it/home/research-computing/cheaha "CHEAHA") center.  For example, the queue names have been set to use queues there.  Queue names may need adjusting (e.g. "medium" changed or "express" changed to valid names for your SLURM installation environment!).

___

###### NOTES on container profiles available 

Several profiles for *containers* are available in the demo config file:
* docker : for using docker containers from dockerhub
* singularity : for using singularity containers after making them from pulled docker images from dockerhub
* singularity_local_sifs : for using singularity images from .sif files downloaded (see file downloads/links above).  ***NOTE*** if errors/problems arise and persist during download/expansion of singularity images from dockerhub downloads, then try this profile!

___

###### NOTES on executor profiles available 

Several profiles for *executors* are available in the demo config file:
* local : for running processes on the same computer which is running singularity.  

NOTE in this case, it is suggested to edit the *queueSize* option in the config file.  Less capable systems should use 1 or 2.  More capable systems (having additional CPUs/cores and RAM) can use higher values (4, 6, or more).   For our 256GB, 24-core server we use 6 here.
* slurm : for using slurm to launch jobs.  

As noted above, depending on the nodes/queues available, the queue names may need to be modified.  The slurm cluster we use (as users, not admins!) have queues available including amd-hdr100,express, & medium.  You may need to edit those names in the config file.  We have set values for "cpus" and "time" which we hope sensible and permissive but not too extreme, but "cpus" and "time" may also need editing for slurm.  We'd like to note too that DCC is the most resource intensive in our experience.  


####  8. Submit the job!

```
sbatch sbatch_me.sh
```

####  9. Look for output.

Per the config file, output will appear in "quick_start_output"

___


## **Notes on Using Your Own Data**

To submit your own data, the command line parameters will need to be modified.  Additionally, the config file will also need to be modified.  Use the .sh and the .config files in this directory as an example!  Additionally, for comparisons see the section below ("For Defining Comparisons") on how to create/edit the cohort_comp_conf and comp_list files.

#### For Input FastqGZ

For Input FastQ data, be sure to name the files of a pair the exact same, but "\_R1\_" for R1 and similarly "\_R2\_" for R2.  Additionally, the files sholdbe named like this : "SRR8383065_R1_f.fastq.gz" matching the pattern:
* alphanumeric (for sample name),
* \_R1\_ for either R1 or R2,
* f.fastq.gz as the suffix.

These file name restrictions are in place instead of a sample sheet.

#### For Defining Comparisons

For comparison configuration use the files "pipe_data/cohort_comp_conf.json" and "pipe_data/comp_list.txt" available in [lumacaftor_small_test_data.circrnaflow.tar.gz](https://zenodo.org/records/7339842/files/lumacaftor_small_test_data.circrnaflow.tar.gz).  Additionally, more comparisons can be made by adding to the lists and structures in those files.  The files comp_list_BIG_example.txt and cohort_comp_conf_BIG_example.json in this directory serve as additional examples.  Notably in the JSON file the following sections are used and described here:
* comments : an array of comment strings
* comparison_counts: a dict-structure of key-value (string -> int) pairs indicating for each sample cohort the number of samples in it
* cohort_comps: a dict-structure of string -> list-of-strings key-value pairs defining samples contained in each cohort
* comp_names: short for "comparison names", a dict of string -> list-of-strings key-value pairs defining short/abbreviated names used in some outputs of the pipeline
* file_samp_ren: a dict of string->string key value pairs assigning a sample name to each input fastq

#### On Large Reference Data

For reference data (e.g. rRNA database, STAR database, GTF file, repeat_file), please see examples using human data in [lumacaftor_small_test_data.circrnaflow.tar.gz](https://zenodo.org/records/7339842/files/lumacaftor_small_test_data.circrnaflow.tar.gz).  The repeat_file is a GTF of repeat regions.  Some small bioinformatics exercises may be necessary to generate this for other organisms.

#### On Small Reference Data

For the circatlas file see an example file in [lumacaftor_small_test_data.circrnaflow.tar.gz](https://zenodo.org/records/7339842/files/lumacaftor_small_test_data.circrnaflow.tar.gz).  Newer data may be available later at the CircAtlas website: [here](https://ngdc.cncb.ac.cn/circatlas/).

For KEGG, used in ClusterProfiler use the organism code "hsa" for human for example or "mmu" for mouse.  The example data here is human.  The kegg_db file is a CSV file of geneID to KEGG ID.  See the file in  [lumacaftor_small_test_data.circrnaflow.tar.gz](https://zenodo.org/records/7339842/files/lumacaftor_small_test_data.circrnaflow.tar.gz) for an example.

For CRAFT we setup data/files as the [CRAFT page](https://github.com/annadalmolin/CRAFT) suggests.    Editing of the params.txt may be desired to customize settings.

For deeptarget a copy of the mirbase FASTA was downloaded and is used.  
