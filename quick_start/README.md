# CircRNAFlow Quick Start

This file offers some instructions for a "quick start" with CircRNAFlow to run it *locally*.  It assumes that you are running Linux, with *docker* installed, as well as the latest version of nextflow.  It also assume that the user has sufficient permission to run docker commands as well.

Follow these steps to run a demonstration analysis.

1. Download the sample/test data from [Zenodo](https://zenodo.org/) to **this** directory at this link : [https://zenodo.org/record/7339842#.Y3uML9LMI5m](https://zenodo.org/record/7339842#.Y3uML9LMI5m) .  Once the file is downloaded, uncompress it with this command :
```
tar -xvzf lumacaftor_small_test_data.circrnaflow.tar.gz
```
2. Build the docker containers using this command run from **this** directory: 

```
cd ../singularity && bash ./build_all_dockers.sh && cd ../quick_start
```
3. SymLink the CircRNAFlow pipline to this directory :
```
ln -vs `find ../*.nf`
```
4. Run the command below (also from this directory) :

```
./run_pipe_demo.sh
```

### Some notes about the test data

Inside the tar file the test data are under sample_data/*.gz .   The command above and the config file refer to it properly and should run as is.

The test data are subsetted to reads that align to a gene named ATP8B4.  Reads in the test data align to it and some circRNA there too.  The test data is small by design , so the pathway analysis (ClusterProfiler) should be considered with some skepticism!  The CircRNA ID is hsa-ATP8B4_0010 and is the CircATLAS ID :  http://159.226.67.237:8080/new/circ_detail.php?ID=hsa-ATP8B4_0010

An output is cp_outputs.zip in the ClusterProfiler step.  It is a zip file containing both input to an output from ClusterProfiler.  The input to it is cp_input.tsv ; adjacent to it (and also in the zip file) is full_annotated_cp_input.tsv which has additional information in a tabular format.  

Also, a premier output of the pipeline is a plot of "differential" expression  (circRNA relative to its gene).  The plot output is named as 15.49876278.49879459.ATP8B4.minus.pdf in the output directory of the plotting step (not the 'plain' plotting step).
