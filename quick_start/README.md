# CircRNAFlow Quick Start

This file offers some instructions for a "quick start" with CircRNAFlow.  It assumes that you are running Linux, with *docker* installed, as well as the latest version of nextflow.  It also assume that the user has sufficient permission to run docker commands as well.

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
