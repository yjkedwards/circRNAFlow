## CircRNAFlow Quick Start for DSL2

🚧This is not supported. See ../quick_start/README.md until this DSL2 is more supported.🚧


This "quick start" for the DSL2 version of the CircRNAflow pipline provides a step-by-step method along with test data to execute a test run of the pipeline.

So far the pipeline is developed and tested with *SLURM* and *local* executors with containers (either *docker* or *singularity*).  This quick-start is intended for using SLURM, but by changing the profile, docker containers or the local executor could be used!

### Quick Start Steps For Using the SLURM executor and singularity containers

1.  Download to this directory the test data from Zenodo using this *link* or this command:

```
wget LINK
```

2.   Unpack the test data to reveal directories for: a) sample data, b) pipeline/configuration data, c) reference data, and d) singularity image files.

```
tar -xvzf NAME
```
3. Create symlinks to the DSL2 code:
```
ln -vs ../modules
ln -vs ../circRNAFlow.DSL2.nf
```
4. Create an SBATCH file which will run the pipeline.  The example below is available to be edited and customized:
```

```
**NOTE**: *be sure* to customize the directories/paths in the file so that they exist on your system and so that they point to this clone of the repo.

5. sdfsadf
6. dfgsadfgsdafg
7. asfsdafsadf
8. asdfasdfsadf
9. asfsadfsdfa

