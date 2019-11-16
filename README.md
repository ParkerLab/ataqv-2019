# ataqv analysis
This repository contains scripts to run analyses and recreate figures from our ataqv manuscript.

## General organization

Pipelines were run on a Linux server. Pipelines are meant to be used with Snakemake and NextFlow (NextFlow v. 19.04.1). The general organization of the repo is as follows:

1. `src` and `bin` directories should be self-explanatory.
2. `data` directory contains our raw data and other 'generic' data such as fasta files, bwa indices, etc. (some elements are distributed in the repository itself, some are created/downloaded in the steps outlined below).
3. `tn5` directory contains the processing and analyses for the ATAC-seq experiments that kept PCR cycles constant across Tn5 concentrations.
4. `cluster_density` directory contains the processing and analyses for the ATAC-seq experiments that did NOT keep PCR cycles constant across Tn5 concentrations.
5. `public_survey` directory contains the processing and analyses for the survey of public ATAC-seq data.


## Dependencies not included in repo
These tools must be in your $PATH:
1. cta (v. 0.1.2; can be downloaded from the Parker Lab GitHub)
2. ataqv (v 1.1.0)

## Setup
1. Install conda environment: conda env create -f environment.yml; source activate ataqv
2. Clone the repository
3. Set a few environmental variables (of course, you may want to add each of these `export` commands to your `.bashrc`):
```bash
export ATAQV_HOME='/path/to/repo'
export PICARD_JAR=/path/to/picard.jar
```
4. Fetch the necessary submodules:

```bash
git submodule init
git submodule update
cd ${ATAQV_HOME}/src/ATACseq-Snakemake && git checkout 0f49ddf && cd $ATAQV_HOME
```

5. Prepare the `data` directory, containing all generic data (e.g., fasta files, bwa indices, etc). When the GEO repository becomes public, this will also download our ATAC-seq fastq files from GEO (for now, this will obviously not happen):
```bash
make data
```
This will run most commands necessary to set up the data/ directory immediately in a consecutive fashion; however, because bwa indices take several hours to put together, it will submit jobs to drmr for the index creation (job names BWAINDEX). You should wait for those jobs to finish before proceeding.

## Running analyses
Change into `tn5` and `cluster_density` and follow the READMEs there to process the original data.
To process the public data, change into public_survey and follow the README there.
