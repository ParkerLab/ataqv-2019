# ataqv analysis
This repository contains scripts to download data, run analyses, and recreate many of the figures from our ataqv manuscript.

## General organization

This has been tested on a Linux server. Pipelines are meant to be used with Snakemake and drmr (can be downloaded from the Parker Lab GitHub). The general organization of the repo is as follows:

1. `src` and `bin` directories should be self-explanatory.
2. `data` directory contains our raw data and other 'generic' data such as fasta files, bwa indices, etc. (some elements are distributed in the repository itself, some are created/downloaded in the steps outlined below).
3. `tn5` directory contains the processing and analyses for the ATAC-seq experiments that kept PCR cycles constant across Tn5 concentrations.
4. `cluster_density` directory contains the processing and analyses for the ATAC-seq experiments that did NOT keep PCR cycles constant across Tn5 concentrations.
5. The `figures` directory will be created as the `make` commands outlined below are run. If figures are recreated using the `make` commands outlined later, they will appear in here (along with some necessary intermediate files).


## Dependencies not included in repo
These tools must be in your $PATH:
1. cta (v. 0.1.2; can be downloaded from the Parker Lab GitHub)
2. ataqv (v 1.0.0)

As stated, many of the pipelines for this analysis utilize `snakemake` and `drmr` to submit jobs to a resource manager (we use SLURM). `drmr` can be downloaded from the Parker Lab GitHub.

## Setup
1. Install conda environment: conda env create -f environment.yml; source activate ataqv
2. Clone the repository
3. Set a few environmental variables (of course, you may want to add each of these `export` commands to your `.bashrc`):
```bash
export ATAQV_HOME='/path/to/repo'
export PICARD_JAR=/path/to/picard.jar
export PYTHONPATH="$PYTHONPATH:${ATAQV_HOME}/src"
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

## Re-creating figures

Once the analyses have been run, this can be done using `make knit`
