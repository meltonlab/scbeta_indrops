# Charting cellular identity during human in vitro β-cell differentiation

## What's this? 

This repository is a copy of the code used in the analyses performed in the "Charting cellular identity during human in vitro β-cell differentiation" publication (Veres *et. al*, Nature 2019). The code is provided here in the interest of reproducibility and transparency, and displays all data processing steps from the inDrops read data to the final figures in the paper. Although most steps can be executed on a personal device, some steps are best carried out on a high-performance environment and we provide example commands for our specific environment (Slurm). The code is organized in sections that follow the sections of the paper. As described in the methods accompanying the paper, much of the overall approach is derived from the Zeisel *et al.* (Cell 2018) publication and the accompanying repository and we make extensive use of Scanpy for pseudotime-related analyses.

## Downloading the data

The data is available for download from GEO (at this time, the GEO submission is being updated to the final version and will be publicly available shortly). Each subdirectory with a `notebooks` subdirectory, should also contain a `data` directory. These `data` components are available and packaged in the following tarball on GEO.

### Data structure

We provide on GEO raw FASTQ files, raw post-alignment counts, processed counts data (with associated metadata) and analysis results files.

#### Raw FASTQ
These are the sequencing reads directly from the sequencer, demultiplexed according to the appropriate libraries. Processing requires using the inDrops pipeline or related tools (such as DropEst, etc).

#### Post-alignment counts
This is the output of the inDrops pipeline (`GSE114412_*.raw_indrops_counts.tsv.gz`), run with default settings and aligned using the GRCh38.88 genome and transcriptome annotations. These counts are further processed and filtered during analysis, as detailed in the online methods of the paper and in the notebooks and code of this repository. 

#### Processed counts data
These are the final filtered counts used `GSE114412_*.processed_counts.tsv.gz`) throughout the analyses presented in Veres *et al.* 2019, provided separately from the raw reads and counts to simplify re-analysis and exploration work. The associated metadata files include tSNE projection coordinates, cell type cluster labels and other information, such as pseudotime values for each dataset. 

#### Complete record of intermediate analysis steps
The `GSE114412_Veres2019_analysis_data.tar` tarball on GEO contains a `data` directory to match each `notebooks` directory in the `01_` to `04_` directories of this repository. They contain `.loom` or `.df.npz` files which can be read by the notebooks in this repository for further access. `05_Figures` contains the code used to produced the final figures of the paper. 

## Setting up the right environment.

The following steps create an environment that can be used to execute the notebooks.

```
```