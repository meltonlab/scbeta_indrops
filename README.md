# Charting cellular identity during human in vitro β-cell differentiation

## What's this? 

This repository is a copy of the code used in the analyses performed in the "Charting cellular identity during human in vitro β-cell differentiation" publication (Veres *et. al*, Nature 2019). The code is provided here in the interest of reproducibility and transparency, and displays all data processing steps from the inDrops read data to the final figures in the paper. Although most steps can be executed on a personal device, some steps are best carried out on a high-performance environment and we provide example commands for our specific environment (Slurm). The code is organized in sections that follow the sections of the paper. As described in the methods accompanying the paper, much of the overall approach is derived from the Zeisel *et al.* (Cell 2018) publication and the accompanying repository and we make extensive use of Scanpy for pseudotime-related analyses.

## Downloading the data

The data is available for download from GEO (at this time, the GEO submission is being updated to the final version and will be publicly available shortly). Each subdirectory with a `notebooks` subdirectory, should also contain a `data` directory. These `data` components are available and packaged in the following tarball on GEO.

## Setting up the right environment.

The following steps create an environment that can be used to execute the notebooks.

```
```