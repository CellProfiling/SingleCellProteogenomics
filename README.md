# SingleCellProteogenomics

This repository contains the code used to perform the single-cell proteogenomic analysis of the human cell cycle. This study was based on immunofluorescence staining of ~200k cells for single-cell analysis of proteomic heterogeneity and ~1k cells for analysis of single-cell RNA variability. These analyses were integrated with MS proteomic studies to investigate the functional importance of transcript-regulated and non-transcript regulated variability.

## Preprint
Please find the preprint manuscript here: https://www.biorxiv.org/content/10.1101/543231v2

## Structure of repository
This repository contains several analysis files. These are listed in order of execution, e.g. "1_", "2_" etc. The output of each script is used in the subsequent script.

The logic for these analyses is contained in the SingleCellProteogenomics python folder/module.

The input files are contained in "input." Output files are added to a folder "output" during the analysis, and figures are added to a folder "figures."

An R-script used to analyze skewness and kurtosis (noted in the Methods of the manuscript) is contained in the other_scripts folder.

## Single cell RNA-Seq analysis

For the `snakemake` workflow used to analyze the scRNA-Seq dataset, including RNA velocity calculations, please see this repository: ...