# SingleCellProteogenomics

This repository contains the code used to perform the [single-cell proteogenomic analysis of the human cell cycle](https://www.nature.com/articles/s41586-021-03232-9). This study was based on immunofluorescence staining of ~200k cells for single-cell analysis of proteomic heterogeneity and ~1k cells for analysis of single-cell RNA variability. These analyses were integrated with other proteomic studies and databases to investigate the functional importance of transcript-regulated and non-transcript regulated variability.

## Structure of repository
The code is listed in order of execution, e.g. "1_", "2_" etc. The output of each script is used in the subsequent script.

The logic for these analyses is contained in the `SingleCellProteogenomics` folder.

The input files are contained in the "input" folder. This folder is linked [here](https://drive.google.com/file/d/1mdQbYcDPqiTOHeiYbv_4RtrxrmlhYMNl/view?usp=sharing) as a zip file, `input.zip`. Expand this folder within the base directory of this repository. If you are looking for the raw imaging proteomic dataset produced after filtering artifacts and such, that is located [here](https://drive.google.com/file/d/11vjsZV-nmzPpFmA7ShbfHzmbrk057b1V/view?usp=sharing).

The output files are added to a folder "output" during the analysis, and figures are added to a folder "figures."

An R-script used to analyze skewness and kurtosis (noted in the Methods of the manuscript) is contained in the other_scripts folder. The `other_scripts/ProteinDisorder.py` script utilizes [IUPRED2A](https://iupred2a.elte.hu/) and a [human UniProt](https://www.uniprot.org/proteomes/UP000005640) database.

## Prerequisites

Prerequisites are listed in `enviro.yaml` file. They can be installed using `conda`:

1. Install Miniconda from https://docs.conda.io/en/latest/miniconda.html

2. From this directory run `conda env create -n fucci -f enviro.yaml; conda activate fucci` to install and activate these prerequisites.

## Single-cell RNA-Seq analysis

The single-cell RNA-Seq data is available at GEO SRA under project number [GSE146773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773). 

The `snakemake` workflow used to analyze the scRNA-Seq dataset, including RNA velocity calculations and louvain unsupervised clustering, can be found in this repository: https://github.com/CellProfiling/FucciSingleCellSeqPipeline.

## Citation

Mahdessian, D.\*; Cesnik, A. J.\*; Gnann, C.; Danielsson, F.; Stenström, L.; Arif, M.; Zhang, C.; Le, T.; Johansson, F.; Shutten, R.; Bäckström, A.; Axelsson, U.; Thul, P.; Cho, N. H.; Carja, O.; Uhlén, M.; Mardinoglu, A.; Stadler, C.; Lindskog, C.; Ayoglu, B.; Leonetti, M. D.; Pontén, F.; Sullivan, D. P.; Lundberg, E. “Spatiotemporal dissection of the cell cycle with single cell proteogenomics.” Nature, 2021, 590, 649–654. \*Contributed equally. https://www.nature.com/articles/s41586-021-03232-9
