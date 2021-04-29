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

The cell cycle phase and FACS intensity data for these ~1,000 cells are contained in the [input folder](https://drive.google.com/file/d/1mdQbYcDPqiTOHeiYbv_4RtrxrmlhYMNl/view?usp=sharing) within the file `ProteinData/WellPlatePhasesLogNormIntensities.csv`:
* The column "Well_Plate" (e.g., A10_355) corresponds to the sample title within GEO SRA (e.g., "Single U2OS cell A10_355"). 
* The column "Stage" corresponds to the phase assigned by FACS gating.
* The "Green530" and "Red585" columns correspond to the log-intensities for the red (CDT1) and green (GMNN) FUCCI markers for the individual cells.

The `snakemake` workflow used to analyze the scRNA-Seq dataset, including RNA velocity calculations and louvain unsupervised clustering, can be found in this repository: https://github.com/CellProfiling/FucciSingleCellSeqPipeline.

The `loom` file containing the results of RNA velocity analysis, including spliced and unspliced counts, can be found in the [input folder](https://drive.google.com/file/d/1mdQbYcDPqiTOHeiYbv_4RtrxrmlhYMNl/view?usp=sharing) under `RNAData/a.loom`, and the observation names used for each cell that match the "Well_Plate" identifiers can be found in `RNAData/a.obs_names.csv`.

## Citation

Mahdessian, D.\*; Cesnik, A. J.\*; Gnann, C.; Danielsson, F.; Stenström, L.; Arif, M.; Zhang, C.; Le, T.; Johansson, F.; Shutten, R.; Bäckström, A.; Axelsson, U.; Thul, P.; Cho, N. H.; Carja, O.; Uhlén, M.; Mardinoglu, A.; Stadler, C.; Lindskog, C.; Ayoglu, B.; Leonetti, M. D.; Pontén, F.; Sullivan, D. P.; Lundberg, E. “Spatiotemporal dissection of the cell cycle with single cell proteogenomics.” Nature, 2021, 590, 649–654. \*Contributed equally. https://www.nature.com/articles/s41586-021-03232-9
