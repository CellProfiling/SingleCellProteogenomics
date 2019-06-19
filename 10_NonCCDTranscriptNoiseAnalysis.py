#%% Imports
from imports import *
import numpy as np

#%% Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
count_or_rpkm = "Rpkms" # so that the results match for cross-gene comparisons
adata, phases_filt = read_counts_and_phases(dd, count_or_rpkm)
qc_filtering(adata, False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Analyze the profiles of gene expression for whether the noise is Gaussian
# Idea: SmartSeq2 noise is Gaussian, so could we look at the expression levels
#      of the non-CCD proteins and 