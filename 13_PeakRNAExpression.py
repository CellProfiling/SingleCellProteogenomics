#%% Imports
from imports import *
import numpy as np
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
import collections
import fucci_plotting
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms" # so that the gene-specific results scales match for cross-gene comparisons
print("reading scRNA-Seq data")
biotype_to_use="protein_coding"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=False, biotype_to_use=biotype_to_use)
qc_filtering(adata, do_log_normalize=False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Make the peak RNA heatmap
# Idea: generate a heatmap similar to the one for protein data
# Execution: use the list of significantly differentially expressed genes to pluck the normalized intensities for those genes
#      and order them by the time of peak expression
# output: heatmap