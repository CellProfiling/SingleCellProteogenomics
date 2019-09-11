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

#%% Look at the percent variance that's explained by protein expression and by RNA expression
# Idea: It'd be nice to show there are genes that are highly regulated by protein and ones that are highly regulated
#     by RNA. And there will be some that aren't really regulated by either.
# Execution: Save the corresponding data from the antibody work in TemporalDelay
# Output: scatter of % var explained by protein vs % var explained by RNA

#%% 
# Idea: Making a circle plot linking protein and RNA peak expression to see where the chords line up. 
#      Are there any themes?
# Execution: ??? circos?
# Output: circle plot