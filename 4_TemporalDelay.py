#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, Loaders, TemporalDelay
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

#%% Read in the protein data
import_dict = Loaders.load_temporal_delay()
u_well_plates, wp_ensg = import_dict["u_well_plates"], import_dict["wp_ensg"]
wp_iscell, wp_isnuc, wp_iscyto = import_dict["wp_iscell"], import_dict["wp_isnuc"], import_dict["wp_iscyto"]
ccd_comp, ccdtranscript = import_dict["ccd_comp"], import_dict["ccdtranscript"] 
pol_sort_well_plate, pol_sort_norm_rev = import_dict["pol_sort_well_plate"], import_dict["pol_sort_norm_rev"]
pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell = import_dict["pol_sort_ab_nuc"], import_dict["pol_sort_ab_cyto"], import_dict["pol_sort_ab_cell"], import_dict["pol_sort_mt_cell"]
var_comp_prot, gini_comp_prot, cv_comp_prot = import_dict["var_comp"], import_dict["gini_comp"], import_dict["cv_comp"]
var_cell_prot, gini_cell_prot, cv_cell_prot = import_dict["var_cell"], import_dict["gini_cell"], import_dict["cv_cell"]

#%% Idea: Make temporal heatmap for peak protein expression, and compare known and novel proteins that peak at similar times
# Execution: plt.imshow makes a heatmap if given a 2D array
# Output: heatmap; correlations of known/novel proteins
highlights = []#'ORC6','DUSP19','BUB1B','DPH2', 'FLI1']
highlights_ensg = []#'ORC6','DUSP19','BUB1B','DPH2', 'FLI1']

protein_heatmap_results = TemporalDelay.protein_heatmap(highlights, highlights_ensg, 
    ccd_comp, u_well_plates, wp_ensg, pol_sort_norm_rev, pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, wp_iscell, wp_isnuc, wp_iscyto)
sorted_maxpol_array, wp_binned_values, wp_max_pol, wp_max_pol_ccd, xvals = protein_heatmap_results

# Correlations of known and novel proteins that peak at similar times
TemporalDelay.peak_expression_correlation_analysis(wp_binned_values, wp_max_pol, wp_ensg, pol_sort_well_plate, u_well_plates)

#%% Create a heatmap of peak RNA expression
highlight_names, highlight_ensg = [],[]
adata, sorted_max_moving_avg_pol_ccd, norm_exp_sort, max_moving_avg_pol = TemporalDelay.rna_heatmap(highlight_names, highlight_ensg, ccdtranscript, xvals)

#%% Compare the variances and time of peak expression between protein and RNA
TemporalDelay.compare_variances_prot_v_rna(adata, norm_exp_sort, wp_ensg, var_comp_prot, gini_comp_prot, cv_comp_prot, var_cell_prot, gini_cell_prot, cv_cell_prot)
TemporalDelay.compare_peak_expression_prot_v_rna(adata, wp_ensg, ccd_comp, ccdtranscript, wp_max_pol, wp_max_pol_ccd, sorted_maxpol_array, max_moving_avg_pol, sorted_max_moving_avg_pol_ccd)

