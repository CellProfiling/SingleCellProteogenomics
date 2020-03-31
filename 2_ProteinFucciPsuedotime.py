#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import FucciPseudotime, ProteinVariability, ProteinBimodality, CellCycleDependence
from SingleCellProteogenomics.Loaders import load_protein_fucci_pseudotime
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

#%% Read in the protein data
import_dict = load_protein_fucci_pseudotime()
u_plate, well_plate, well_plate_imgnb, u_well_plates = import_dict["u_plate"], import_dict["well_plate"], import_dict["well_plate_imgnb"], import_dict["u_well_plates"]
ab_nuc, ab_cyto, ab_cell, mt_cell = import_dict["ab_nuc"], import_dict["ab_cyto"], import_dict["ab_cell"], import_dict["mt_cell"]
area_cell, area_nuc =  import_dict["area_cell"], import_dict["area_nuc"]
wp_ensg, wp_ab = import_dict["wp_ensg"], import_dict["wp_ab"]
green_fucci, red_fucci = import_dict["green_fucci"], import_dict["red_fucci"]
log_green_fucci_zeroc, log_red_fucci_zeroc = import_dict["log_green_fucci_zeroc"], import_dict["log_red_fucci_zeroc"]
log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale = import_dict["log_green_fucci_zeroc_rescale"], import_dict["log_red_fucci_zeroc_rescale"]
wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp = import_dict["wp_comp_kruskal_gaussccd_adj"], import_dict["wp_pass_kruskal_gaussccd_bh_comp"]
fucci_data = import_dict["fucci_data"]
wp_iscell, wp_isnuc, wp_iscyto = import_dict["wp_iscell"], import_dict["wp_isnuc"], import_dict["wp_iscyto"]
        
#%% 
# Idea: Calculate the polar coordinates and other stuff
# Exec: Devin's calculations
# Output: fucci plot with polar coordinates

pseudotime_result = FucciPseudotime.pseudotime_protein(fucci_data, 
                           ab_nuc,ab_cyto,ab_cell,mt_cell,area_cell, area_nuc,
                           well_plate,well_plate_imgnb, log_red_fucci_zeroc_rescale,log_green_fucci_zeroc_rescale)
pol_sort_well_plate, pol_sort_norm_rev, pol_sort_well_plate_imgnb, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell, pol_sort_area_cell, pol_sort_area_nuc, pol_sort_fred, pol_sort_fgreen = pseudotime_result

#%% Calculate measures of variance of protein abundance in single cells
# Idea: Calculate measures of variance, and show them in plots
# Execution: Now that we already have the data filtered for variability, this is just descriptive.
# Output: scatters of antibody vs microtubule variances by different measures of variaibility

use_log = False # toggle for using log-transformed intensities; we decided to use natural intensities
calculate_variaton_result = ProteinVariability.calculate_variation(use_log, u_well_plates, wp_iscell, wp_isnuc, wp_iscyto, 
                                                pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, pol_sort_well_plate_imgnb)
mean_mean_comp, var_comp, gini_comp, cv_comp, var_cell, gini_cell, cv_cell, var_mt, gini_mt, cv_mt = calculate_variaton_result

# Compare variances for protein and microtubules, the internal control for each image
general_boxplot((var_comp, var_mt), ("Protein", "Microtubules"), "", "Variance", "", False, f"figures/ProteinMicrotubuleVariances.pdf")
general_boxplot((cv_comp, gini_mt), ("Protein", "Microtubules"), "", "CV", "", False, f"figures/ProteinMicrotubuleCVs.pdf")
general_boxplot((gini_comp, gini_mt), ("Protein", "Microtubules"), "", "Gini", "", False, f"figures/ProteinMicrotubuleGinis.pdf")
print(f"{scipy.stats.kruskal(var_comp, var_mt)[1]}: p-value for difference between protein and microtubule variances")
print(f"{scipy.stats.kruskal(cv_comp, gini_mt)[1]}: p-value for difference between protein and microtubule CVs")
print(f"{scipy.stats.kruskal(gini_comp, gini_mt)[1]}: p-value for difference between protein and microtubule Gini indices")

#%% Gaussian clustering to identify biomodal intensity distributions
bimodal_results = ProteinBimodality.identify_bimodal_intensity_distributions(u_well_plates, wp_ensg,
             pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
             wp_iscell, wp_isnuc, wp_iscyto)
wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs, wp_isbimodal_generally, wp_bimodal_fcmaxmin = bimodal_results

#%% Determine cell cycle dependence for each protein
use_log_ccd = False
do_remove_outliers = True
alphaa = 0.05

# Determine cell cycle dependence for proteins
ccd_results = CellCycleDependence.cell_cycle_dependence_protein(
        u_well_plates, wp_ensg, use_log_ccd, do_remove_outliers,
        pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
        pol_sort_area_cell, pol_sort_area_nuc,
        wp_iscell, wp_isnuc, wp_iscyto,
        wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs, wp_comp_kruskal_gaussccd_adj)
wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_ccd_unibimodal, wp_comp_ccd_gauss, perc_var_comp, mean_diff_from_rng, wp_comp_eq_percvar_adj, mean_diff_from_rng_clust1, wp_comp_eq_percvar_adj_clust1, mean_diff_from_rng_clust2, wp_comp_eq_percvar_adj_clust2, folder = ccd_results

# Move the temporal average plots to more informative places
CellCycleDependence.copy_mvavg_plots_protein(folder, wp_ensg, wp_comp_ccd_difffromrng, wp_isbimodal_fcpadj_pass, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_ccd_unibimodal, wp_comp_ccd_gauss)
CellCycleDependence.global_plots_protein(alphaa, u_well_plates, wp_ccd_unibimodal, perc_var_comp, mean_mean_comp, gini_comp, cv_comp, mean_diff_from_rng, wp_comp_eq_percvar_adj, wp_comp_kruskal_gaussccd_adj)
CellCycleDependence.analyze_ccd_variation_protein(
    folder, u_well_plates, wp_ensg, wp_ab, wp_iscell, wp_isnuc, wp_iscyto,
    wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, 
    var_comp, gini_comp, 
    mean_diff_from_rng, wp_comp_kruskal_gaussccd_adj, wp_comp_eq_percvar_adj, 
    mean_diff_from_rng_clust1, wp_comp_eq_percvar_adj_clust1, mean_diff_from_rng_clust2, wp_comp_eq_percvar_adj_clust2,
    wp_isbimodal_fcpadj_pass, wp_isbimodal_generally, wp_ccd_unibimodal, wp_bimodal_fcmaxmin, wp_comp_ccd_gauss)

