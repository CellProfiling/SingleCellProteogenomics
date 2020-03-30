#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import ProteinDataPreparation, ProteinGaussianClustering, ProteinVariability, ProteinFucciPseudotime, ProteinBimodality
from SingleCellProteogenomics import CellCycleDependence
import numpy as np
import sklearn.mixture
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

#%% Read in the protein data (methods in ProteinDataPreparation.py)
my_df = ProteinDataPreparation.read_raw_data()
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = ProteinDataPreparation.read_sample_info(my_df)
wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = ProteinDataPreparation.previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% Idea: Filter the raw data (methods in ProteinDataPreparation.py)
# Execution: Use manual annotations and nucleus size to filter samples and images
# Output: Filtered dataframe
my_df_filtered = ProteinDataPreparation.apply_manual_filtering(my_df, result_dict, ab_dict)
my_df_filtered = ProteinDataPreparation.apply_big_nucleus_filter(my_df_filtered)
my_df_filtered = ProteinDataPreparation.apply_cell_count_filter(my_df_filtered)
my_df_filtered.to_csv("input/processed/python/nuc_predicted_prob_phases_filtered.csv")
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = ProteinDataPreparation.read_sample_info(my_df_filtered)
wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = ProteinDataPreparation.previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% Idea: Filter for variation and get compartments (methods in ProteinDataPreparation.py)
# Execution: Use annotated variation and compartment information
# Output: Number of cells filtered
my_df_filtered_variation, my_df_filtered_novariation = ProteinDataPreparation.apply_variation_filter(my_df_filtered, result_dict, my_df)

## Uncomment to output these dataframes (used for skewness / kurtosis analysis)
# my_df_filtered_variation.to_csv("output/nuc_predicted_prob_phases_filtered_variation.csv")
# my_df_filtered_novariation.to_csv("output/nuc_predicted_prob_phases_filtered_novariation.csv")

# filter out the ones missing compartment information; these are localized to mitotic structures and handled differently
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = ProteinDataPreparation.read_sample_info(my_df_filtered_variation)
wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation = ProteinDataPreparation.metacompartments(u_well_plates, compartment_dict, my_df_filtered_variation)

# demonstrate that there are no more missing compartment information
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = ProteinDataPreparation.read_sample_info(my_df_filtered_compartmentvariation)
wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation = ProteinDataPreparation.metacompartments(u_well_plates, compartment_dict, my_df_filtered_compartmentvariation)
wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = ProteinDataPreparation.previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% Idea: Get and process intensities (methods in ProteinDataPreparation.py and ProteinGaussianClustering.py)
# Execution: get intensities; zero center fucci intensities
# Output: Fucci plot
ab_nuc, ab_cyto, ab_cell, mt_cell, green_fucci, red_fucci = ProteinDataPreparation.read_sample_data(my_df_filtered_compartmentvariation)
log_green_fucci, log_red_fucci, log_green_fucci_zeroc, log_red_fucci_zeroc, log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale, fucci_data = ProteinGaussianClustering.zero_center_fucci(green_fucci, red_fucci, u_plate, well_plate, plate)
    
plt.hist2d(log_green_fucci_zeroc_rescale,log_red_fucci_zeroc_rescale,bins=200)
plt.xlabel("Log10 Green Fucci Intensity")
plt.ylabel("Log10 Red Fucci Intensity")
plt.savefig("figures/FucciPlotProteinIFData_unfiltered.png")
plt.show()
plt.close()

# General picture of antibody intensity density
sbn.distplot(ab_cell, hist=False)
plt.xlabel("Mean Intensity")
plt.ylabel("Density")
plt.savefig("figures/antibody_cell_intensity.pdf")
plt.show()
plt.close()

#%% Idea: Gaussian clustering per plate to identify G1/S/G2 and do kruskal test for variance
# Exec: sklearn.mixture.GaussianMixture & scipy.stats.kruskal
# Output: FDR for cell cycle variation per well per compartment
cluster_labels = ProteinGaussianClustering.gaussian_clustering(log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale)

# NB! The cluster labels can change if any prior analysis changes. Inspect the plots so that top-left FUCCI cluster is G1, top-right is S, bottom-right is G2.
g1, sph, g2 = cluster_labels == 2, cluster_labels == 1, cluster_labels == 0
alpha_gauss, doGenerateBoxplotsPerGene = 0.05, True
wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp, wp_mt_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_mt =  ProteinGaussianClustering.gaussian_clustering_analysis(alpha_gauss, doGenerateBoxplotsPerGene, g1, sph, g2, wp_ensg, well_plate, u_well_plates, ab_cell, ab_nuc, ab_cyto, mt_cell, wp_iscell, wp_isnuc, wp_iscyto)

# General look at replicates in mock-bulk analysis
ProteinGaussianClustering.address_replicates(alpha_gauss, wp_pass_kruskal_gaussccd_bh_comp, wp_ensg, wp_ab, u_well_plates)

#%% 
# Idea: Calculate the polar coordinates and other stuff
# Exec: Devin's calculations
# Output: fucci plot with polar coordinates

ProteinFucciPseudotime.fucci_polar_coordinate_calculations(fucci_data, 
                           ab_nuc,ab_cyto,ab_cell,mt_cell,area_cell, area_nuc,well_plate,well_plate_imgnb, log_red_fucci_zeroc_rescale,log_green_fucci_zeroc_rescale)

#%% Calculate measures of variance of protein abundance in single cells
# Idea: Calculate measures of variance, and show them in plots
# Execution: Now that we already have the data filtered for variability, this is just descriptive.
# Output: scatters of antibody vs microtubule variances by different measures of variaibility

use_log = False # toggle for using log-transformed intensities; we decided to use natural intensities
mean_mean_comp, var_comp, gini_comp, cv_comp, var_cell, gini_cell, cv_cell = ProteinVariability.calculate_variation(use_log, 
                                                u_well_plates, wp_iscell, wp_isnuc, wp_iscyto, 
                                                pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, pol_sort_well_plate_imgnb)

# Compare variances for protein and microtubules, the internal control for each image
general_boxplot((var_comp, var_mt), ("Protein", "Microtubules"), "", "Variance", "", False, f"figures/ProteinMicrotubuleVariances.pdf")
general_boxplot((cv_comp, gini_mt), ("Protein", "Microtubules"), "", "CV", "", False, f"figures/ProteinMicrotubuleCVs.pdf")
general_boxplot((gini_comp, gini_mt), ("Protein", "Microtubules"), "", "Gini", "", False, f"figures/ProteinMicrotubuleGinis.pdf")
print(f"{scipy.stats.kruskal(var_comp, var_mt)[1]}: p-value for difference between protein and microtubule variances")
print(f"{scipy.stats.kruskal(cv_comp, gini_mt)[1]}: p-value for difference between protein and microtubule CVs")
print(f"{scipy.stats.kruskal(gini_comp, gini_mt)[1]}: p-value for difference between protein and microtubule Gini indices")

#%% read in reliability scores
# ab_scores = list(np.zeros(wp_ab.shape, dtype=str))
# with open("input/processed/manual/reliabilityscores.txt") as file:
#     for line in file:
#         if line.startswith("Antibody RRID"): continue
#         score = line.split('\t')[1].strip()
#         ablist = line.split('\t')[0].replace(":","").replace(",","").split()
#         for ab in ablist:
#             if ab in wp_ab:
#                 ab_scores[wp_ab_list.index(ab)] = score
# compartmentstring = np.array(["Cell"] * len(wp_iscell))
# compartmentstring[wp_iscyto] = "Cyto" 
# compartmentstring[wp_isnuc] = "Nuc"
# pd.DataFrame({"ENSG":wp_ensg,"ab":wp_ab, "ab_hpa_scores": ab_scores, "compartment": compartmentstring}).to_csv("output/antibody_list.csv",index=False)

#%% Gaussian clustering to identify biomodal intensity distributions
wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs = ProteinBimodality.identify_bimodal_intensity_distributions(u_well_plates, 
             pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
             wp_iscell, wp_isnuc, wp_iscyto)

#%% Determine cell cycle dependence for each protein
use_log_ccd = False
do_remove_outliers = True
alphaa = 0.05

# Determine cell cycle dependence for proteins
wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_comp_ccd_gauss, folder = CellCycleDependence.cell_cycle_dependence_protein(
                    use_log_ccd, do_remove_outliers,
                    pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
                    wp_iscell, wp_isnuc, wp_iscyto,
                    wp_isbimodal_fcpadj_pass)

# Move the temporal average plots to more informative places
CellCycleDependence.copy_mvavg_plots_protein(folder, wp_comp_ccd_difffromrng, wp_isbimodal_fcpadj_pass, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_comp_ccd_gauss)
CellCycleDependence.global_plots_protein(alphaa, u_well_plates, perc_var_comp, mean_mean_comp, mean_diff_from_rng, wp_comp_eq_percvar_adj)
CellCycleDependence.analyze_ccd_variation_protein(u_well_plates, wp_ensg, wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_ccd_unibimodal)

ccdstring = np.array(["No                 "] * len(ccd_comp))
ccdstring[ccd_comp] = "Pseudotime"
ccdstring[np.isin(wp_ensg, bioccd)] = "Mitotic"
ccdstring[ccd_comp & np.isin(wp_ensg, bioccd)] = "Pseudotime&Mitotic"
pd.DataFrame({
    "well_plate" : u_well_plates, 
    "ENSG": wp_ensg,
    "antibody": wp_ab,
    "variance_comp":var_comp,
    "gini_comp":gini_comp,
    "known_by_GoReactomeCyclebaseNcbi":np.isin(wp_ensg, np.concatenate((knownccd1, knownccd2, knownccd3))),
    "mean_percvar_diff_from_random":mean_diff_from_rng,
    "wp_comp_kruskal_gaussccd_adj":wp_comp_kruskal_gaussccd_adj,
    "log_10_pval_eq_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj, wp_comp_eq_percvar_adj+1)),
    "pass_median_diff":wp_comp_ccd_difffromrng,
    "pass_gauss":wp_comp_ccd_gauss,
    "CCD_COMP":ccd_comp,
    "ccd_reason":ccdstring,
    "nonccd_comp":nonccd_comp,
    "wp_prev_ccd":wp_prev_ccd,
    
    # bimodal significance testing
    "ccd_unimodal":wp_comp_ccd_difffromrng,
    "ccd_clust1":wp_comp_ccd_clust1,
    "clust1_difffromrng":mean_diff_from_rng_clust1,
    "clust1_log10pval_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj_clust1, wp_comp_eq_percvar_adj_clust1+1)),
    "ccd_clust2":wp_comp_ccd_clust2,
    "ccd_clust2_difffromrng":mean_diff_from_rng_clust2,
    "ccd_clust2_log10pval_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj_clust2, wp_comp_eq_percvar_adj_clust2+1)),
    }).to_csv("output/CellCycleVariationSummary.csv", index=False)

#%%