# -*- coding: utf-8 -*-
"""
Analysis of protein abundance by cell cycle phase.
-  Uses Gaussian clustering to group cell expression into cell cycle phases.
-  This is referred to as the mock-bulk analysis in the paper.

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import ProteinDataPreparation, ProteinGaussianClustering, ProteinVariability, ProteinBimodality, ProteinCellCycleDependence
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
