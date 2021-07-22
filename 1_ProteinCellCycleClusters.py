# -*- coding: utf-8 -*-
"""
Analysis of protein abundance by cell cycle phase.
-  Uses Gaussian clustering to group cell expression into cell cycle phases.
-  This is referred to as the "mock-bulk" analysis in the paper.

@author: Anthony J. Cesnik, cesnik [at] kth.se
"""

from SingleCellProteogenomics import (ProteinDataPreparation,
                                      ProteinGaussianClustering)

# Read in the raw protein imaging data
protein_data = ProteinDataPreparation.ProteinData()
protein_data.read_raw_data()
protein_data.read_sample_info(protein_data.raw_protein_data)
protein_data.previous_results()

# Filter the raw data using manual annotations and nucleus size
protein_data.apply_manual_filtering()
protein_data.apply_big_nucleus_filter()
protein_data.apply_cell_count_filter()
protein_data.filtered_protein_data.to_csv("output/filtered_protein_df.csv.gz")
protein_data.read_sample_info(protein_data.filtered_protein_data)
protein_data.previous_results()

# Filter for variation and get compartment information
protein_data.apply_variation_filter(output_dataframes=False)

# Filter out missing compartment information; these are localized to mitotic structures and handled differently
# Then demonstrate there is no more missing compartment information
protein_data.read_sample_info(protein_data.my_df_filtered_variation)
protein_data.metacompartments(protein_data.my_df_filtered_variation)
protein_data.read_sample_info(protein_data.my_df_filtered_compartmentvariation)
protein_data.metacompartments(protein_data.my_df_filtered_compartmentvariation)
protein_data.previous_results()

# Get and process single cell proteomic intensities
protein_data.read_sample_data(protein_data.my_df_filtered_compartmentvariation)
protein_clustering = ProteinGaussianClustering.ProteinGaussianClustering(protein_data)
protein_clustering.zero_center_fucci()
protein_clustering.fucci_plots()

# Gaussian clustering per plate to identify G1/S/G2 and do kruskal test for variance
# Exec: sklearn.mixture.GaussianMixture & scipy.stats.kruskal
# Output: FDR for cell cycle variation per well per compartment
# NB! The cluster labels can change if any prior analysis changes. Inspect the plots so that top-left FUCCI cluster is G1, top-right is S, bottom-right is G2.
protein_clustering.make_cluster_labels(
    g1_idx=1, 
    sph_idx=2, 
    g2_idx=0
)
protein_clustering.gaussian_clustering_analysis(
    alpha_gauss=0.05, 
    doGenerateBoxplotsPerGene=False
)
protein_clustering.address_replicates(
    alpha_gauss=0.05
)

