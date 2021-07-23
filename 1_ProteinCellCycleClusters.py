# -*- coding: utf-8 -*-
"""
Single cell proteogenomics pipeline

@author: Anthony J. Cesnik, cesnik [at] kth.se
"""

import argparse
from SingleCellProteogenomics import (
    ProteinDataPreparation,
    ProteinGaussianClustering,
    ProteinCellCycleDependence,
    ProteinCellCycleDependenceAnalysis,
)

description = "Single cell proteogenomic analysis scripts"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--all', action='store_true', help="Perform plotting and side analyses not needed for HPA releases.")
parser.add_argument('--outputFilteredDataframe', action='store_true', help="Output the protein imaging post-processing data after filtering artifacts and such.")
parser.add_argument('--useLogProteinIntensities', action='store_true', help="Toggle for using log-transformed intensities. We decided to use natural intensities instead.")
parser.add_argument('--keepOutlierIntensities', action='store_true', help="Keep the intensities that are 5 SD from the mean; we decided to remove these.")
args = parser.parse_args()

do_plotting = args.all

"""
Protein Cell Cycle Clusters

Analysis of protein abundance by cell cycle phase.
-  Uses Gaussian clustering to group cell expression into cell cycle phases.
-  This is referred to as the "mock-bulk" analysis in the paper.
"""

ALPHA_GAUSS = 0.05

# Read in the raw protein imaging data
protein_data = ProteinDataPreparation.ProteinData(do_plotting)
protein_data.read_raw_data()
protein_data.read_sample_info(protein_data.raw_protein_data)
protein_data.previous_results()

# Filter the raw data using manual annotations and nucleus size
protein_data.apply_manual_filtering()
protein_data.apply_big_nucleus_filter()
protein_data.apply_cell_count_filter()
if args.outputFilteredDataframe:
    protein_data.filtered_protein_data.to_csv("output/filtered_protein_df.csv.gz")
protein_data.read_sample_info(protein_data.filtered_protein_data)
protein_data.previous_results()

# Filter for variation and get compartment information
protein_data.apply_variation_filter(output_dataframes=args.outputFilteredDataframe)

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
if do_plotting:
    protein_clustering.fucci_plots()

# Gaussian clustering per plate to identify G1/S/G2 and do kruskal test for variance
# Exec: sklearn.mixture.GaussianMixture & scipy.stats.kruskal
# Output: FDR for cell cycle variation per well per compartment
# NB! The cluster labels can change if any prior analysis changes. Inspect the plots so that top-left FUCCI cluster is G1, top-right is S, bottom-right is G2.
protein_clustering.make_cluster_labels(g1_idx=1, sph_idx=2, g2_idx=0)
protein_clustering.gaussian_clustering_analysis(
    alpha_gauss=ALPHA_GAUSS, do_plotting=do_plotting
)
if args.all:
    protein_clustering.address_replicates(alpha_gauss=ALPHA_GAUSS)

"""
Protein Fucci Pseudotime

Analysis of protein abundance in individual cells over cell division time.
-  Cell division time is measured with FUCCI markers and modeled in log-log space using polar coordinates.
-  The cell division time modeling is referred to as the "pseudotime" analysis in the paper.
-  Protein abundance is measured in individual asynchronous cells using immunofluorescence and antibody staining.
"""

ALPHA_PROTEIN_FUCCI = 0.01

# Calculate polar coordinates and fucci pseudotime
protein_ccd = ProteinCellCycleDependence.ProteinCellCycleDependence(protein_data, protein_clustering, do_plotting)
protein_ccd.fucci_polar_coords.pseudotime_protein(protein_ccd)

# Calculate measures of variance of protein abundance in single cells
protein_ccd.protein_variability.calculate_variation(protein_ccd, args.useLogProteinIntensities)

# Gaussian clustering to identify biomodal intensity distributions
protein_ccd.protein_bimodality.identify_bimodal_intensity_distributions(protein_ccd)

# Determine cell cycle dependence for each protein
protein_ccd.calculate_cell_cycle_dependence_protein(ALPHA_PROTEIN_FUCCI, args.useLogProteinIntensities, args.keepOutlierIntensities)
protein_ccd.analyze_ccd_variation_protein()

# Make a dataframe for plotting on the HPA website
ProteinCellCycleDependenceAnalysis.make_plotting_dataframe(protein_ccd)

if args.all:
    ProteinCellCycleDependenceAnalysis.additional_ccd_analysis(protein_ccd)
    ProteinCellCycleDependenceAnalysis.global_plots_protein(protein_ccd, ALPHA_PROTEIN_FUCCI)

    # Perform comparison to LASSO for finding CCD proteins
    ProteinCellCycleDependenceAnalysis.compare_to_lasso_analysis(protein_ccd)
    
    # Generate UMAPs to illustrate cutoffs and stability
    ProteinCellCycleDependenceAnalysis.generate_protein_umaps(protein_ccd)
