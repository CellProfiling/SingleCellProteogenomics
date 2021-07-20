# -*- coding: utf-8 -*-
"""
Analysis of transcript abundance in individual cells over cell division time.
-  Cell division time is measured with FUCCI markers and modeled in log-log space using polar coordinates.
-  RNA abundance was measured with single-cell RNA sequencing.
-  FUCCI marker intensities are measured for each individual cell with fluorescence assisted cell sorting (FACS)

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import (FucciPseudotime,
                                      RNACellCycleDependence,
                                      RNADataPreparation,
                                      RNAVelocity,
                                      utils)
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import shutil
import os

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

bioccd = np.genfromtxt(
    "input/ProteinData/BiologicallyDefinedCCD.txt", dtype="str"
)  # from mitotic structures
wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)
u_plates = ["355","356","357"]

# #%% Check for new RNA inputs
# newinputsfolder = "newinputs/RNAData/"
# if os.path.exists(newinputsfolder):
    # for file in os.listdir(newinputsfolder):
        # filepath=f"{newinputsfolder}{file}"
        # targetfile=f"input/RNAData/{os.path.basename(file)}"
        # targetfilepc=f"{targetfile}.protein_coding.csv"
        # if os.path.exists(targetfile): os.remove(targetfile)
        # if os.path.exists(targetfilepc): os.remove(targetfilepc)
        # shutil.copyfile(filepath, targetfile)

#%% Convert FACS intensities for FUCCI markers to pseudotime using the same polar coordinate methods as for protein
# Idea: Use the polar coordinate pseudotime calculations to calculate the pseudotime for each cell
# Execution: Adapt Devin's code for the cells sorted for RNA-Seq
# Output: Make log-log fucci intensity plots for the cells analyzed by RNA-Seq; Plot of all fucci pseudotimes; table of pseudotimes for each cell
adata, phases_filt = RNADataPreparation.read_counts_and_phases(
    "Counts", False, "protein_coding", u_plates
)
adata = RNADataPreparation.zero_center_fucci(adata)
FucciPseudotime.pseudotime_rna(adata)

#%% Single cell RNA-Seq data preparation and general analysis
RNADataPreparation.general_plots(u_plates)
RNADataPreparation.analyze_noncycling_cells(u_plates)
adata, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=False
)
adata = RNADataPreparation.zero_center_fucci(adata)
RNADataPreparation.plot_pca_for_batch_effect_analysis(adata, "BeforeRemovingNoncycling")

#%% Idea: Similar to mock-bulk analysis for proteins, we can evaluate each gene bundled by phase across cells
# Execution: Make boxplots of RNA expression by phase
# Output: boxplots for each gene
valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_plates
)
adata, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=True
)
adata = RNADataPreparation.zero_center_fucci(adata)
# RNADataPreparation.plot_markers_vs_reads(adata)
# RNADataPreparation.plot_pca_for_batch_effect_analysis(adata, "AfterRemovingNoncycling")
g1 = adata.obs["phase"] == "G1"
s = adata.obs["phase"] == "S-ph"
g2 = adata.obs["phase"] == "G2M"
do_make_boxplots = False
if do_make_boxplots:
    for iii, ensg in enumerate(adata.var_names):
        maxtpm = np.max(
            np.concatenate((adata.X[g1, iii], adata.X[s, iii], adata.X[g2, iii]))
        )
        RNACellCycleDependence.boxplot_result(
            adata.X[g1, iii] / maxtpm,
            adata.X[s, iii] / maxtpm,
            adata.X[g2, iii] / maxtpm,
            "figures/RNABoxplotByPhase",
            ensg,
        )

#%% Idea: Display general RNA expression patterns in single cells using UMAP dimensionality reduction, and display with FUCCI pseudotime overlayed
# FucciPseudotime.pseudotime_umap(adata)  # Generate a UMAP with the pseudotime overlayed

# We can also show that the cycle pattern remains when the curated CCD genes or CCD proteins are removed,
# demonstrating that there's still valuable information about cell cycling beyond what was called CCD
# RNADataPreparation.demonstrate_umap_cycle_without_ccd(adata)

# Read in the curated CCD genes / CCD proteins from the present work / Non-CCD genes from the present work; filter for genes that weren't filtered in QC of RNA-Seq
ccd_regev_filtered, ccd_filtered, nonccd_filtered = utils.ccd_gene_lists(adata)
adata_ccdprotein, adata_nonccdprotein, adata_regevccdgenes = RNADataPreparation.is_ccd(
    adata, wp_ensg, ccd_comp, nonccd_comp, bioccd, ccd_regev_filtered
)

# Generate plots with expression of genes overlayed
expression_data = adata.X
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:, None]).T

# # Log-log FUCCI plot with RNA expression overlayed
# RNACellCycleDependence.plot_expression_facs(
    # adata,
    # wp_ensg[np.isin(wp_ensg, adata.var_names)],
    # normalized_exp_data,
    # "figures/GeneExpressionFucci",
# )

# Cluster the expression into phases and analyze it that way
bulk_phase_tests = RNACellCycleDependence.analyze_ccd_variation_by_phase_rna(
    adata, normalized_exp_data, biotype_to_use
)
# RNACellCycleDependence.plot_expression_boxplots(adata, wp_ensg[np.isin(wp_ensg, adata.var_names)], bulk_phase_tests, "figures/GeneExpressionBoxplots")
# RNACellCycleDependence.plot_expression_umap(adata, wp_ensg[np.isin(wp_ensg, adata.var_names)], "figures/GeneExpressionUmap")

#%% Moving average calculations and randomization analysis for RNA
rna_ccd_analysis_results = RNACellCycleDependence.analyze_ccd_variation_by_mvavg_rna(
    adata,
    wp_ensg,
    ccd_comp,
    bioccd,
    adata_nonccdprotein,
    adata_regevccdgenes,
    biotype_to_use,
)
percent_ccd_variance, total_gini, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms, ccdtranscript, ccdprotein, mvpercs = (
    rna_ccd_analysis_results
)

# RNACellCycleDependence.plot_umap_ccd_cutoffs(adata, mean_diff_from_rng)
# RNACellCycleDependence.figures_ccd_analysis_rna(
    # adata,
    # percent_ccd_variance,
    # mean_diff_from_rng,
    # pass_meandiff,
    # eq_percvar_adj,
    # wp_ensg,
    # ccd_comp,
    # ccd_regev_filtered,
# )
# RNACellCycleDependence.plot_overall_and_ccd_variances(
    # adata,
    # biotype_to_use,
    # total_gini,
    # percent_ccd_variance,
    # pass_meandiff,
    # adata_ccdprotein,
    # adata_nonccdprotein,
    # adata_regevccdgenes,
# )
RNACellCycleDependence.make_plotting_dataframe(
    adata,
    ccdtranscript,
    wp_ensg,
    bioccd,
    norm_exp_sort[np.argsort(fucci_time_inds), :],
    mvavg_xvals,
    moving_averages,
    mvpercs,
)
# RNACellCycleDependence.compare_to_lasso_analysis(adata, ccdtranscript)
# RNACellCycleDependence.analyze_cnv_calls(adata, ccdtranscript)

#%% Moving average calculations and randomization analysis for the spike-in internal controls
adata_spikeins, phases_spikeins = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spike_ins=True, biotype_to_use="", u_plates=u_plates
)
sc.pp.filter_genes(adata_spikeins, min_cells=100)
print(f"data shape after filtering: {adata_spikeins.X.shape}")

RNACellCycleDependence.ccd_analysis_of_spikeins(adata_spikeins, perms)

#%% Analyze RNA velocity
valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_plates, load_velocities=True
)
adata_blobless, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=False
)
RNAVelocity.analyze_rna_velocity(adata_blobless, mean_diff_from_rng)

#%% Analyze isoforms
adata_isoform, ccdtranscript_isoform = RNACellCycleDependence.analyze_isoforms(
    adata, ccdtranscript, wp_ensg, ccd_comp, nonccd_comp, u_plates
)
RNACellCycleDependence.compare_genes_to_isoforms(
    adata,
    ccdprotein,
    ccdtranscript,
    adata_nonccdprotein,
    adata_isoform,
    ccdtranscript_isoform,
)
