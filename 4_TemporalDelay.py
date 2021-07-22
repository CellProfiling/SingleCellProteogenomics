# -*- coding: utf-8 -*-
"""
Comparison of the times of peak expression for protein and RNA for each gene
-  The peak expression for each protein and transcript were determined using the FUCCI pseudotime analysis
-  This is the first demonstration of the temporal delay between protein and RNA on the single cell level

@author: Anthony J. Cesnik, cesnik [at] kth.se
"""

from SingleCellProteogenomics import (Loaders, RNADataPreparation,
                                      TemporalDelay)
import matplotlib.pyplot as plt

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

#%% Read in the protein data
import_dict = Loaders.load_temporal_delay()
u_well_plates, wp_ensg = import_dict["u_well_plates"], import_dict["wp_ensg"]
wp_iscell, wp_isnuc, wp_iscyto = (
    import_dict["wp_iscell"],
    import_dict["wp_isnuc"],
    import_dict["wp_iscyto"],
)
ccd_comp, ccdtranscript, ccdtranscript_isoform = (
    import_dict["ccd_comp"],
    import_dict["ccdtranscript"],
    import_dict["ccdtranscript_isoform"],
)
pol_sort_well_plate, pol_sort_norm_rev = (
    import_dict["pol_sort_well_plate"],
    import_dict["pol_sort_norm_rev"],
)
pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell = (
    import_dict["pol_sort_ab_nuc"],
    import_dict["pol_sort_ab_cyto"],
    import_dict["pol_sort_ab_cell"],
    import_dict["pol_sort_mt_cell"],
)
var_comp_prot, gini_comp_prot, cv_comp_prot = (
    import_dict["var_comp"],
    import_dict["gini_comp"],
    import_dict["cv_comp"],
)
var_cell_prot, gini_cell_prot, cv_cell_prot = (
    import_dict["var_cell"],
    import_dict["gini_cell"],
    import_dict["cv_cell"],
)

#%% Idea: Make temporal heatmap for peak protein expression, and compare known and novel proteins that peak at similar times
# Execution: plt.imshow makes a heatmap if given a 2D array
# Output: heatmap; correlations of known/novel proteins
highlights = []  #'ORC6','DUSP19','BUB1B','DPH2', 'FLI1']
highlights_ensg = []  #'ORC6','DUSP19','BUB1B','DPH2', 'FLI1']

nbins = 20
protein_heatmap_results = TemporalDelay.protein_heatmap(
    nbins,
    highlights,
    highlights_ensg,
    ccd_comp,
    u_well_plates,
    wp_ensg,
    pol_sort_norm_rev,
    pol_sort_well_plate,
    pol_sort_ab_cell,
    pol_sort_ab_nuc,
    pol_sort_ab_cyto,
    pol_sort_mt_cell,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
)
sorted_maxpol_array, wp_binned_values, wp_max_pol, wp_max_pol_ccd, xvals = (
    protein_heatmap_results
)

# Correlations of known and novel proteins that peak at similar times
TemporalDelay.peak_expression_correlation_analysis(
    wp_binned_values, wp_max_pol, wp_ensg, pol_sort_well_plate, u_well_plates
)

#%% Create a heatmap of peak RNA expression
highlight_names, highlight_ensg = [], []
u_rna_plates = ["355","356","357"]

# Read in RNA-Seq data; use TPMs so that the gene-specific results scales match for cross-gene comparisons
valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_rna_plates
)
adata, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=True
)
adata = RNADataPreparation.zero_center_fucci(adata)
sorted_max_moving_avg_pol_ccd, norm_exp_sort, max_moving_avg_pol, sorted_rna_binned_norm = TemporalDelay.rna_heatmap(
    adata, highlight_names, highlight_ensg, ccdtranscript, xvals
)

# Analyze isoforms
adata_isoform, phases_isoform = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_rna_plates, use_isoforms=True,
)
adata_isoform, phasesfilt_isoform = RNADataPreparation.qc_filtering(
    adata_isoform, do_log_normalize=True, do_remove_blob=True
)
adata_isoform = RNADataPreparation.zero_center_fucci(adata_isoform)
sorted_max_moving_avg_pol_ccd_isoform, norm_exp_sort_isoform, max_moving_avg_pol_isoform, sorted_rna_binned_norm_isoform = TemporalDelay.rna_heatmap(
    adata_isoform,
    highlight_names,
    highlight_ensg,
    ccdtranscript_isoform,
    xvals,
    isIsoformData=True,
)
pearsonCorrelations = TemporalDelay.analyze_ccd_isoform_correlations(
    adata, adata_isoform, ccdtranscript, ccdtranscript_isoform, xvals, u_rna_plates
)

#%% Compare the variances and time of peak expression between protein and RNA
TemporalDelay.compare_variances_prot_v_rna(
    adata,
    norm_exp_sort,
    wp_ensg,
    var_comp_prot,
    gini_comp_prot,
    cv_comp_prot,
    var_cell_prot,
    gini_cell_prot,
    cv_cell_prot,
)

TemporalDelay.compare_peak_expression_prot_v_rna(
    adata,
    wp_ensg,
    ccd_comp,
    ccdtranscript,
    wp_max_pol,
    wp_max_pol_ccd,
    sorted_maxpol_array,
    max_moving_avg_pol,
    sorted_max_moving_avg_pol_ccd,
)
