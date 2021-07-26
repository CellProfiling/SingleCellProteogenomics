# -*- coding: utf-8 -*-
"""
Analysis of protein abundance in individual cells over cell division time.
-  Cell division time is measured with FUCCI markers and modeled in log-log space using polar coordinates.
-  The cell division time modeling is referred to as the "pseudotime" analysis in the paper.
-  Protein abundance is measured in individual asynchronous cells using immunofluorescence and antibody staining.

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import (FucciPseudotime, Loaders,
                                      ProteinBimodality,
                                      ProteinCellCycleDependence,
                                      ProteinVariability, utils)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy

ALL_PLOTS_AND_ANALYSES = False

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

#%% Read in the protein data
import_dict = Loaders.load_protein_fucci_pseudotime()
u_plate, well_plate, well_plate_imgnb, well_plate_imgnb_objnb, u_well_plates = (
    import_dict["u_plate"],
    import_dict["well_plate"],
    import_dict["well_plate_imgnb"],
    import_dict["well_plate_imgnb_objnb"],
    import_dict["u_well_plates"],
)
ab_nuc, ab_cyto, ab_cell, mt_cell = (
    import_dict["ab_nuc"],
    import_dict["ab_cyto"],
    import_dict["ab_cell"],
    import_dict["mt_cell"],
)
area_cell, area_nuc = import_dict["area_cell"], import_dict["area_nuc"]
wp_ensg, wp_ab = import_dict["wp_ensg"], import_dict["wp_ab"]
green_fucci, red_fucci = import_dict["green_fucci"], import_dict["red_fucci"]
log_green_fucci_zeroc, log_red_fucci_zeroc = (
    import_dict["log_green_fucci_zeroc"],
    import_dict["log_red_fucci_zeroc"],
)
log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale = (
    import_dict["log_green_fucci_zeroc_rescale"],
    import_dict["log_red_fucci_zeroc_rescale"],
)
wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp = (
    import_dict["wp_comp_kruskal_gaussccd_adj"],
    import_dict["wp_pass_kruskal_gaussccd_bh_comp"],
)
fucci_data = import_dict["fucci_data"]
wp_iscell, wp_isnuc, wp_iscyto = (
    import_dict["wp_iscell"],
    import_dict["wp_isnuc"],
    import_dict["wp_iscyto"],
)
curr_wp_phases, mockbulk_phases = (
    import_dict["curr_wp_phases"],
    import_dict["mockbulk_phases"],
)

#%%
# Idea: Calculate the polar coordinates and other stuff
# Exec: Devin's calculations
# Output: fucci plot with polar coordinates
pseudotime_result = FucciPseudotime.pseudotime_protein(
    fucci_data,
    ab_nuc,
    ab_cyto,
    ab_cell,
    mt_cell,
    area_cell,
    area_nuc,
    well_plate,
    well_plate_imgnb,
    well_plate_imgnb_objnb,
    log_red_fucci_zeroc_rescale,
    log_green_fucci_zeroc_rescale,
    mockbulk_phases,
)
pol_sort_well_plate, pol_sort_norm_rev, pol_sort_well_plate_imgnb, pol_sort_well_plate_imgnb_objnb, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell, pol_sort_area_cell, pol_sort_area_nuc, pol_sort_fred, pol_sort_fgreen, pol_sort_mockbulk_phases = (
    pseudotime_result
)

#%% Calculate measures of variance of protein abundance in single cells
# Idea: Calculate measures of variance, and show them in plots
# Execution: Now that we already have the data filtered for variability, this is just descriptive.
# Output: scatters of antibody vs microtubule variances by different measures of variaibility

# toggle for using log-transformed intensities
# we decided to use natural intensities
use_log = False
calculate_variaton_result = ProteinVariability.calculate_variation(
    use_log,
    u_well_plates,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
    pol_sort_well_plate,
    pol_sort_ab_cell,
    pol_sort_ab_nuc,
    pol_sort_ab_cyto,
    pol_sort_mt_cell,
    pol_sort_well_plate_imgnb,
)
mean_mean_comp, var_comp, gini_comp, cv_comp, var_cell, gini_cell, cv_cell, var_mt, gini_mt, cv_mt = (
    calculate_variaton_result
)

# Compare variances for protein and microtubules, the internal control for each image
if ALL_PLOTS_AND_ANALYSES:
    removeThese = pd.read_csv("input/ProteinData/ReplicatesToRemove.txt", header=None)[0] # make these independent samples for one-sided Kruskal-Wallis tests
    utils.general_boxplot(
        (var_comp[~np.isin(u_well_plates, removeThese)], var_mt[~np.isin(u_well_plates, removeThese)]),
        ("Protein", "Microtubules"),
        "",
        "Variance",
        "",
        False,
        f"figures/ProteinMicrotubuleVariances.pdf",
    )
    utils.general_boxplot(
        (cv_comp[~np.isin(u_well_plates, removeThese)], gini_mt[~np.isin(u_well_plates, removeThese)]),
        ("Protein", "Microtubules"),
        "",
        "CV",
        "",
        False,
        f"figures/ProteinMicrotubuleCVs.pdf",
    )
    utils.general_boxplot(
        (gini_comp[~np.isin(u_well_plates, removeThese)], gini_mt[~np.isin(u_well_plates, removeThese)]),
        ("Protein", "Microtubules"),
        "",
        "Gini",
        "",
        False,
        f"figures/ProteinMicrotubuleGinis.pdf",
    )
    p_varProt_varMt = 2*scipy.stats.kruskal(var_comp[~np.isin(u_well_plates, removeThese)], var_mt[~np.isin(u_well_plates, removeThese)])[1]
    p_cvProt_cvMt = 2*scipy.stats.kruskal(cv_comp[~np.isin(u_well_plates, removeThese)], cv_mt[~np.isin(u_well_plates, removeThese)])[1]
    p_giniProt_giniMt = 2*scipy.stats.kruskal(gini_comp[~np.isin(u_well_plates, removeThese)], gini_mt[~np.isin(u_well_plates, removeThese)])[1]
    print(
        f"{p_varProt_varMt}: p-value for difference between protein and microtubule variances"
    )
    print(
        f"{p_cvProt_cvMt}: p-value for difference between protein and microtubule CVs"
    )
    print(
        f"{p_giniProt_giniMt}: p-value for difference between protein and microtubule Gini indices"
    )

#%% Gaussian clustering to identify biomodal intensity distributions
bimodal_results = ProteinBimodality.identify_bimodal_intensity_distributions(
    u_well_plates,
    wp_ensg,
    pol_sort_well_plate,
    pol_sort_norm_rev,
    pol_sort_ab_cell,
    pol_sort_ab_nuc,
    pol_sort_ab_cyto,
    pol_sort_mt_cell,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
    ALL_PLOTS_AND_ANALYSES
)
wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs, wp_isbimodal_generally, wp_bimodal_fcmaxmin = (
    bimodal_results
)

#%% Determine cell cycle dependence for each protein
use_log_ccd = False
do_remove_outliers = True
alphaa = 0.05

# Determine cell cycle dependence for proteins
ccd_results = ProteinCellCycleDependence.cell_cycle_dependence_protein(
    u_well_plates,
    wp_ensg,
    wp_ab,
    use_log_ccd,
    do_remove_outliers,
    pol_sort_well_plate,
    pol_sort_norm_rev,
    pol_sort_ab_cell,
    pol_sort_ab_nuc,
    pol_sort_ab_cyto,
    pol_sort_mt_cell,
    pol_sort_fred,
    pol_sort_fgreen,
    pol_sort_mockbulk_phases,
    pol_sort_area_cell,
    pol_sort_area_nuc,
    pol_sort_well_plate_imgnb,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
    wp_isbimodal_fcpadj_pass,
    wp_bimodal_cluster_idxs,
    wp_comp_kruskal_gaussccd_adj,
    ALL_PLOTS_AND_ANALYSES
)
wp_comp_ccd_difffromrng, mean_diff_from_rng_mt, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_ccd_unibimodal, wp_comp_ccd_gauss, perc_var_comp, mean_diff_from_rng, wp_comp_eq_percvar_adj, mean_diff_from_rng_clust1, wp_comp_eq_percvar_adj_clust1, mean_diff_from_rng_clust2, wp_comp_eq_percvar_adj_clust2, mvavgs_x, mvavgs_comp, curr_pols, curr_ab_norms, mvperc_comps, curr_freds, curr_fgreens, curr_mockbulk_phases, curr_area_cell, curr_ab_nuc, curr_well_plate_imgnb, folder = (
    ccd_results
)

# Move the temporal average plots to more informative places
if ALL_PLOTS_AND_ANALYSES:
    ProteinCellCycleDependence.copy_mvavg_plots_protein(
        folder,
        wp_ensg,
        wp_comp_ccd_difffromrng,
        wp_isbimodal_fcpadj_pass,
        wp_comp_ccd_clust1,
        wp_comp_ccd_clust2,
        wp_ccd_unibimodal,
        wp_comp_ccd_gauss,
    )
    ProteinCellCycleDependence.global_plots_protein(
        alphaa,
        u_well_plates,
        wp_ccd_unibimodal,
        perc_var_comp,
        mean_mean_comp,
        gini_comp,
        cv_comp,
        mean_diff_from_rng,
        wp_comp_eq_percvar_adj,
        wp_comp_kruskal_gaussccd_adj,
    )

# Analyze the CCD results and save them
ccd_comp, nonccd_comp, bioccd = ProteinCellCycleDependence.analyze_ccd_variation_protein(
    folder,
    u_well_plates,
    wp_ensg,
    wp_ab,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
    wp_comp_ccd_difffromrng,
    wp_comp_ccd_clust1,
    wp_comp_ccd_clust2,
    var_comp,
    gini_comp,
    perc_var_comp,
    mean_diff_from_rng,
    wp_comp_kruskal_gaussccd_adj,
    wp_comp_eq_percvar_adj,
    mean_diff_from_rng_clust1,
    wp_comp_eq_percvar_adj_clust1,
    mean_diff_from_rng_clust2,
    wp_comp_eq_percvar_adj_clust2,
    wp_isbimodal_fcpadj_pass,
    wp_isbimodal_generally,
    wp_ccd_unibimodal,
    wp_bimodal_fcmaxmin,
    wp_comp_ccd_gauss,
)

# Make a dataframe for plotting on the HPA website
ProteinCellCycleDependence.make_plotting_dataframe(
    wp_ensg,
    wp_ab,
    u_well_plates,
    wp_iscell,
    wp_iscyto,
    wp_isnuc,
    ccd_comp,
    bioccd,
    curr_pols,
    curr_ab_norms,
    curr_freds,
    curr_fgreens,
    curr_mockbulk_phases,
    mvavgs_x,
    mvavgs_comp,
    mvperc_comps,
    gini_comp,
    perc_var_comp,
)

# Perform comparison to LASSO for finding CCD proteins
if ALL_PLOTS_AND_ANALYSES:
    ProteinCellCycleDependence.compare_to_lasso_analysis(
        u_well_plates,
        pol_sort_norm_rev,
        pol_sort_well_plate,
        pol_sort_ab_cell,
        pol_sort_ab_nuc,
        pol_sort_ab_cyto,
        pol_sort_mt_cell,
        pol_sort_fred,
        pol_sort_fgreen,
        wp_iscell,
        wp_isnuc,
        wp_iscyto,
        wp_ensg,
        ccd_comp,
    )
    
    # Generate UMAPs to illustrate cutoffs and stability
    ProteinCellCycleDependence.generate_protein_umaps(
        u_well_plates,
        pol_sort_norm_rev,
        pol_sort_well_plate,
        pol_sort_ab_cell,
        pol_sort_ab_nuc,
        pol_sort_ab_cyto,
        pol_sort_mt_cell,
        wp_iscell,
        wp_isnuc,
        wp_iscyto,
        mean_diff_from_rng,
    )
