# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:51:50 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, MovingAverages

WINDOW = 10 # Number of points for moving average window
PERMUTATIONS = 10000

def clust_to_wp(clust, clust_idx):
    wp_clust = np.array([False] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust

def clust_to_wp_doub(clust, clust_idx):
    wp_clust = np.array([0.0] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust

def remove_outliers_idx(values):
    '''Returns indices of outliers to remove'''
    max_cutoff = np.mean(values) + 5 * np.std(values)
    min_cutoff = np.mean(values) - 5 * np.std(values)
    return (values < max_cutoff) & (values > min_cutoff)

def remove_outliers(values, return_values):
    '''Remove outliers on "values" and return "return_values" based on that filter'''
    return return_values[remove_outliers_idx(values)]

def permutation_analysis_protein(idx, curr_pol, curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm,
                                 perc_var_cell_val, perc_var_nuc_val, perc_var_cyto_val,
                                 wp_iscell, wp_isnuc, wp_iscyto,
                                 mvavg_cell, mvavg_nuc, mvavg_cyto):
    '''Randomization analysis of cell cycle dependence: permute cell order and calculate percent variance due to the cell cycle'''
    perms = np.asarray([np.random.permutation(len(curr_pol)) for nnn in np.arange(PERMUTATIONS)])
    curr_comp_norm = np.asarray(curr_ab_cell_norm if wp_iscell[idx] else curr_ab_nuc_norm if wp_isnuc[idx] else curr_ab_cyto_norm)
    curr_comp_percvar = np.asarray(perc_var_cell_val if wp_iscell[idx] else perc_var_nuc_val if wp_isnuc[idx] else perc_var_cyto_val)
    curr_comp_mvavg = np.asarray(mvavg_cell if wp_iscell[idx] else mvavg_nuc if wp_isnuc[idx] else mvavg_cyto)
    curr_comp_perm = np.asarray([curr_comp_norm[perm] for perm in perms])
    curr_mt_perm = np.asarray([curr_mt_cell_norm[perm] for perm in perms])
    curr_mvavg_rng_comp = np.apply_along_axis(MovingAverages.mvavg, 1, curr_comp_perm, WINDOW)
    curr_mvavg_rng_mt = np.apply_along_axis(MovingAverages.mvavg, 1, curr_mt_perm, WINDOW)
    curr_percvar_rng_comp = np.var(curr_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm, axis=1)
    curr_percvar_rng_mt = np.var(curr_mvavg_rng_mt, axis=1) / np.var(curr_mt_perm, axis=1)
    return (curr_comp_norm, curr_comp_percvar, curr_comp_mvavg, curr_comp_perm, curr_mt_perm, 
        curr_mvavg_rng_comp, curr_mvavg_rng_mt, curr_percvar_rng_comp, curr_percvar_rng_mt)
    
def cell_cycle_dependence_protein(u_well_plates, wp_ensg, use_log_ccd, do_remove_outliers,
                                  pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
                                  pol_sort_area_cell, pol_sort_area_nuc,
                                  wp_iscell, wp_isnuc, wp_iscyto,
                                  wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs):
    '''
    Use a moving average model of protein expression over the cell cycle to determine cell cycle dependence.
    Generates plots for each gene
    '''
    xvals = np.linspace(0, 1, num=21)
    perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = [],[],[],[] # percent variance attributed to cell cycle (mean POI intensities)
    perc_var_comp_rng, perc_var_mt_rng = [],[] # randomized in pseudotime; percent variances
    perc_var_comp_clust1, perc_var_comp_clust2, mvavgs_comp_clust1, mvavgs_comp_clust2, perc_var_comp_clust1_rng, perc_var_comp_clust2_rng, mvavgs_x_clust1, mvavgs_x_clust2 = [],[],[],[],[],[],[],[] # percent variances for bimodal
    mvavgs_cell, mvavgs_nuc, mvavgs_cyto, mvavgs_mt, mvavgs_x = [],[],[],[],[] # moving average y values & x value
    mvavgs_cell_rng, mvavgs_nuc_rng, mvavgs_cyto_rng, mvavgs_mt_rng = [],[],[],[] # moving average y values from randomization
    cell_counts = []
    
    analysis = "MeanRng"
    folder = f"figures/TemporalMovingAverages{analysis}191205"
    folder_mt = f"figures/TemporalMovingAverages{analysis}191205_mt"
    folder_rng = f"figures/TemporalMovingAverageRandomizationExamples"
    if not os.path.exists(folder): os.mkdir(folder)
    if not os.path.exists(folder_mt): os.mkdir(folder_mt)
    if not os.path.exists(folder_rng): os.mkdir(folder_rng)
    fileprefixes = np.array([f"{ensg}_{sum(wp_ensg[:ei] == ensg)}" for ei, ensg in enumerate(wp_ensg)])
    
    for i, well in enumerate(u_well_plates):
    #    print(well)
        plt.close('all')
        if i % 100 == 0: print(f"well {i} of {len(u_well_plates)}")
    #    well = 'H05_55405991'#GMNN well, used for testing
        curr_well_inds = pol_sort_well_plate==well # the reversal isn't really helpful here
        curr_pol = pol_sort_norm_rev[curr_well_inds]
        curr_ab_cell = pol_sort_ab_cell[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_cell[curr_well_inds])
        curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_nuc[curr_well_inds])
        curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_cyto[curr_well_inds])
        curr_mt_cell = pol_sort_mt_cell[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_mt_cell[curr_well_inds])
        if do_remove_outliers:
            curr_comp = curr_ab_cell if wp_iscell[i] else curr_ab_nuc if wp_isnuc[i] else curr_ab_cyto
            curr_pol = remove_outliers(curr_comp, curr_pol)
            curr_ab_cell = remove_outliers(curr_comp,curr_ab_cell)
            curr_ab_nuc = remove_outliers(curr_comp,curr_ab_nuc)
            curr_ab_cyto = remove_outliers(curr_comp,curr_ab_cyto)
            curr_mt_cell = remove_outliers(curr_comp,curr_mt_cell)
    
        # Normalize mean intensities, normalized for display
        curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
        curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
        curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
        curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
        curr_area_cell_norm = pol_sort_area_cell[curr_well_inds] / np.max(pol_sort_area_cell[curr_well_inds])
        curr_area_nuc_norm = pol_sort_area_nuc[curr_well_inds] / np.max(pol_sort_area_nuc[curr_well_inds])
            
        # Original method from Devin's work
        perc_var_cell_val, mvavg_cell = MovingAverages.mvavg_perc_var(curr_ab_cell_norm, WINDOW)
        perc_var_nuc_val, mvavg_nuc = MovingAverages.mvavg_perc_var(curr_ab_nuc_norm, WINDOW)
        perc_var_cyto_val, mvavg_cyto = MovingAverages.mvavg_perc_var(curr_ab_cyto_norm, WINDOW)
        perc_var_mt_val, mvavg_mt = MovingAverages.mvavg_perc_var(curr_mt_cell_norm, WINDOW)
        mvavg_xvals = MovingAverages.mvavg(curr_pol, WINDOW)
        
        # Permutation analysis
        permutation_result = permutation_analysis_protein(i, curr_pol, curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm,
                                 perc_var_cell_val, perc_var_nuc_val, perc_var_cyto_val, wp_iscell, wp_isnuc, wp_iscyto, mvavg_cell, mvavg_nuc, mvavg_cyto)
        curr_comp_norm, curr_comp_percvar, curr_comp_mvavg, curr_comp_perm, curr_mt_perm, curr_mvavg_rng_comp, curr_mvavg_rng_mt, curr_percvar_rng_comp, curr_percvar_rng_mt = permutation_result
        perc_var_comp_rng.append(curr_percvar_rng_comp)
        perc_var_mt_rng.append(curr_percvar_rng_mt)
        
        # Assess CCD in high- and low-expressing cell populations of proteins determined to have bimodal population intensity
        clusters = None
        if wp_isbimodal_fcpadj_pass[i]:
            clust1_idx, clust2_idx = wp_bimodal_cluster_idxs[i]
            if do_remove_outliers:
                clust1_idx, clust2_idx = remove_outliers(curr_comp, clust1_idx), remove_outliers(curr_comp, clust2_idx)
            perc_var_comp_clust1_val, mvavg_clust1 = MovingAverages.mvavg_perc_var(curr_comp_norm[clust1_idx], WINDOW)
            perc_var_comp_clust2_val, mvavg_clust2 = MovingAverages.mvavg_perc_var(curr_comp_norm[clust2_idx], WINDOW)
            mvavgs_x_clust1.append(MovingAverages.mvavg(curr_pol[clust1_idx], WINDOW))
            mvavgs_x_clust2.append(MovingAverages.mvavg(curr_pol[clust2_idx], WINDOW))
            perc_var_comp_clust1.append(perc_var_comp_clust1_val)
            perc_var_comp_clust2.append(perc_var_comp_clust2_val)
            mvavgs_comp_clust1.append(mvavg_clust1)
            mvavgs_comp_clust2.append(mvavg_clust2)
            
            clust1gt = np.mean(mvavg_clust1) > np.mean(mvavg_clust2)
            clusters = np.array([clust1gt] * len(clust1_idx))
            clusters[clust2_idx] = not clust1gt
            
            perms1 = np.asarray([np.random.permutation(sum(clust1_idx)) for nnn in np.arange(PERMUTATIONS)])
            perms2 = np.asarray([np.random.permutation(sum(clust2_idx)) for nnn in np.arange(PERMUTATIONS)])
            curr_comp_perm1 = np.asarray([curr_comp_norm[clust1_idx][perm] for perm in perms1])
            curr_comp_perm2 = np.asarray([curr_comp_norm[clust2_idx][perm] for perm in perms2])
            curr_clust1_mvavg_rng_comp = np.apply_along_axis(MovingAverages.mvavg, 1, curr_comp_perm1, WINDOW)
            curr_clust2_mvavg_rng_comp = np.apply_along_axis(MovingAverages.mvavg, 1, curr_comp_perm2, WINDOW)
            curr_clust1_percvar_rng_comp = np.var(curr_clust1_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm1, axis=1)
            curr_clust2_percvar_rng_comp = np.var(curr_clust2_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm2, axis=1)
            perc_var_comp_clust1_rng.append(curr_clust1_percvar_rng_comp)
            perc_var_comp_clust2_rng.append(curr_clust2_percvar_rng_comp)
            
            windows1 = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(sum(clust1_idx) - WINDOW + 1)])
            windows2 = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(sum(clust2_idx) - WINDOW + 1)])
            MovingAverages.temporal_mov_avg_protein(curr_pol[clust1_idx], curr_comp_norm[clust1_idx], mvavgs_x_clust1[-1], mvavg_clust1, windows1, None, folder, fileprefixes[i] + "_clust1")
            MovingAverages.temporal_mov_avg_protein(curr_pol[clust2_idx], curr_comp_norm[clust2_idx], mvavgs_x_clust2[-1], mvavg_clust2, windows2, None, folder, fileprefixes[i] + "_clust2")
        
        # Make example plots for the randomization trials, but no need for the bimodal ones
        elif np.mean(curr_comp_percvar - curr_percvar_rng_comp) > 0.05:
            median_rng_idx = np.argsort(curr_percvar_rng_comp)
            for iii, idx in enumerate([0,len(curr_percvar_rng_comp)//4,len(curr_percvar_rng_comp)//2,3*len(curr_percvar_rng_comp)//4,len(curr_percvar_rng_comp)-1]):
                MovingAverages.temporal_mov_avg_randomization_example_protein(curr_pol, curr_comp_norm, curr_comp_perm[median_rng_idx[idx]], 
                                         mvavg_xvals, curr_comp_mvavg, curr_mvavg_rng_comp[median_rng_idx[idx]], 
                                         folder_rng, f"{fileprefixes[i]}_withrandomization_{iii}")
    
        # Test for equal variances of the moving averages and raw values
        perc_var_cell.append(perc_var_cell_val)
        perc_var_nuc.append(perc_var_nuc_val)
        perc_var_cyto.append(perc_var_cyto_val)
        perc_var_mt.append(perc_var_mt_val)
        
        mvavgs_cell.append(mvavg_cell)
        mvavgs_nuc.append(mvavg_nuc)
        mvavgs_cyto.append(mvavg_cyto)
        mvavgs_mt.append(mvavg_mt)
        mvavgs_x.append(mvavg_xvals)
        
        cell_counts.append(len(curr_pol))
        percvar = perc_var_cell_val if wp_iscell[i] else perc_var_nuc_val if wp_isnuc[i] else perc_var_cyto_val
        
        # Uncomment to make the plots (takes 10 mins)
        windows = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(len(curr_pol) - WINDOW + 1)])
        MovingAverages.temporal_mov_avg_protein(curr_pol, curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm, mvavg_xvals,
             mvavg_cell if wp_iscell[i] else mvavg_nuc if wp_isnuc[i] else mvavg_cyto,
             windows, None, folder, fileprefixes[i])
        
        # Uncomment to make the plots for microtubules (takes 10 mins)
        MovingAverages.temporal_mov_avg_protein(curr_pol, curr_mt_cell_norm, mvavg_xvals, mvavg_mt, windows, None, folder_mt, f"{fileprefixes[i]}_mt")
        
        if clusters is not None:
            MovingAverages.temporal_mov_avg_protein(curr_pol, curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm, mvavg_xvals,
                 mvavg_cell if wp_iscell[i] else mvavg_nuc if wp_isnuc[i] else mvavg_cyto,
                 windows, clusters, folder, fileprefixes[i] + "_clust1&2")
        
    alpha_ccd = 0.01
    perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = np.array(perc_var_cell),np.array(perc_var_nuc),np.array(perc_var_cyto),np.array(perc_var_mt) # percent variance attributed to cell cycle (mean POI intensities)
    perc_var_mt_rng, perc_var_comp_rng = np.array(perc_var_mt_rng), np.array(perc_var_comp_rng) 
    
    # Let's check out which percent variances are greater than the permuted values
    perc_var_comp = values_comp(perc_var_cell, perc_var_nuc, perc_var_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    perc_var_comp_withbimodal = np.concatenate((perc_var_comp, perc_var_comp_clust1, perc_var_comp_clust2))
    perc_var_comp_rng_withbimodal = np.concatenate((perc_var_comp_rng, perc_var_comp_clust1_rng, perc_var_comp_clust2_rng))
    ccd_var_comp_rng_wilcoxp_withbimodal = np.apply_along_axis(scipy.stats.wilcoxon, 1, (perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T, None, "wilcox", False, "greater").T[1].T
    ccd_var_mt_rng_wilcoxp = np.apply_along_axis(scipy.stats.wilcoxon, 1, (perc_var_mt - perc_var_mt_rng.T).T, None, "wilcox", False, "greater").T[1].T
    
    # randomization tests, try being a bit more stringent, try drawing the cutoff based on microtubules per sample
    wp_comp_eq_percvar_adj_withbimodal, wp_comp_pass_eq_percvar_adj_withbimodal = bonf(alpha_ccd, ccd_var_comp_rng_wilcoxp_withbimodal)
    wp_comp_gtpass_eq_percvar_adj_withbimodal = wp_comp_pass_eq_percvar_adj_withbimodal & (perc_var_comp_withbimodal > np.median(perc_var_comp_rng_withbimodal, axis=1))
    
    # median differences from random
    MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM = 0.08
    print(f"Requiring {MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM*100}% additional percent variance explained than random.")
    mean_diff_from_rng_withbimodal = np.mean((perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T, 1)
    wp_comp_ccd_difffromrng_withbimodal = mean_diff_from_rng_withbimodal >= MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM
    
    ###### Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
    # separate unimodal from bimodal again
    wp_comp_eq_percvar_adj = wp_comp_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
    wp_comp_pass_eq_percvar_adj = wp_comp_pass_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
    wp_comp_gtpass_eq_percvar_adj = wp_comp_gtpass_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
    mean_diff_from_rng = mean_diff_from_rng_withbimodal[:len(perc_var_comp)]
    wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng_withbimodal[:len(perc_var_comp)]

    wp_comp_pass_eq_percvar_adj_clust1 = clust_to_wp(wp_comp_pass_eq_percvar_adj_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
    wp_comp_eq_percvar_adj_clust1 = clust_to_wp(wp_comp_eq_percvar_adj_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
    mean_diff_from_rng_clust1 = clust_to_wp_doub(mean_diff_from_rng_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
    wp_comp_ccd_difffromrng_clust1 = clust_to_wp(wp_comp_ccd_difffromrng_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
    wp_comp_pass_eq_percvar_adj_clust2 = clust_to_wp(wp_comp_pass_eq_percvar_adj_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
    wp_comp_eq_percvar_adj_clust2 = clust_to_wp(wp_comp_eq_percvar_adj_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
    mean_diff_from_rng_clust2 = clust_to_wp_doub(mean_diff_from_rng_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
    wp_comp_ccd_difffromrng_clust2 = clust_to_wp(wp_comp_ccd_difffromrng_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
    
    wp_normal_randompercvar_p = np.apply_along_axis(scipy.stats.normaltest, 1, (perc_var_comp - perc_var_comp_rng.T).T).T[1].T
    wp_randompercvarnorm_adj, wp_randompercvarnorm_pass = benji_hoch(0.05, wp_normal_randompercvar_p)
    print(f"{sum(wp_randompercvarnorm_pass)}: number of genes with randomized percvars that form normal distributions")
    
    wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng
    print(f"{sum(wp_comp_ccd_difffromrng)}: # proteins showing CCD variation unimodally, comp, percvar rng median diff")
    print(f"{sum(wp_comp_ccd_difffromrng) / len(wp_comp_ccd_difffromrng)}: fraction of variable proteins showing CCD variation, comp, percvar rng median diff")
    wp_comp_ccd_gauss = wp_comp_kruskal_gaussccd_adj <= alpha_ccd
    print(f"{sum(wp_comp_ccd_gauss)}: # proteins showing CCD variation, comp, gaussian analysis")
    print(f"{sum(wp_comp_ccd_gauss) / len(wp_comp_ccd_gauss)}: fraction of variable proteins showing CCD variation, comp, gaussian analysis")
    print(f"{sum(wp_comp_ccd_gauss & wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, CCD and gaussian")
    print(f"{sum(wp_comp_ccd_gauss & ~wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, not CCD and gaussian")
    print(f"{sum(~wp_comp_ccd_gauss & wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, CCD and not gaussian")
    
    # Address bimodal ones
    # 1) The number that are bimodal in one cluster (number of those that are also CCD as unimodal)
    # 2) The number that are bimodal in both clusters (number of those that are also CCD as unimodal)
    wp_comp_ccd_clust1 = wp_comp_ccd_difffromrng_clust1
    wp_comp_ccd_clust2 = wp_comp_ccd_difffromrng_clust2
    print(f"{sum(wp_isbimodal_fcpadj_pass)}: samples with bimodal antibody intensities")
    print(f"{sum(wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2)}: bimodal samples with one CCD cluster ({sum((wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2) & wp_comp_ccd_difffromrng)}: also CCD unimodally)")
    print(f"{sum(wp_comp_ccd_clust1 & wp_comp_ccd_clust2)}: bimodal samples with two CCD clusters ({sum((wp_comp_ccd_clust1 & wp_comp_ccd_clust2) & wp_comp_ccd_difffromrng)}: also CCD unimodally)")
    
    wp_comp_ccd_use = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2 # percvar randomization, not like in original manuscript
    
    plt.figure(figsize=(10,10))
    plt.scatter(perc_var_comp_withbimodal, mean_diff_from_rng_withbimodal, c=wp_comp_ccd_difffromrng_withbimodal, cmap="bwr_r")
    plt.hlines(MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM, np.min(perc_var_comp), np.max(perc_var_comp), color="gray")
    plt.xlabel("Percent Variance Explained by Cell Cycle")
    plt.ylabel("Mean Difference from Random")
    plt.savefig("figures/MedianDiffFromRandom.png")
    plt.savefig("figures/MedianDiffFromRandom.pdf")
    plt.show()
    plt.close()
    
    pervar_adj_withbimodal_nextafter = np.nextafter(wp_comp_eq_percvar_adj_withbimodal, wp_comp_eq_percvar_adj_withbimodal + 1)
    plt.figure(figsize=(10,10))
    plt.scatter(mean_diff_from_rng_withbimodal, -np.log10(pervar_adj_withbimodal_nextafter), c=wp_comp_ccd_difffromrng_withbimodal, cmap="bwr_r")
    plt.vlines(MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM, np.min(-np.log10(pervar_adj_withbimodal_nextafter)), np.max(-np.log10(wp_comp_eq_percvar_adj_withbimodal)), color="gray")
    plt.xlabel("Mean Difference from Random")
    plt.ylabel("-log10 adj p-value from randomization")
    plt.savefig("figures/MedianDiffFromRandomVolcano.png")
    plt.savefig("figures/MedianDiffFromRandomVolcano.pdf")
    plt.show()
    plt.close()
    
    return wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_comp_ccd_gauss, folder

def copy_mvavg_plots_protein(folder, wp_comp_ccd_difffromrng, wp_isbimodal_fcpadj_pass, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_comp_ccd_gauss):
    # Copy profiles to the right place:
    # 1) CCD Unimodal
    # 2) CCD Bimodal
    # 3) Non-CCD
    # 3) Gaussian analysis and CCD (unimodal or bimodal)
    # 4) Gaussian analysis and non-CCD (unimodal or bimodal)
    ccdunifolder = f"figures/CCDUnimodal"
    ccdunibifolder = f"figures/CCDUnimodalAndBimodal"
    ccdpbifolder = f"figures/CCDBimodal"
    ccdgaussccdfolder = f"figures/GaussAndCCD"
    ccdgaussnonccdfolder = f"figures/GaussAndNonCCD"
    nongaussccdfolder = f"figures/CCDAndNonGauss"
    nonccdfolder = f"figures/NonCCD"
    bimodalnonccdfolder = f"figures/NonCCDBimodal"
    replicatefolder = f"figures/Replicates"
    examplesfolder = f"figures/Examples"
    for f in [ccdunifolder,ccdunibifolder,ccdpbifolder,ccdgaussccdfolder,ccdgaussnonccdfolder,nongaussccdfolder,nonccdfolder,bimodalnonccdfolder,replicatefolder,examplesfolder]:
        if not os.path.exists(f): os.mkdir(f)
    
    # CCD Unimodal
    for ensg in fileprefixes[wp_comp_ccd_difffromrng]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdunifolder, ensg +'_mvavg.pdf'))
        
    # CCD Unimodal and Bimodal
    for ensg in fileprefixes[wp_comp_ccd_difffromrng & (wp_comp_ccd_clust1 | wp_comp_ccd_clust2)]:
        shutil.copy(os.path.join(folder, ensg+'_clust1&2_mvavg.pdf'), os.path.join(ccdunibifolder, ensg +'_clust1&2_mvavg.pdf'))
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdunibifolder, ensg +'_mvavg.pdf'))
    for ensg in fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
        shutil.copy(os.path.join(folder, ensg+'_clust1_mvavg.pdf'), os.path.join(ccdunibifolder, ensg+'_clust1_mvavg.pdf'))
    for ensg in fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
        shutil.copy(os.path.join(folder, ensg+'_clust2_mvavg.pdf'), os.path.join(ccdunibifolder, ensg+'_clust2_mvavg.pdf'))
        
    # CCD Bimodal
    for ensg in fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
        shutil.copy(os.path.join(folder, ensg+'_clust1&2_mvavg.pdf'), os.path.join(ccdpbifolder, ensg +'_clust1&2_mvavg.pdf'))
        shutil.copy(os.path.join(folder, ensg+'_clust1_mvavg.pdf'), os.path.join(ccdpbifolder, ensg+'_clust1_mvavg.pdf'))
    for ensg in fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
        shutil.copy(os.path.join(folder, ensg+'_clust1&2_mvavg.pdf'), os.path.join(ccdpbifolder, ensg +'_clust1&2_mvavg.pdf'))
        shutil.copy(os.path.join(folder, ensg+'_clust2_mvavg.pdf'), os.path.join(ccdpbifolder, ensg+'_clust2_mvavg.pdf'))
    
    # Non-CCD 
    for ensg in fileprefixes[~wp_comp_ccd_difffromrng & ~wp_comp_ccd_clust1 & ~wp_comp_ccd_clust2]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(nonccdfolder, ensg+'_mvavg.pdf'))
        
    # Non-CCD Bimodal
    for ensg in fileprefixes[wp_isbimodal_fcpadj_pass & ~wp_comp_ccd_clust1 & ~wp_comp_ccd_clust2]:
        shutil.copy(os.path.join(folder, ensg+'_clust1&2_mvavg.pdf'), os.path.join(bimodalnonccdfolder, ensg +'_clust1&2_mvavg.pdf'))
        shutil.copy(os.path.join(folder, ensg+'_clust1_mvavg.pdf'), os.path.join(bimodalnonccdfolder, ensg+'_clust1_mvavg.pdf'))
        shutil.copy(os.path.join(folder, ensg+'_clust2_mvavg.pdf'), os.path.join(bimodalnonccdfolder, ensg+'_clust2_mvavg.pdf'))
        
    # Gauss and CCD
    wp_ccd_unibimodal = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
    for ensg in fileprefixes[wp_comp_ccd_gauss & wp_ccd_unibimodal]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdgaussccdfolder, ensg+'_mvavg.pdf'))
    # Gauss and Non-CCD
    for ensg in fileprefixes[wp_comp_ccd_gauss & ~wp_ccd_unibimodal]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdgaussnonccdfolder, ensg+'_mvavg.pdf'))
    # Non-Gauss and CCD
    for ensg in fileprefixes[~wp_comp_ccd_gauss & wp_ccd_unibimodal]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(nongaussccdfolder, ensg+'_mvavg.pdf'))
        
    # Examples
    for ensg in fileprefixes[np.isin(wp_ensg,["ENSG00000011426","ENSG00000111605","ENSG00000102908","ENSG00000091651","ENSG00000083812","ENSG00000162999","ENSG00000134057","ENSG00000178999","ENSG00000156970","ENSG00000132768","ENSG00000138801","ENSG00000156239","ENSG00000019144","ENSG00000151702","ENSG00000123607","ENSG00000173599","ENSG00000109814"])]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(examplesfolder, ensg+'_mvavg.pdf'))

def global_plots_protein(alphaa, u_well_plates, perc_var_comp, mean_mean_comp, mean_diff_from_rng, wp_comp_eq_percvar_adj):
    '''Illustrate the CCD variances of all proteins'''
    u_well_plates_list = list(u_well_plates)
    utils.general_scatter(perc_var_comp, mean_mean_comp, "percent variance", "mean mean intensity", f"figures/PercVarVsMeanMeanIntensity_comp.png")
    utils.general_scatter(mean_mean_comp, mean_diff_from_rng, "Mean Mean Intensity", "Mean Additional Percent Variance Explained than Random", f"figures/IntensityVsMeanDiff.png")
    utils.general_scatter_color(gini_comp, perc_var_comp, "Gini of Protein Expression", "Fraction of Variance Due to Cell Cycle",
                                -np.log10(wp_comp_kruskal_gaussccd_adj), "FDR for Cell Cycle Dependence", True,
                                "Compartment - Fraction of Variance Due to Cell Cycle", "figures/CompartmentProteinFractionVariance.png")
    utils.general_scatter_color(gini_comp, perc_var_comp, "Gini of Protein Expression", "Fraction of Variance Due to Cell Cycle",
                                wp_comp_ccd_use, "", False, 
                                "Compartment - Fraction of Variance Due to Cell Cycle", "figures/CompartmentProteinFractionVarianceTF.png",
                                "bwr_r", 0.5)
    utils.general_scatter_color(cv_comp, perc_var_comp, "CV of Protein Expression", "Fraction of Variance Due to Cell Cycle",
                                -np.log10(wp_comp_kruskal_gaussccd_adj), "FDR for Cell Cycle Dependence", True, 
                                "Compartment - Fraction of Variance Due to Cell Cycle", "figures/CompartmentCVProteinFractionVariance.png")
    
    plt.figure(figsize=(10,10))
    plt.scatter(perc_var_comp, -np.log10(wp_comp_eq_percvar_adj))
    plt.xlabel("percent variance new")
    plt.ylabel("-log10 FDR for CCD")
    plt.hlines(-np.log10(alphaa), np.min(perc_var_comp), np.max(perc_var_comp))
    plt.savefig(f"figures/PercVarVsLog10FdrCCD_comp.png")
    plt.show()
    plt.close()
    
def analyze_replicates(wp_ccd_with_replicates, wp_ensg, analysis_tag):
    wp_ensg_counts = np.array([sum([eeee == ensg for eeee in wp_ensg]) for ensg in wp_ensg])
    ensg_is_duplicated = wp_ensg_counts > 1
    duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
    duplicated_ensg_ccd = np.array([sum(wp_ccd_with_replicates[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(wp_ccd_with_replicates[~ensg_is_duplicated])}: number of CCD proteins (no replicate, {analysis_tag})")
    print(f"{sum(~wp_ccd_with_replicates[~ensg_is_duplicated])}: number of non-CCD proteins (no replicate, {analysis_tag})")
    print(f"{sum(duplicated_ensg_ccd == 2)}: number of replicated stainings shown to be CCD in both replicates, {analysis_tag}")
    print(f"{sum(duplicated_ensg_ccd == 1)}: number of replicated stainings shown to be CCD in just one replicate, {analysis_tag}")
    print(f"{sum(duplicated_ensg_ccd == 0)}: number of replicated stainings shown to be non-CCD in both replicate,  {analysis_tag}")
    
def analyze_ccd_variation_protein(u_well_plates, wp_ensg, wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2, wp_ccd_unibimodal, wp_isbimodal_fcpadj_pass):
    n_tot_variable = len(u_well_plates)
    ccd_comp = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
    nonccd_comp = ~ccd_comp
    print(f"{n_tot_variable}: # total proteins showing variation")
    print(f"{sum(ccd_comp)}: CCD variable proteins (before addressing redundancy and mitotic structures)")
    print(f"{sum(nonccd_comp)}: non-CCD variable proteins")
    
    ### address gene redundancy
    wp_ccd_unibimodal = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
    wp_ccd_bimodalonecluster = wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2
    wp_ccd_bimodaltwocluster = wp_comp_ccd_clust1 & wp_comp_ccd_clust2
    wp_ensg_counts = np.array([sum([eeee == ensg for eeee in wp_ensg]) for ensg in wp_ensg])
    ensg_is_duplicated = wp_ensg_counts > 1
    duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
    duplicated_ensg_pairs = [u_well_plates[wp_ensg == ensg] for ensg in duplicated_ensg]
    print(f"{sum(wp_ccd_unibimodal[~ensg_is_duplicated])}: number of CCD proteins (no replicate, unimodal and bimodal)")
    print(f"{sum(~wp_ccd_unibimodal[~ensg_is_duplicated])}: number of non-CCD proteins (no replicate, unimodal and bimodal)")
    print(f"{sum(wp_ccd_bimodalonecluster[~ensg_is_duplicated])}: bimodal samples with one CCD cluster ({sum((wp_ccd_bimodalonecluster & wp_comp_ccd_difffromrng)[~ensg_is_duplicated])}: also CCD unimodally), no replicate")
    print(f"{sum(wp_ccd_bimodaltwocluster[~ensg_is_duplicated])}: bimodal samples with two CCD clusters ({sum((wp_ccd_bimodaltwocluster & wp_comp_ccd_difffromrng)[~ensg_is_duplicated])}: also CCD unimodally), no replicate")
    duplicated_ensg_ccd = np.array([sum(wp_ccd_unibimodal[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(duplicated_ensg_ccd == 2)}: number of replicated stainings shown to be CCD in both replicates, unimodal and bimodal")
    print(f"{sum(duplicated_ensg_ccd == 1)}: number of replicated stainings shown to be CCD in just one replicate, unimodal and bimodal")
    print(f"{sum(duplicated_ensg_ccd == 0)}: number of replicated stainings shown to be non-CCD in both replicate, unimodal and bimodal")
    duplicated_ensg_ccd_bi1 = np.array([sum(wp_ccd_bimodalonecluster[wp_ensg == ensg]) for ensg in duplicated_ensg])
    duplicated_ensg_ccd_plusunimodal = np.array([sum((wp_comp_ccd_difffromrng & wp_ccd_bimodalonecluster)[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(duplicated_ensg_ccd_bi1 == 2)}: number of replicated stainings shown to be CCD in both replicates, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 2)} also unimodally)")
    print(f"{sum(duplicated_ensg_ccd_bi1 == 1)}: number of replicated stainings shown to be CCD in just one replicate, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 1)} also unimodally)")
    print(f"{sum(duplicated_ensg_ccd_bi1 == 0)}: number of replicated stainings shown to be non-CCD in both replicate, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 0)} also unimodally)")
    duplicated_ensg_ccd_bi2 = np.array([sum(wp_ccd_bimodaltwocluster[wp_ensg == ensg]) for ensg in duplicated_ensg])
    duplicated_ensg_ccd_plusunimodal = np.array([sum((wp_comp_ccd_difffromrng & wp_ccd_bimodaltwocluster)[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(duplicated_ensg_ccd_bi2 == 2)}: number of replicated stainings shown to be CCD in both replicates, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 2)} also unimodally)")
    print(f"{sum(duplicated_ensg_ccd_bi2 == 1)}: number of replicated stainings shown to be CCD in just one replicate, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 1)} also unimodally)")
    print(f"{sum(duplicated_ensg_ccd_bi2 == 0)}: number of replicated stainings shown to be non-CCD in both replicate, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 0)} also unimodally)")
    
    bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
    protein_ct = len(np.unique(np.concatenate((wp_ensg, bioccd))))
    ccd_protein_ct = sum(wp_ccd_unibimodal[~ensg_is_duplicated]) + sum(duplicated_ensg_ccd == 2)
    nonccd_protein_ct = sum(~wp_ccd_unibimodal[~ensg_is_duplicated]) + sum(duplicated_ensg_ccd == 0)
    duplicated_ensg_bimodal_generally = np.array([sum(wp_isbimodal_generally[wp_ensg == ensg]) for ensg in duplicated_ensg])
    unimodal_generally_protein_ct = sum(~wp_isbimodal_generally[~ensg_is_duplicated]) + sum(duplicated_ensg_bimodal_generally < 2)
    bimodal_generally_protein_ct = sum(wp_isbimodal_generally[~ensg_is_duplicated]) + sum(duplicated_ensg_bimodal_generally == 2)
    print(f"{ccd_protein_ct}: number of ccd proteins; addressed replicates; not including mitotic structures")
    
    # Decision: remove proteins that didn't pass CCD in both replicates
    ccd_comp[np.isin(wp_ensg, duplicated_ensg[duplicated_ensg_ccd == 1])] = False
    nonccd_comp[np.isin(wp_ensg, duplicated_ensg[duplicated_ensg_ccd == 1])] = False
    
    # Analyze replicates for CCD proteins (both mock-bulk and pseudotime analysis)
    analyze_replicates(wp_comp_ccd_gauss, wp_ensg, "gauss")
    analyze_replicates(wp_comp_ccd_gauss & wp_comp_ccd_difffromrng, wp_ensg, "CCD and gaussian")
    analyze_replicates(wp_comp_ccd_gauss & ~wp_comp_ccd_difffromrng, wp_ensg, "not CCD and gaussian")
    analyze_replicates(~wp_comp_ccd_gauss & wp_comp_ccd_difffromrng, wp_ensg, "CCD and not gaussian")
    
    pd.DataFrame({
            "ENSG":duplicated_ensg,
            "well_plate_pair":[",".join(pair) for pair in duplicated_ensg_pairs],
            "sum(ccd)":duplicated_ensg_ccd,
            "sum(bulk_ccd)":np.array([sum(wp_comp_ccd_gauss[wp_ensg == ensg]) for ensg in duplicated_ensg]),
            "ccd_pair":[",".join([str(wp_comp_ccd_use[u_well_plates == wp][0]) for wp in pair]) for pair in duplicated_ensg_pairs]
        }).to_csv("output/ReplicatedCellCycleDependentProteins.csv")
    
    # Replicates
    for ensg in fileprefixes[ensg_is_duplicated]:
        for file in glob.glob(os.path.join(folder, f"{ensg}*_mvavg.png")):
            shutil.copy(file, os.path.join(replicatefolder, os.path.basename(file)))
            
    # Accounting for biologically CCD ones
    knownccd1 = np.genfromtxt("input/processed/manual/knownccd.txt", dtype='str') # from gene ontology, reactome, cyclebase 3.0, NCBI gene from mcm3
    knownccd2 = np.genfromtxt("input/processed/manual/known_go_ccd.txt", dtype='str') # from GO cell cycle
    knownccd3 = np.genfromtxt("input/processed/manual/known_go_proliferation.txt", dtype='str') # from GO proliferation
    print(f"{len(bioccd)}: number of mitotic structure proteins")
    
    ccd_prots_withmitotic = np.unique(np.concatenate((wp_ensg[ccd_comp], bioccd)))
    total_proteins_minusmitotic = len(ccd_prots_withmitotic) - sum(np.isin(wp_ensg[nonccd_comp], bioccd))
    total_ccd_proteins_withmitotic = len(ccd_prots_withmitotic)
    total_nonccd_proteins_minusmitotic = nonccd_protein_ct - sum(np.isin(wp_ensg[nonccd_comp], bioccd))
    
    overlapping_knownccd1 = sum(np.isin(ccd_prots_withmitotic, np.concatenate((knownccd1, knownccd2, knownccd3))))
    overlapping_knownccd2 = sum(np.isin(ccd_prots_withmitotic, knownccd2))
    overlapping_knownccd3 = sum(np.isin(ccd_prots_withmitotic, knownccd3))

    with open("output/figuresofmerit.txt", "w") as file:
        fom = "--- protein pseudotime\n\n"
        fom += f"{len(np.unique(wp_ensg))} proteins that were expressed and exhibited variations in the U-2 OS cell line were selected" + "\n\n"
        fom += f"present the first evidence of cell cycle association for {overlapping_knownccd1} proteins" + "\n\n"
        fom += f"Based on this analysis, we identified {ccd_protein_ct} out of {len(np.unique(wp_ensg))} proteins ({100 * ccd_protein_ct / len(np.unique(wp_ensg))}%) to have variance in expression levels temporally correlated to cell cycle progression, and for which the cell-cycle explained 8% or more variance in expression than random." + "\n\n"
        fom += f"majority of the proteins analyzed ({100 * total_nonccd_proteins_minusmitotic / len(np.unique(wp_ensg))}%) showed cell-to-cell variations that were largely unexplained by cell cycle progression" + "\n\n"
        fom += f"Of the {total_ccd_proteins_withmitotic} proteins ({ccd_protein_ct} in interphase, {len(bioccd)} in mitotic structures, and {-(total_ccd_proteins_withmitotic - ccd_protein_ct - len(bioccd))} in both sets) identified to correlate to cell cycle progression, {overlapping_knownccd1} ({100 * overlapping_knownccd1 / total_ccd_proteins_withmitotic}%) had a known association to the cell cycle as determined either by a GO BP term ... The remaining {total_ccd_proteins_withmitotic - overlapping_knownccd1} proteins ({100* (total_ccd_proteins_withmitotic - overlapping_knownccd1) / total_ccd_proteins_withmitotic}%)," + "\n\n"
        fom += f"The patterns of variability were investigated for these {len(wp_ensg)} proteins for the population of cells measured for each protein. The mean fold change between the highest and lowest expressing cells per protein was {np.mean(wp_bimodal_fcmaxmin)}." + "\n\n"
        fom += f"We determined that {unimodal_generally_protein_ct} proteins ({100 * unimodal_generally_protein_ct / len(np.unique(wp_ensg))}%) had unimodal intensity distributions, and {bimodal_generally_protein_ct} proteins ({100 * bimodal_generally_protein_ct / len(np.unique(wp_ensg))}%) were found to display bimodality" + "\n\n"
        fom += f"Of {sum(wp_isbimodal_fcpadj_pass)} bimodal samples that were analyzed for cell cycle dependence, {sum(wp_ccd_bimodalonecluster)} were CCD in one cluster ({sum(wp_ccd_bimodalonecluster & wp_comp_ccd_difffromrng)} of these were CCD when analyzed unimodally), and {sum(wp_ccd_bimodaltwocluster)} were CCD in both clusters ({sum(wp_ccd_bimodaltwocluster & wp_comp_ccd_difffromrng)} were also CCD when analyzed unimodally), and the remaining {sum(wp_isbimodal_fcpadj_pass & ~wp_ccd_bimodalonecluster & ~wp_ccd_bimodaltwocluster)} were non-CCD in both clusters."
        print(fom)
        file.write(fom)
       
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
    
    # pickle the results
    utils.np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp) # removed ones passing in only one replicate
    utils.np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp) # removed ones passing in only one replicate
    