# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:51:50 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
import SingleCellProteogenomics.MovingAverages

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

def permutation_analysis_protein(curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm,
                                 perc_var_cell_val, perc_var_nuc_val, perc_var_cyto_val,
                                 wp_iscell, wp_isnuc, wp_iscyto,
                                 mvavg_cell, mvavg_nuc, mvavg_cyto):
    '''Randomization analysis of cell cycle dependence: permute cell order and calculate percent variance due to the cell cycle'''
    perms = np.asarray([np.random.permutation(len(curr_pol)) for nnn in np.arange(PERMUTATIONS)])
    curr_comp_norm = np.asarray(curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm)
    curr_comp_percvar = np.asarray(perc_var_cell_val if wp_iscell[i] else perc_var_nuc_val if wp_isnuc[i] else perc_var_cyto_val)
    curr_comp_mvavg = np.asarray(mvavg_cell if wp_iscell[i] else mvavg_nuc if wp_isnuc[i] else mvavg_cyto)
    curr_comp_perm = np.asarray([curr_comp_norm[perm] for perm in perms])
    curr_mt_perm = np.asarray([curr_mt_cell_norm[perm] for perm in perms])
    curr_mvavg_rng_comp = np.apply_along_axis(MovingAverages.mvavg, 1, curr_comp_perm, WINDOW)
    curr_mvavg_rng_mt = np.apply_along_axis(MovingAverages.mvavg, 1, curr_mt_perm, WINDOW)
    curr_percvar_rng_comp = np.var(curr_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm, axis=1)
    curr_percvar_rng_mt = np.var(curr_mvavg_rng_mt, axis=1) / np.var(curr_mt_perm, axis=1)
    return curr_percvar_rng_comp, curr_percvar_rng_mt
    
def cell_cycle_dependence_protein(use_log_ccd, do_remove_outliers,
                                  pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
                                  wp_iscell, wp_isnuc, wp_iscyto,
                                  wp_isbimodal_fcpadj_pass):
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
        perc_var_comp_rng, perc_var_mt_rng = permutation_analysis_protein(curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm,
                                 perc_var_cell_val, perc_var_nuc_val, perc_var_cyto_val, wp_iscell, wp_isnuc, wp_iscyto, mvavg_cell, mvavg_nuc, mvavg_cyto)
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
            MovingAverages.temporal_mov_avg(curr_pol[clust1_idx], curr_comp_norm[clust1_idx], mvavgs_x_clust1[-1], mvavg_clust1, windows1, None, folder, fileprefixes[i] + "_clust1")
            MovingAverages.temporal_mov_avg(curr_pol[clust2_idx], curr_comp_norm[clust2_idx], mvavgs_x_clust2[-1], mvavg_clust2, windows2, None, folder, fileprefixes[i] + "_clust2")
        
        # Make example plots for the randomization trials, but no need for the bimodal ones
        elif np.mean(curr_comp_percvar - curr_percvar_rng_comp) > 0.05:
            median_rng_idx = np.argsort(curr_percvar_rng_comp)
            for iii, idx in enumerate([0,len(curr_percvar_rng_comp)//4,len(curr_percvar_rng_comp)//2,3*len(curr_percvar_rng_comp)//4,len(curr_percvar_rng_comp)-1]):
                MovingAverages.temporal_mov_avg_randomization_example(curr_pol, curr_comp_norm, curr_comp_perm[median_rng_idx[idx]], 
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
        MovingAverages.temporal_mov_avg(curr_pol, curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm, mvavg_xvals,
             mvavg_cell if wp_iscell[i] else mvavg_nuc if wp_isnuc[i] else mvavg_cyto,
             windows, None, folder, fileprefixes[i])
        
        # Uncomment to make the plots for microtubules (takes 10 mins)
        MovingAverages.temporal_mov_avg(curr_pol, curr_mt_cell_norm, mvavg_xvals, mvavg_mt, windows, None, folder_mt, f"{fileprefixes[i]}_mt")
        
        if clusters is not None:
            MovingAverages.temporal_mov_avg(curr_pol, curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm, mvavg_xvals,
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
    
    return wp_comp_ccd_difffromrng, wp_comp_ccd_clust1, wp_comp_ccd_clust2