# -*- coding: utf-8 -*-
"""
Methods for assessing bimodality of protein intensity distributions.

Distinct high- and low-expressing cell populations that had no correlation to cell division time
were evaluated separately for cell cycle dependence in subsequent analysis.

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
import sklearn.mixture

def identify_bimodal_intensity_distributions(u_well_plates, wp_ensg,
             pol_sort_well_plate, pol_sort_norm_rev, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell,
             wp_iscell, wp_isnuc, wp_iscyto):
    '''
    Some proteins display bimodal intensity distributions. 
    This method seeks to identify distributions with high- and low-expressing cells, 
        so that they may be assessed for CCD independently in `ProteinCellCycleDependence.py`.
    '''
    wp_bimodal_cluster_idxs = []
    wp_bimodal_diffmeans = []
    wp_bimodal_fcmeans = []
    wp_bimodal_fcmaxmin = []
    wp_bimodal_clusterlabels = []
    wp_isbimodal_p = []
    wp_timebimodal_p = []
    wp_intensities = []

    # Use Gaussian clustering to investigate if there is bimodality
    gaussian = sklearn.mixture.GaussianMixture(n_components=2, random_state=1, max_iter=500)
    for i, well in enumerate(u_well_plates):
        curr_well_inds = pol_sort_well_plate==well # the reversal isn't really helpful here
        curr_pol = pol_sort_norm_rev[curr_well_inds]
        curr_ab_cell = pol_sort_ab_cell[curr_well_inds]
        curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds]
        curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds]
        curr_mt_cell = pol_sort_mt_cell[curr_well_inds]
    
        # Normalize mean intensities, normalized for display
        curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
        curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
        curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
        curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
        curr_comp_norm = np.asarray(curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm)
        wp_intensities.append(curr_comp_norm)
    
        cluster_labels = gaussian.fit_predict(curr_comp_norm.reshape(1, -1).T)
    #    cluster_labels = gaussian.fit_predict(np.array([curr_pol, curr_comp_norm]).T)
        wp_bimodal_clusterlabels.append(cluster_labels)
        c1 = cluster_labels == 0
        c2 = cluster_labels == 1
        wp_bimodal_cluster_idxs.append([c1, c2])
        wp_bimodal_diffmeans.append(np.mean(curr_comp_norm[c2]) - np.mean(curr_comp_norm[c1]))
        wp_bimodal_fcmeans.append(np.mean(curr_comp_norm[c2]) / np.mean(curr_comp_norm[c1]))
        wp_bimodal_fcmaxmin.append(np.max(curr_comp_norm) / np.min(curr_comp_norm))
        
        # Use a kruskal-wallis test to assess whether there's a significant difference of intensities between clusters
        k, p = scipy.stats.kruskal(curr_comp_norm[c1], curr_comp_norm[c2])
        wp_isbimodal_p.append(p)
        
        # Use a kruskal-wallis test to assess whether there's (not) a sigificant difference in pseudotime between clusters,
        # since strongly CCD proteins will produce bimodal intensity distributions that should still be assessed as one population
        k, p = scipy.stats.kruskal(curr_pol[c1], curr_pol[c2])
        wp_timebimodal_p.append(p)
    
    # Multiple testing corrections
    wp_isbimodal_padj, wp_isbimodal_pass = benji_hoch(0.01, wp_isbimodal_p)
    wp_timebimodal_padj, wp_timebimodal_pass = benji_hoch(0.01, wp_timebimodal_p)
    
    wp_enoughcellsinbothclusters = np.array([sum(c1[0]) > 50 and sum(c1[1]) > 50 for c1 in wp_bimodal_cluster_idxs])
    wp_isbimodal_generally = (np.abs(np.log(wp_bimodal_fcmeans) / np.log(2)) > 1) & wp_isbimodal_pass
    wp_isbimodal_fcpadj_pass = (np.abs(np.log(wp_bimodal_fcmeans) / np.log(2)) > 1) & wp_isbimodal_pass & ~wp_timebimodal_pass & wp_enoughcellsinbothclusters
    print(f"{sum(~wp_isbimodal_generally)}: number of proteins displaying unimodal distributions ({sum(~wp_isbimodal_generally)/len(wp_isbimodal_generally)}%)")
    print(f"{sum(wp_isbimodal_generally)}: number of proteins displaying bimodal distributions ({sum(wp_isbimodal_generally)/len(wp_isbimodal_generally)}%)")
    
    # Show that the intensity measurements are reasonable for these bimodal samples
    plt.hist(np.concatenate(np.array(wp_intensities)[wp_isbimodal_generally]))
    plt.xlabel("Mean intensity")
    plt.ylabel("Count")
    plt.title("Intensities of Cells within Bimodal Distributions\nAre Similar to those Overall")
    plt.show();plt.close()

    print("Illustrate the significantly distinct high- and low-expressing cell populations")
    plt.scatter(np.log(wp_bimodal_fcmeans) / np.log(2), -np.log10(wp_isbimodal_padj), c=wp_isbimodal_generally, alpha=0.5, cmap="bwr_r")
    plt.xlabel("Log2 Fold Change Between Gaussian Clusters")
    plt.ylabel("-Log10 Adj. p-Value for Difference Between Clusters")
    plt.savefig("figures/BimodalSignificance_GeneralBimodality.png")
    plt.show();plt.close()

    print("Illustrate the significantly distinct high- and low-expressing cell populations")
    print("with no difference in pseudotime. These are evaluated separately for CCD.")
    plt.scatter(np.log(wp_bimodal_fcmeans) / np.log(2), -np.log10(wp_isbimodal_padj), c=wp_isbimodal_fcpadj_pass, alpha=0.5, cmap="bwr_r")
    plt.xlabel("Log2 Fold Change Between Gaussian Clusters")
    plt.ylabel("-Log10 Adj. p-Value for Difference Between Clusters")
    plt.savefig("figures/BimodalSignificance.png")
    plt.savefig("figures/BimodalSignificance.pdf")
    plt.show();plt.close()
    
    print("Illustrate the samples with sufficient cell count for CCD evaluation of high- and low-expressing cell populations.")
    plt.scatter([sum(c1[0]) for c1 in wp_bimodal_cluster_idxs], [sum(c1[1]) for c1 in wp_bimodal_cluster_idxs], c=wp_enoughcellsinbothclusters, alpha=0.5, cmap="bwr_r")
    plt.xlabel("Cell Count, Cluster 1")
    plt.ylabel("Cell Count, Cluster 2")
    plt.savefig("figures/BimodalCellCount.png")
    plt.show();plt.close()
    
    return wp_isbimodal_fcpadj_pass, wp_bimodal_cluster_idxs, wp_isbimodal_generally, wp_bimodal_fcmaxmin