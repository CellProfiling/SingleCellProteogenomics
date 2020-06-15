# -*- coding: utf-8 -*-
"""
Clustering the cells into cell cycle phases using FUCCI cell cycle markers:
    - The log intensities of FUCCI markers, CDT1 and GMNN, are centered in log-log space to remove batch effects
    - Gaussian clustering into three phases (G1, S-transition, G2) are used to distinguish phases
    - Differences in protein abundance evaluated for significance using Kruskal-Wallis tests
    - This constitutes a "mock-bulk" analysis using single cell data

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
import numpy as np
import sklearn.mixture
import seaborn as sbn
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

def zero_center_fucci(green_fucci, red_fucci, u_plate, well_plate, plate):
    '''Zero center and rescale FUCCI data in the log space'''
    log_green_fucci, log_red_fucci = np.log10(green_fucci), np.log10(red_fucci)
    wp_p_dict = dict([(str(p), plate == p) for p in u_plate])
    logmed_green_fucci_p = dict([(str(p), np.log10(np.median(green_fucci[wp_p_dict[str(p)]]))) for p in u_plate])
    logmed_red_fucci_p = dict([(str(p), np.log10(np.median(red_fucci[wp_p_dict[str(p)]]))) for p in u_plate])
    logmed_green_fucci = np.array([logmed_green_fucci_p[wp.split("_")[1]] for wp in well_plate])
    logmed_red_fucci = np.array([logmed_red_fucci_p[wp.split("_")[1]] for wp in well_plate])
    log_green_fucci_zeroc = np.array(log_green_fucci) - logmed_green_fucci
    log_red_fucci_zeroc = np.array(log_red_fucci) - logmed_red_fucci
    log_green_fucci_zeroc_rescale = (log_green_fucci_zeroc - np.min(log_green_fucci_zeroc)) / np.max(log_green_fucci_zeroc)
    log_red_fucci_zeroc_rescale = (log_red_fucci_zeroc - np.min(log_red_fucci_zeroc)) / np.max(log_red_fucci_zeroc)
    fucci_data = np.column_stack([log_green_fucci_zeroc_rescale,log_red_fucci_zeroc_rescale])
    result = (log_green_fucci, log_red_fucci,
              log_green_fucci_zeroc, log_red_fucci_zeroc,
              log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale,
              fucci_data)
    
    # Pickle the results
    utils.np_save_overwriting("output/pickles/log_green_fucci_zeroc.npy", log_green_fucci_zeroc)
    utils.np_save_overwriting("output/pickles/log_red_fucci_zeroc.npy", log_red_fucci_zeroc)
    utils.np_save_overwriting("output/pickles/log_green_fucci_zeroc_rescale.npy", log_green_fucci_zeroc_rescale)
    utils.np_save_overwriting("output/pickles/log_red_fucci_zeroc_rescale.npy", log_red_fucci_zeroc_rescale)
    utils.np_save_overwriting("output/pickles/fucci_data.npy", fucci_data)

    return result

def gaussian_boxplot_result(g1, s, g2, outfolder, ensg):
    '''Boxplot for intensities within each cell cycle phase'''
    if not os.path.exists(f"{outfolder}_png"): os.mkdir(f"{outfolder}_png")
    if not os.path.exists(f"{outfolder}_pdf"): os.mkdir(f"{outfolder}_pdf")
    mmmm = np.concatenate((g1, s, g2))
    cccc = (["G1"] * len(g1))
    cccc.extend(["G1/S"] * len(s))
    cccc.extend(["G2"] * len(g2))
    boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=False, color="grey")
    boxplot.set_xlabel("", size=36,fontname='Arial')
    boxplot.set_ylabel("Normalized Mean Intensity", size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    plt.ylim(0,1)
    plt.title("")
    plt.savefig(f"{outfolder}_png/GaussianClusteringProtein_{ensg}.png")
    plt.savefig(f"{outfolder}_pdf/GaussianClusteringProtein_{ensg}.pdf")
    plt.close()
    
def gaussian_clustering(log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale):
    '''Perform gaussian clustering of FUCCI data into 3 phases: G1, S, G2'''
    gaussian = sklearn.mixture.GaussianMixture(n_components=3, random_state=1, max_iter=500)
    cluster_labels = gaussian.fit_predict(np.array([log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale]).T)
    clusternames = ["G2","S-ph","G1"]
    for cluster in range(3):
        plt.hist2d(log_green_fucci_zeroc_rescale[cluster_labels == cluster],log_red_fucci_zeroc_rescale[cluster_labels == cluster],bins=200)
        plt.title(f"Gaussian clustered data, {clusternames[cluster]}")
        plt.xlabel("Log10 Green Fucci Intensity")
        plt.ylabel("Log10 Red Fucci Intensity")
        plt.savefig(f"figures/FucciPlotProteinIFData_unfiltered_Gauss{cluster}.png")
        plt.show()
        plt.close()
    return cluster_labels

def get_phase_strings(is_g1, is_sph, is_g2):
    '''Make strings to represent the metacompartment'''
    phasestring = np.array(["G1"] * len(is_g1))
    phasestring[is_sph] = "S" 
    phasestring[is_g2] = "G2"
    return phasestring

def gaussian_clustering_analysis(alpha_gauss, doGeneratePlots, g1, sph, g2, 
             wp_ensg, well_plate, u_well_plates, ab_cell, ab_nuc, ab_cyto, mt_cell, wp_iscell, wp_isnuc, wp_iscyto):
    '''Analyze the results of Gaussian clustering of FUCCI data for each protein antibody staining'''
    wp_cell_kruskal, wp_nuc_kruskal, wp_cyto_kruskal, wp_mt_kruskal = [],[],[],[]
    curr_wp_phases = []
    mockbulk_phases = np.array(["  "] * len(ab_cell))
    fileprefixes = np.array([f"{ensg}_{sum(wp_ensg[:ei] == ensg)}" for ei, ensg in enumerate(wp_ensg)])
    for iii, wp in enumerate(u_well_plates):
        curr_well_inds = well_plate==wp
        curr_wp_g1 = curr_well_inds & g1
        curr_wp_sph = curr_well_inds & sph
        curr_wp_g2 = curr_well_inds & g2
        mockbulk_phases[np.arange(len(ab_cell))[curr_well_inds]] = get_phase_strings(g1[curr_well_inds], sph[curr_well_inds], g2[curr_well_inds])
        curr_wp_phases.append(get_phase_strings(g1[curr_well_inds], sph[curr_well_inds], g2[curr_well_inds]))
        wp_cell_kruskal.append(scipy.stats.kruskal(ab_cell[curr_wp_g1], ab_cell[curr_wp_sph], ab_cell[curr_wp_g2])[1])
        wp_nuc_kruskal.append(scipy.stats.kruskal(ab_nuc[curr_wp_g1], ab_nuc[curr_wp_sph], ab_nuc[curr_wp_g2])[1])
        wp_cyto_kruskal.append(scipy.stats.kruskal(ab_cyto[curr_wp_g1], ab_cyto[curr_wp_sph], ab_cyto[curr_wp_g2])[1])
        wp_mt_kruskal.append(scipy.stats.kruskal(mt_cell[curr_wp_g1], mt_cell[curr_wp_sph], mt_cell[curr_wp_g2])[1])
        max_val_for_norm = np.max(ab_cell[curr_well_inds] if wp_iscell[iii] else ab_nuc[curr_well_inds] if wp_isnuc[iii] else ab_cyto[curr_well_inds])
        max_mt_for_norm = np.max(mt_cell[curr_well_inds])
        if doGeneratePlots:
            gaussian_boxplot_result(
                    (ab_cell[curr_wp_g1] if wp_iscell[iii] else ab_nuc[curr_wp_g1] if wp_isnuc[iii] else ab_cyto[curr_wp_g1]) / max_val_for_norm,
                    (ab_cell[curr_wp_sph] if wp_iscell[iii] else ab_nuc[curr_wp_sph] if wp_isnuc[iii] else ab_cyto[curr_wp_sph]) / max_val_for_norm,
                    (ab_cell[curr_wp_g2] if wp_iscell[iii] else ab_nuc[curr_wp_g2] if wp_isnuc[iii] else ab_cyto[curr_wp_g2]) / max_val_for_norm,
                    "figures/GaussianBoxplots", fileprefixes[iii])
            gaussian_boxplot_result(
                mt_cell[curr_wp_g1] / max_mt_for_norm,
                mt_cell[curr_wp_sph] / max_mt_for_norm,
                mt_cell[curr_wp_g2] / max_mt_for_norm,
                "figures/GaussianBoxplots_mt", f"{fileprefixes[iii]}_mt")
        
    # multiple testing correction for protein of interest
    wp_comp_kruskal_gaussccd_p = utils.values_comp(wp_cell_kruskal, wp_nuc_kruskal, wp_cyto_kruskal, wp_iscell, wp_isnuc, wp_iscyto)
    wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp = utils.benji_hoch(alpha_gauss, wp_comp_kruskal_gaussccd_p)
    utils.np_save_overwriting("output/pickles/wp_comp_kruskal_gaussccd_adj.npy", wp_comp_kruskal_gaussccd_adj)
    utils.np_save_overwriting("output/pickles/wp_pass_kruskal_gaussccd_bh_comp.npy", wp_pass_kruskal_gaussccd_bh_comp)

    # multiple testing correction for microtubules
    wp_mt_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_mt = utils.benji_hoch(alpha_gauss, wp_mt_kruskal) 
    utils.np_save_overwriting("output/pickles/wp_mt_kruskal_gaussccd_adj.npy", wp_mt_kruskal_gaussccd_adj)
    utils.np_save_overwriting("output/pickles/wp_pass_gaussccd_bh_mt.npy", wp_pass_gaussccd_bh_mt)
    
    # save the phase information
    utils.np_save_overwriting("output/pickles/curr_wp_phases.npy", np.array(curr_wp_phases))
    utils.np_save_overwriting("output/pickles/mockbulk_phases.npy", np.array(mockbulk_phases))

    print(f"{len(wp_pass_kruskal_gaussccd_bh_comp)}: number of genes tested")
    print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp)}: number of passing genes at {alpha_gauss*100}% FDR in compartment")

    return wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp, wp_mt_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_mt

def address_replicates(alpha_gauss, wp_pass_kruskal_gaussccd_bh_comp, wp_ensg, wp_ab, u_well_plates):
    '''Look for replicated protein samples and antibody stainings'''
    # address gene redundancy
    wp_ensg_counts = np.array([sum([1 for eeee in wp_ensg if eeee == ensg]) for ensg in wp_ensg])
    ensg_is_duplicated = wp_ensg_counts > 1
    duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
    duplicated_ensg_pairs = [u_well_plates[wp_ensg == ensg] for ensg in duplicated_ensg]
    print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp[~ensg_is_duplicated])}: number of passing genes at {alpha_gauss*100}% FDR in compartment (no replicate)")
    duplicated_ensg_ccd = np.array([sum(wp_pass_kruskal_gaussccd_bh_comp[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(duplicated_ensg_ccd == 2)}: number of CCD genes shown to be CCD in both replicates")
    print(f"{sum(duplicated_ensg_ccd == 1)}: number of CCD genes shown to be CCD in just one replicate")
    print(f"{sum(duplicated_ensg_ccd == 0)}: number of CCD genes shown to be non-CCD in both replicate")
    
    # any antibody redundancy?
    wp_ab_counts = np.array([sum([1 for aaaa in wp_ab if aaaa == ab]) for ab in wp_ab])
    ab_is_duplicated = wp_ab_counts > 1
    duplicated_ab = np.unique(wp_ab[ab_is_duplicated])
    print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp[~ab_is_duplicated])}: number of passing antibodies at {alpha_gauss*100}% FDR in compartment (no replicate)")
    duplicated_ab_ccd = np.array([sum(wp_pass_kruskal_gaussccd_bh_comp[wp_ab == ab]) for ab in duplicated_ab])
    print(f"{sum(duplicated_ab_ccd == 2)}: number of duplicated antibodies shown to be CCD in both replicates")
    print(f"{sum(duplicated_ab_ccd == 1)}: number of duplicated antibodies shown to be CCD in just one replicate")
    print(f"{sum(duplicated_ab_ccd == 0)}: number of duplicated antibodies shown to be non-CCD in both replicate")