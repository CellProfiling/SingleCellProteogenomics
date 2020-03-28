# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 19:28:05 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
import scipy.optimize
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

def plot_average_intensities_by_batch(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, mean_mean_mt, wp_iscell, wp_isnuc, wp_iscyto):
    '''Plot average intensities by batch to look at agreement'''
    mean_mean_comp = values_comp(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    firstbatch = np.asarray([not str(p).split("_")[1].startswith("67") for p in u_well_plates])
    plt.hist(np.array(mean_mean_comp)[firstbatch], bins=200, label="firstbatch")
    plt.hist(np.array(mean_mean_comp)[~firstbatch], bins=200, label="secondbatch")
    plt.legend()
    plt.savefig("figures/MeanMeancComp.png")
    plt.show()
    plt.close()
    
    firstbatch = np.asarray([not str(p).split("_")[1].startswith("67") for p in u_well_plates])
    plt.hist(np.array(mean_mean_mt)[firstbatch], bins=200, label="firstbatch")
    plt.hist(np.array(mean_mean_mt)[~firstbatch], bins=200, label="secondbatch")
    plt.legend()
    plt.savefig("figures/MeanMeanMt.png")
    plt.show()
    plt.close()
    

def calculate_variation(use_log, u_well_plates, wp_iscell, wp_isnuc, wp_iscyto, 
                        pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, pol_sort_well_plate_imgnb):
    '''Calculate overall variation of protein staining intensity in single cells'''
    var_cell, var_nuc, var_cyto, var_mt = [],[],[],[] # mean intensity variances per antibody
    cv_cell, cv_nuc, cv_cyto, cv_mt = [], [], [], []
    gini_cell, gini_nuc, gini_cyto, gini_mt = [],[],[],[] # mean intensity ginis per antibody
    var_cell_test_p, var_nuc_test_p, var_cyto_test_p = [],[],[]
    mean_mean_cell, mean_mean_nuc, mean_mean_cyto, mean_mean_mt = [],[],[],[] # mean mean-intensity
    cell_counts = []
    
    wpi_img = []
    var_cell_img, var_nuc_img, var_cyto_img, var_mt_img = [],[],[],[] # mean intensity variances per field of view
    cv_cell_img, cv_nuc_img, cv_cyto_img, cv_mt_img = [], [], [], []
    
    # The variance needs to be calculated separately for each well because they all have different numbers of cells
    for well in u_well_plates:
        curr_well_inds = pol_sort_well_plate==well
        curr_ab_cell = pol_sort_ab_cell[curr_well_inds] if not use_log else np.log10(pol_sort_ab_cell[curr_well_inds])
        curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds] if not use_log else np.log10(pol_sort_ab_nuc[curr_well_inds])
        curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds] if not use_log else np.log10(pol_sort_ab_cyto[curr_well_inds])
        curr_mt_cell = pol_sort_mt_cell[curr_well_inds] if not use_log else np.log10(pol_sort_mt_cell[curr_well_inds])
        curr_fuccigreen=0
        curr_fuccired=0
        
        cell_counts.append(len(curr_ab_cell))
        
        var_cell.append(np.var(curr_ab_cell))
        var_nuc.append(np.var(curr_ab_nuc))
        var_cyto.append(np.var(curr_ab_cyto))
        var_mt.append(np.var(curr_mt_cell))
        
        cv_cell.append(scipy.stats.variation(curr_ab_cell))
        cv_nuc.append(scipy.stats.variation(curr_ab_nuc))
        cv_cyto.append(scipy.stats.variation(curr_ab_cyto))
        cv_mt.append(scipy.stats.variation(curr_mt_cell))
        
        gini_cell.append(gini(curr_ab_cell))
        gini_nuc.append(gini(curr_ab_nuc))
        gini_cyto.append(gini(curr_ab_cyto))
        gini_mt.append(gini(curr_mt_cell))
        
        # Save the mean mean intensities
        mean_mean_cell.append(np.mean(curr_ab_cell))
        mean_mean_nuc.append(np.mean(curr_ab_nuc)) 
        mean_mean_cyto.append(np.mean(curr_ab_cyto))
        mean_mean_mt.append(np.mean(curr_mt_cell))
        
        curr_well_plate_imgnbs = pol_sort_well_plate_imgnb[curr_well_inds]
        curr_wpi_img = []
        curr_var_cell_img, curr_var_nuc_img, curr_var_cyto_img, curr_var_mt_img = [],[],[],[] # mean intensity variances per field of view
        curr_cv_cell_img, curr_cv_nuc_img, curr_cv_cyto_img, curr_cv_mt_img = [], [], [], []
        for wpi in np.unique(curr_well_plate_imgnbs):
            curr_wpis = pol_sort_well_plate_imgnb == wpi
            curr_ab_cell = pol_sort_ab_cell[curr_wpis] if not use_log else np.log10(pol_sort_ab_cell[curr_wpis])
            curr_ab_nuc = pol_sort_ab_nuc[curr_wpis] if not use_log else np.log10(pol_sort_ab_nuc[curr_wpis])
            curr_ab_cyto = pol_sort_ab_cyto[curr_wpis] if not use_log else np.log10(pol_sort_ab_cyto[curr_wpis])
            curr_mt_cell = pol_sort_mt_cell[curr_wpis] if not use_log else np.log10(pol_sort_mt_cell[curr_wpis])
            
            curr_wpi_img.append(wpi)
            
            curr_var_cell_img.append(np.var(curr_ab_cell))
            curr_var_nuc_img.append(np.var(curr_ab_nuc))
            curr_var_cyto_img.append(np.var(curr_ab_cyto))
            curr_var_mt_img.append(np.var(curr_mt_cell))
            
            curr_cv_cell_img.append(scipy.stats.variation(curr_ab_cell))
            curr_cv_nuc_img.append(scipy.stats.variation(curr_ab_nuc))
            curr_cv_cyto_img.append(scipy.stats.variation(curr_ab_cyto))
            curr_cv_mt_img.append(scipy.stats.variation(curr_mt_cell))
        
        wpi_img.append(curr_wpi_img)
        var_cell_img.append(curr_var_cell_img)
        var_nuc_img.append(curr_var_nuc_img)
        var_cyto_img.append(curr_var_cyto_img)
        var_mt_img.append(curr_var_mt_img)
        
        cv_cell_img.append(curr_cv_cell_img)
        cv_nuc_img.append(curr_cv_nuc_img)
        cv_cyto_img.append(curr_cv_cyto_img)
        cv_mt_img.append(curr_cv_mt_img)
    
    print("Plotting average intensities of proteins and microtubules by batch.")
    plot_average_intensities_by_batch(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, mean_mean_mt, wp_iscell, wp_isnuc, wp_iscyto)
    
    print("Making general plots for variance, CV, and gini by compartment")
    var_cell, var_nuc, var_cyto, var_mt = np.array(var_cell), np.array(var_nuc), np.array(var_cyto), np.array(var_mt)
    gini_cell, gini_nuc, gini_cyto, gini_mt = np.array(gini_cell), np.array(gini_nuc), np.array(gini_cyto), np.array(gini_mt)
    cv_cell, cv_nuc, cv_cyto, cv_mt = np.array(cv_cell), np.array(cv_nuc), np.array(cv_cyto), np.array(cv_mt)
    general_boxplot((var_cell, var_cyto, var_nuc, var_mt), ("var_cell", "var_cyto", "var_nuc", "var_mt"),
                "Metacompartment", f"Variance using {'log' if use_log else 'natural'} intensity values",
                "", "figures/VarianceBoxplot.png")
    general_boxplot((cv_cell, cv_cyto, cv_nuc, cv_mt), ("cv_cell", "cv_cyto", "cv_nuc", "cv_mt"),
                "Metacompartment", f"Coeff. of Var. using {'log' if use_log else 'natural'} intensity values",
                "", "figures/CVBoxplot.png")
    general_boxplot((gini_cell, gini_cyto, gini_nuc, gini_mt), ("gini_cell", "gini_cyto", "gini_nuc", "gini_mt"),
                "Metacompartment", f"Gini using {'log' if use_log else 'natural'} intensity values",
                "", "figures/GiniBoxplot.png")

    print("Making general plots for variance, CV, and gini in the compartment the protein localizes to")
    mean_mean_comp = values_comp(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    cv_comp = values_comp(cv_cell, cv_nuc, cv_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    gini_comp = values_comp(gini_cell, gini_nuc, gini_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    var_comp = values_comp(var_cell, var_nuc, var_cyto, wp_iscell, wp_isnuc, wp_iscyto)
    general_scatter(var_comp, var_mt, f"var_comp", f"var_mt", f"figures/var_comp_mt.png")
    general_scatter(cv_comp, cv_mt, f"cv_comp", f"cv_mt", f"figures/cv_comp_mt.png")
    general_scatter(gini_comp, gini_mt, f"gini_comp", f"gini_mt", f"figures/gini_comp_mt.png")
    general_scatter(var_comp, mean_mean_comp, f"var_comp", f"{'log10' if use_log else 'natural'} intensity", f"figures/VarianceVsIntensityComp.png")

    print("Comparing image to sample variance")
    var_comp_img = values_comp(var_cell_img, var_nuc_img, var_cyto_img, wp_iscell, wp_isnuc, wp_iscyto)
    gini_comp_img = values_comp(gini_cell_img, gini_nuc_img, gini_cyto_img, wp_iscell, wp_isnuc, wp_iscyto)
    cv_comp_img = values_comp(cv_cell_img, cv_nuc_img, cv_cyto_img, wp_iscell, wp_isnuc, wp_iscyto)
    general_scatter(np.concatenate([[var_comp[i]] * len(vvv) for i, vvv in enumerate(var_comp_img)]), np.concatenate(var_comp_img), "variance within compartment", "variance for each image", "figures/VarianceByImage.png")
    general_scatter(np.concatenate([[gini_comp[i]] * len(vvv) for i, vvv in enumerate(gini_comp_img)]), np.concatenate(gini_comp_img), "gini within compartment", "gini for each image", "figures/GiniByImage.png")
    general_scatter(np.concatenate([[cv_comp[i]] * len(vvv) for i, vvv in enumerate(cv_comp_img)]), np.concatenate(cv_comp_img), "cv within compartment", "cv for each image", "figures/CVByImage.png")
    print(np.concatenate(wpi_img)[np.argmax(np.concatenate(var_comp_img))] + ": the image with the max variance")

    plt.hist(np.concatenate([vvv / var_comp[i] for i, vvv in enumerate(var_comp_img)]))
    plt.show();plt.close();
    high_var_img = np.concatenate(wpi_img)[np.concatenate([vvv > 4 * var_comp[i] for i, vvv in enumerate(var_comp_img)])]
    print(f"{high_var_img}: the images with greater than 4x the variance of the whole sample")
    
    norm_cv_img = np.concatenate([vvv / cv_comp[i] for i, vvv in enumerate(cv_comp_img)])
    plt.hist(norm_cv_img)
    plt.show();plt.close();
    cutoff = np.mean(norm_cv_img) + 3 * np.std(norm_cv_img)
    high_cv_img = np.concatenate(wpi_img)[norm_cv_img > cutoff]
    print(f"{high_cv_img}: the images with greater than 4x the variance of the whole sample")
    
    np.intersect1d(high_var_img, high_cv_img)
    
    # Pickle and return main results
    np_save_overwriting("output/pickles/mean_mean_comp.npy", mean_mean_comp)
    np_save_overwriting("output/pickles/cv_comp.npy", cv_comp)
    np_save_overwriting("output/pickles/gini_comp.npy", gini_comp)
    np_save_overwriting("output/pickles/var_comp.npy", var_comp)
    np_save_overwriting("output/pickles/cv_cell.npy", cv_cell)
    np_save_overwriting("output/pickles/gini_cell.npy", gini_cell)
    np_save_overwriting("output/pickles/var_cell.npy", var_cell)
    return mean_mean_comp, var_comp, gini_comp, cv_comp, var_cell, gini_cell, cv_cell

