# -*- coding: utf-8 -*-
"""
Methods for calculating moving averages of single-cell protein and RNA expression over the cell division cycle.

@author: devinsullivan
@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
from SingleCellProteogenomics.FucciCellCycle import FucciCellCycle

TOT_LEN = FucciCellCycle().TOT_LEN  # Object representing FUCCI cell cycle phase durations

## MOVING AVERAGE FUNCTIONS

def mvavg(yvals, mv_window):
    '''Calculate the moving average'''
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')
def mvpercentiles(yvals_binned):
    '''Calculate moving percentiles given moving-window binned values'''
    return np.percentile(yvals_binned, [10, 25, 50, 75, 90], axis=1)
def mvavg_perc_var(yvals,mv_window):
    '''Calculate moving average and the percent variance that can be attributed to cell cycle'''
    yval_avg = np.convolve(yvals,np.ones((mv_window,))/mv_window, mode='valid')
    return np.var(yval_avg)/np.var(yvals), yval_avg
def mvmed(yvals_binned):
    '''Calculate the moving median given moving-window binned values'''
    return np.median(yvals_binned, axis=1)
def mvmed_perc_var(yvals, windows):
    '''Calculate moving median and the percent variance that can be attributed to the cell cycle'''
    yval_avg = mvmed(yvals[windows])
    return np.var(yval_avg) / np.var(yvals), yval_avg

def remove_outliers_idx(values):
    '''Returns indices of outliers to keep'''
    max_cutoff = np.mean(values) + 5 * np.std(values)
    min_cutoff = np.mean(values) - 5 * np.std(values)
    return (values < max_cutoff) & (values > min_cutoff)

def remove_outliers(values, return_values):
    '''Remove outliers on "values" and return "return_values" based on that filter'''
    return return_values[remove_outliers_idx(values)]

## MOVING AVERAGE PLOTS

def temporal_mov_avg_protein(curr_pol, curr_ab_norm, mvavg_xvals, mvavg_yvals, mvperc, clusters, folder, fileprefix):
    '''
    Generates a moving average plot for one protein
    Input: Antibody intensity measurements for the current protein
    Output: Moving average plot with scatter of protein abundance measurement for each cell over the cell cycle
    '''
    plt.close()
    outfile = os.path.join(folder,fileprefix+'_mvavg.pdf')
    if os.path.exists(outfile): return
    plt.figure(figsize=(5,5))
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[0], mvperc[-1], color="lightsteelblue", label="10th & 90th Percentiles")
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[1], mvperc[-2], color="steelblue", label="25th & 75th Percentiles")
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_yvals, color="blue", label="Mean Intensity")
    plt.scatter(curr_pol * TOT_LEN, curr_ab_norm, c=(clusters if clusters is not None else 'b'), cmap="bwr_r" if clusters is not None else None, alpha=0.3)
    plt.xlabel('Cell Cycle Time, hrs')
    plt.ylabel(fileprefix + ' Protein Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
#    plt.legend(fontsize=14)
    plt.ylim(0, 1)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)
    plt.close()
        
def temporal_mov_avg_rna(fucci_time, curr_rna_norm, mvavg_xvals, mvavg_yvals, mvperc, folder, fileprefix):
    '''
    Generates a moving average plot for one transcript
    Input: Antibody intensity measurements for the current gene
    Output: Moving average plot with scatter of gene expression measurement for each cell over the cell cycle
    '''
    plt.close()
    outfile = os.path.join(folder,fileprefix+'_mvavg.pdf')
    if os.path.exists(outfile): return
    plt.figure(figsize=(5,5))
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[0], mvperc[-1], color="bisque", label="10th & 90th Percentiles")
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[1], mvperc[-2], color="orange", alpha=0.7, label="25th & 75th Percentiles")
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_yvals, color="tab:orange", label="Mean Intensity")
    filteridx = remove_outliers_idx(curr_rna_norm)
    plt.scatter(fucci_time[filteridx] * TOT_LEN, curr_rna_norm[filteridx], c='darkorange', alpha=0.2)
    plt.xlabel('Cell Cycle Time, hrs')
    plt.ylabel(fileprefix + ' RNA Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.ylim(0, 1)
    #    plt.legend(fontsize=14)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)
    plt.close()
    
def temporal_mov_avg_randomization_example_protein(curr_pol, curr_ab_norm, curr_ab_norm_rng, mvavg_xvals, mvavg_yvals, mvavg_yvals_rng, folder, fileprefix, rna_or_protein = "Protein"):
    '''
    Generates a moving average plot for one protein illustrating the protein measurements (blue) and randomization of those measurements (red)
    Input: Antibody intensity measurements for the current gene; results from randomization of the cell order
    Output: Moving average plot with original measurements and results of randomization
    '''
    outfile = os.path.join(folder,fileprefix+'_mvavg.png')
    if os.path.exists(outfile): 
        return
    plt.figure(figsize=(5,5))
    sample_color = "blue" if rna_or_protein == "Protein" else "tab:orange"
    randomization_color = "red" if rna_or_protein == "Protein" else "purple"
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_yvals, color=sample_color, label="Mean Intensity")
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_yvals_rng, color=randomization_color, label="Mean Intensity, Randomized")
    plt.scatter(curr_pol * TOT_LEN, curr_ab_norm, c=sample_color, alpha=0.2, label="Normalized Intensity")
    plt.scatter(curr_pol * TOT_LEN, curr_ab_norm_rng, c=randomization_color, alpha=0.2, label="Normalized Intensity, Randomized")
    plt.xlabel('Cell Cycle Time, hrs')
    plt.ylabel(fileprefix.split("_")[0] + f' {rna_or_protein} Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
#    plt.legend(fontsize=14)
    plt.ylim(0, 1)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)
    plt.close()

def fix_nans(binned_values):
    '''Custom fix for nan values during binning.'''
    for i,val in enumerate(binned_values):
        do_fix = np.isnan(val)
        if do_fix:
            is_end = i == len(binned_values)-1
            if is_end:
                prevval = i-1
                nextval = i+1
            else:
                prevval = binned_values[i-1]
                nextval = binned_values[i+1]
            prevval = np.nan
            jind = i
            count = 1
            while np.isnan(prevval) and count<len(binned_values):
                if jind>0:
                    jind = jind-1
                    prevval = binned_values[jind]
                elif jind==0:
                    jind = len(binned_values)-1
                    prevval = binned_values[jind]
                count = count+1
            nextval = np.nan
            jind = i
            count = 1
            while np.isnan(nextval) and count<len(binned_values):
                if jind<(len(binned_values)-1):
                    jind = jind+1
                    nextval = binned_values[jind]
                elif jind==(len(binned_values)-1):
                    # print('doing the end of list')
                    jind = 0
                    nextval = binned_values[jind]
                count = count+1
            binned_values[i] = np.mean([prevval,nextval])
    return binned_values

def bin_values(nbins, u_well_plates, pol_sort_norm_rev, pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, 
               pol_sort_ab_cyto, pol_sort_mt_cell, wp_iscell, wp_isnuc, wp_iscyto, do_normalize=True):
    '''Compute protein expression values, binned over pseudotime into `nbins` number of bins.'''
    xvals = np.linspace(0,1,num=nbins)
    wp_max_pol = []
    wp_binned_values = []
    for i, well in enumerate(u_well_plates):
        curr_well_inds = pol_sort_well_plate==well
        curr_pol = pol_sort_norm_rev[curr_well_inds]
    #    curr_fred = pol_sort_fred[curr_well_inds]
    #    curr_fgreen = pol_sort_fgreen[curr_well_inds]
        curr_ab_cell, curr_ab_nuc, curr_ab_cyto, curr_mt_cell = pol_sort_ab_cell[curr_well_inds], pol_sort_ab_nuc[curr_well_inds],pol_sort_ab_cyto[curr_well_inds], pol_sort_mt_cell[curr_well_inds]
    
        # Normalize FUCCI colors & mean intensities, normalized for display
        if do_normalize:
        #    curr_fred_norm = curr_fred / np.max(curr_fred)
        #    curr_fgreen_norm = curr_fgreen / np.max(curr_fgreen)
            curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
            curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
            curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
            curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
        else:
            # curr_fred_norm, curr_fgreen_norm = curr_fred, curr_fgreen
            curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm = curr_ab_cell, curr_ab_nuc, curr_ab_cyto, curr_mt_cell
        
        # Compute binned values
        binned_values = []
        prev_xval = 0
        for xval in xvals:
            if xval==0:
                prev_xval = xval
                continue
            curr_ab_norm = curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm
            binned_values.append(np.median(curr_ab_norm[(curr_pol < xval) & (curr_pol >= prev_xval)]))
            prev_xval = xval
            
        if do_normalize:
            binned_values = binned_values/np.nanmax(binned_values)
        binned_values = fix_nans(binned_values)
        max_loc = np.nanargmax(binned_values)
        if np.isnan(xvals[max_loc]): print('what')
        
        wp_max_pol.append(xvals[max_loc])
        wp_binned_values.append(binned_values)
    return wp_max_pol, wp_binned_values, xvals
