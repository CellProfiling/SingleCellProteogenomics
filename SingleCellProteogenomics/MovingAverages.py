# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 21:21:59 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
from SingleCellProteogenomics.FucciCellCycle import FucciCellCycle

TOT_LEN = FucciCellCycle().TOT_LEN

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
    '''Returns indices of outliers to remove'''
    max_cutoff = np.mean(values) + 5 * np.std(values)
    min_cutoff = np.mean(values) - 5 * np.std(values)
    return (values < max_cutoff) & (values > min_cutoff)

def remove_outliers(values, return_values):
    '''Remove outliers on "values" and return "return_values" based on that filter'''
    return return_values[remove_outliers_idx(values)]

## MOVING AVERAGE PLOTS

def temporal_mov_avg_protein(curr_pol, curr_ab_norm, mvavg_xvals, mvavg_yvals, windows, clusters, folder, fileprefix):
    plt.close()
    outfile = os.path.join(folder,fileprefix+'_mvavg.pdf')
    if os.path.exists(outfile): return
    plt.figure(figsize=(5,5))
    mvperc = mvpercentiles(curr_ab_norm[windows])
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
        
def temporal_mov_avg_rna(fucci_time, curr_ab_norm, mvavg_xvals, mvavg_yvals, windows, folder, fileprefix):
    plt.close()
    outfile = os.path.join(folder,fileprefix+'_mvavg.pdf')
    if os.path.exists(outfile): return
    plt.figure(figsize=(5,5))
    mvperc = mvpercentiles(curr_ab_norm[windows])
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[0], mvperc[-1], color="bisque", label="10th & 90th Percentiles")
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc[1], mvperc[-2], color="orange", alpha=0.7, label="25th & 75th Percentiles")
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_yvals, color="tab:orange", label="Mean Intensity")
    filteridx = remove_outliers_idx(curr_ab_norm)
    plt.scatter(fucci_time[filteridx] * TOT_LEN, curr_ab_norm[filteridx], c='darkorange', alpha=0.2)
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

    
