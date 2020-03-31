# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 18:52:36 2020

@author: antho
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcolors
import matplotlib.patches as mpatches
import scanpy as sc
import os
import glob
import shutil
import scipy
import scipy.stats
import seaborn as sbn



def ccd_gene_names(id_list_like):
    '''Convert gene ID list to gene name list'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    return gene_info[(gene_info["gene_id"].isin(id_list_like))]["name"]

def values_comp(values_cell, values_nuc, values_cyto, wp_iscell, wp_isnuc, wp_iscyto):
    '''Get the values for the annotated compartment'''
    values_comp = np.empty_like(values_cell)
    values_comp[wp_iscell] = np.array(values_cell)[wp_iscell]
    values_comp[wp_isnuc] = np.array(values_nuc)[wp_isnuc]
    values_comp[wp_iscyto] = np.array(values_cyto)[wp_iscyto]
    return np.array(values_comp) 

def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    # Written by: Olivia Guest github.com/oliviaguest/gini/blob/master/gini.py
    array = array.flatten()
    if np.amin(array) < 0: 
        array -= np.amin(array) # Values cannot be negative
    array = np.sort(array + 0.0000001) # Values must be sorted and nonzero
    index = np.arange(1, array.shape[0] + 1) # Index per array element
    n = array.shape[0] # Number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) # Gini coefficient

def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection'''
    nobs = len(x)
    return np.arange(1, nobs + 1)/float(nobs)

def benji_hoch(alpha, pvals):
    '''benjimini-hochberg multiple testing correction
    source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html'''
    pvals_array = np.array(pvals)
    pvals_array[np.isnan(pvals_array)] = 1 # fail the ones with not enough data
    pvals_sortind = np.argsort(pvals_array)
    pvals_sorted = np.take(pvals_array, pvals_sortind)
    ecdffactor = _ecdf(pvals_sorted)
    reject = pvals_sorted <= ecdffactor*alpha
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    pvals_corrected_BH = np.empty_like(pvals_corrected)

    # deal with sorting
    pvals_corrected_BH[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_BH = np.empty_like(reject)
    reject_BH[pvals_sortind] = reject
    return pvals_corrected_BH, reject_BH

def bonf(alpha, pvals):
    '''Bonferroni multiple testing correction'''
    pvalsarr = np.array(pvals)
    pvalsarr[np.isnan(pvalsarr)] = 1 # fail the ones with not enough data    
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvalsarr, pvals_sortind)
    alphaBonf = alpha / float(len(pvalsarr))
    rejectBonf = pvals_sorted <= alphaBonf
    pvals_correctedBonf = pvals_sorted * float(len(pvalsarr))
    pvals_correctedBonf_unsorted = np.empty_like(pvals_correctedBonf) 
    pvals_correctedBonf_unsorted[pvals_sortind] = pvals_correctedBonf
    rejectBonf_unsorted = np.empty_like(rejectBonf)
    rejectBonf_unsorted[pvals_sortind] = rejectBonf
    return pvals_correctedBonf_unsorted, rejectBonf_unsorted

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

def np_save_overwriting(fn, arr):
    '''Helper function to always overwrite numpy pickles'''
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)
        
def general_boxplot(group_values, group_labels, xlabel, ylabel, title, showfliers, outfile):
    '''Make a boxplot given equal length group_values and group_labels'''
    if len(group_values) != len(group_labels): 
        print("Error: general_boxplot() requires equal length group_values and group_labels.")
        exit(1)
    mmmm = np.concatenate(group_values)
    cccc = np.concatenate([[label] * len(group_values[iii]) for iii, label in enumerate(group_labels)])
    boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=showfliers, color="grey")
    boxplot.set_xlabel(xlabel, size=36,fontname='Arial')
    boxplot.set_ylabel(ylabel, size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    plt.title(title)
    plt.savefig(outfile)
    plt.show()
    plt.close()  

def general_scatter(x, y, xlabel, ylabel, outfile):
    '''Make a general scatterplot with matplotlib'''
    plt.figure(figsize=(10,10))
    plt.scatter(x, y, label="all")
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.legend()
    plt.savefig(outfile)
    plt.show()
    plt.close()
    
def general_scatter_color(x, y, xlabel, ylabel, c, clabel, show_color_bar, title, outfile, cmap="viridis", alpha=1):
    '''Make a general scatterplot with color using matplotlib'''
    plt.figure(figsize=(10,10))
    plt.scatter(x, y, c=c, cmap=cmap, alpha=alpha)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if show_color_bar:
        cb = plt.colorbar()
        cb.set_label(clabel)
    plt.title(title)
    plt.savefig(outfile)
    plt.show()
    plt.close()

def general_histogram(x, xlabel, ylabel, alpha, outfile):
    '''Make a general histogram'''
    plt.hist(x, alpha=alpha)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.show()
    plt.close()

def format_p(p):
    '''3 decimal places, scientific notation'''
    return '{:0.3e}'.format(p)