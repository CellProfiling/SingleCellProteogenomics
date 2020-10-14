# -*- coding: utf-8 -*-
"""
Utility functions for:
    - Statistics
    - Saving results
    - Plotting
    - Gene ID - Gene name conversions

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

def values_comp(values_cell, values_nuc, values_cyto, wp_iscell, wp_isnuc, wp_iscyto):
    '''Get the values for the annotated compartment'''
    values_comp = np.empty_like(values_cell)
    values_comp[wp_iscell] = np.array(values_cell)[wp_iscell]
    values_comp[wp_isnuc] = np.array(values_nuc)[wp_isnuc]
    values_comp[wp_iscyto] = np.array(values_cyto)[wp_iscyto]
    return np.array(values_comp) 

def np_save_overwriting(fn, arr):
    '''Helper function to always overwrite numpy pickles'''
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)

## STATISTICS HELPERS

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

## PLOTTING HELPERS

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

def general_boxplot_setup(group_values, group_labels, xlabel, ylabel, title, showfliers, ylim=()):
    '''Set up a boxplot given equal length group_values and group_labels'''
    if len(group_values) != len(group_labels): 
        print("Error: general_boxplot() requires equal length group_values and group_labels.")
        exit(1)
    mmmm = np.concatenate(group_values)
    cccc = np.concatenate([[label] * len(group_values[iii]) for iii, label in enumerate(group_labels)])
    boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=showfliers, color="grey")
    boxplot.set_xlabel(xlabel, size=36,fontname='Arial')
    boxplot.set_ylabel(ylabel, size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    if len(ylim) > 0: boxplot.set(ylim=ylim)
    plt.title(title)
    return cccc, mmmm, boxplot

def general_boxplot(group_values, group_labels, xlabel, ylabel, title, showfliers, outfile, ylim=()):
    '''Make a boxplot given equal length group_values and group_labels'''
    general_boxplot_setup(group_values, group_labels, xlabel, ylabel, title, showfliers, ylim)
    plt.savefig(outfile)
    plt.show()
    plt.close()

def boxplot_with_stripplot(group_values, group_labels, xlabel, ylabel, title, showfliers, outfile, alpha=0.3, size=5, jitter=0.25, ylim=()):
    plt.figure(figsize=(10,10))
    cccc, mmmm, ax = general_boxplot_setup(group_values, group_labels, xlabel, ylabel, title, showfliers, ylim)
    boxplot = sbn.stripplot(x=cccc, y=mmmm, alpha=alpha, color=".3", size=size, jitter=jitter)
    plt.savefig(outfile)
    plt.show()
    plt.close()

def general_scatter(x, y, xlabel, ylabel, outfile, showLegend=True):
    '''Make a general scatterplot with matplotlib'''
    plt.figure(figsize=(10,10))
    plt.scatter(x, y, label="all")
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    if showLegend: 
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

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values. 
    https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)

## GENE NAME - ENSG CONVERSIONS

def getGeneNameDict():
    '''Make dictionary of IDs to names'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    geneIdNameDict = dict([(ggg[0], ggg[1]) for idx, ggg in gene_info.iterrows()])
    return geneIdNameDict
    
def ccd_gene_names(id_list_like, geneIdNameDict):
    '''Convert gene ID list to gene name list'''
    return np.unique([geneIdNameDict[ggg] for ggg in id_list_like if ggg in geneIdNameDict])

def ccd_gene_names_gapped(id_list_like, geneIdNameDict):
    '''Convert gene ID list to gene name list'''
    return [geneIdNameDict[idd] if idd in geneIdNameDict else "" for idd in id_list_like]

def getHgncDict():
    '''Make dictionary of IDs to HGNC symbols'''
    geneIdToHgncTable = pd.read_csv("input/processed/python/ENSGToHGNC.csv", index_col=False, header=0)
    geneIdToHgncDict = dict([(ggg[1], ggg[0]) for idx, ggg in geneIdToHgncTable.iterrows()])
    return geneIdToHgncDict

def geneIdToHngc(id_list_like, geneDict):
    '''Convert gene ID list to HNGC symbol if it exists'''
    return np.unique([geneDict[ggg] for ggg in id_list_like if ggg])

def geneIdToHngc_withgaps(id_list_like, geneDict):
    '''Convert gene ID list to HNGC symbol if it exists'''
    return [geneDict[ggg] for ggg in id_list_like]

def ccd_gene_lists(adata):
    '''Read in the published CCD genes / Diana's CCD / Non-CCD genes'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    ccd_regev=pd.read_csv("input/processed/manual/ccd_regev.txt")   
    wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
    ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
    nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)
    ccd=wp_ensg[ccd_comp]
    nonccd=wp_ensg[nonccd_comp]
    ccd_regev_filtered = list(gene_info[(gene_info["name"].isin(ccd_regev["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    ccd_filtered = list(ccd[np.isin(ccd, adata.var_names)])
    nonccd_filtered = list(nonccd[np.isin(nonccd, adata.var_names)])
    return ccd_regev_filtered, ccd_filtered, nonccd_filtered

def save_category(genelist, filename):
    pd.DataFrame({"gene" : genelist}).to_csv(filename, index=False, header=False)

def save_gene_names_by_category(adata, wp_ensg, ccd_comp, nonccd_comp, ccdtranscript):
    '''Save files containing the gene names for each category of CCD proteins/transcripts'''
    ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)
    genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])
    bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
    knownccd1 = np.genfromtxt("input/processed/manual/knownccd.txt", dtype='str') # from gene ontology, reactome, cyclebase 3.0, NCBI gene from mcm3
    knownccd2 = np.genfromtxt("input/processed/manual/known_go_ccd.txt", dtype='str') # from GO cell cycle
    knownccd3 = np.genfromtxt("input/processed/manual/known_go_proliferation.txt", dtype='str') # from GO proliferation
    knownccd = np.concatenate((knownccd1, knownccd2, knownccd3))

    # Get the ENSG symbols for lists for GO analysis
    ensg_ccdtranscript = np.unique(adata.var_names[ccdtranscript])
    ensg_nonccdtranscript = np.unique(adata.var_names[~ccdtranscript])
    ensg_ccdprotein = np.unique(np.concatenate((wp_ensg[ccd_comp], bioccd)))
    ensg_nonccdprotein = np.unique(wp_ensg[nonccd_comp & ~np.isin(wp_ensg, bioccd)])
    ensg_ccdprotein_treg = np.unique(ensg_ccdprotein[np.isin(ensg_ccdprotein, ensg_ccdtranscript)])
    ensg_ccdprotein_nontreg = np.unique(ensg_ccdprotein[~np.isin(ensg_ccdprotein, ensg_ccdtranscript)])
    ensg_knownccdprotein = ensg_ccdprotein[np.isin(ensg_ccdprotein, knownccd)]
    ensg_novelccdprotein = ensg_ccdprotein[~np.isin(ensg_ccdprotein, knownccd)]
    
    # Get the HGNC symbols for lists for GO analysis
    geneIdToHgncDict = getHgncDict()
    hgnc_ccdtranscript = geneIdToHngc(ensg_ccdtranscript, geneIdToHgncDict)
    hgnc_ccdprotein_transcript_regulated = geneIdToHngc(ensg_ccdprotein_treg, geneIdToHgncDict)
    hgnc_ccdprotein_nontranscript_regulated = geneIdToHngc(ensg_ccdprotein_nontreg, geneIdToHgncDict)
    hgnc_nonccdprotein = geneIdToHngc(ensg_nonccdprotein, geneIdToHgncDict)
    hgnc_ccdprotein = geneIdToHngc(ensg_ccdprotein, geneIdToHgncDict)
    
    # Convert to gene names and store them as such
    geneIdNameDict = getGeneNameDict()
    names_ccdtranscript = ccd_gene_names(ensg_ccdtranscript, geneIdNameDict)
    names_nonccdtranscript = ccd_gene_names(ensg_nonccdtranscript, geneIdNameDict)
    names_ccdprotein = ccd_gene_names(ensg_ccdprotein, geneIdNameDict)
    names_nonccdprotein = ccd_gene_names(ensg_nonccdprotein, geneIdNameDict)
    names_ccdprotein_transcript_regulated = ccd_gene_names(ensg_ccdprotein_treg, geneIdNameDict)
    names_ccdprotein_nontranscript_regulated = ccd_gene_names(ensg_ccdprotein_nontreg, geneIdNameDict)
    names_genes_analyzed = ccd_gene_names(genes_analyzed, geneIdNameDict)
    names_ccd_regev_filtered = ccd_gene_names(ccd_regev_filtered, geneIdNameDict)
    names_ccd_filtered = ccd_gene_names(ccd_filtered, geneIdNameDict)
    
    # Save the HGNC gene names for each category
    save_category(hgnc_ccdtranscript, "output/hgnc_ccdtranscript.csv")
    save_category(hgnc_ccdprotein_transcript_regulated, "output/hgnc_ccdprotein_transcript_regulated.csv")
    save_category(hgnc_ccdprotein_nontranscript_regulated, "output/hgnc_ccdprotein_nontranscript_regulated.csv")
    save_category(hgnc_nonccdprotein, "output/hgnc_nonccdprotein.csv")
    save_category(hgnc_ccdprotein, "output/hgnc_ccdprotein.csv")

    # Save the geneIds for each category
    save_category(ensg_ccdtranscript, "output/ensg_ccdtranscript.csv")
    save_category(ensg_nonccdtranscript, "output/ensg_nonccdtranscript.csv")
    save_category(ensg_ccdprotein_treg, "output/ensg_ccdprotein_transcript_regulated.csv")
    save_category(ensg_ccdprotein_nontreg, "output/ensg_ccdprotein_nontranscript_regulated.csv")
    save_category(ensg_nonccdprotein, "output/ensg_nonccdprotein.csv")
    save_category(ensg_ccdprotein, "output/ensg_ccdprotein.csv")
    save_category(ensg_knownccdprotein, "output/ensg_knownccdprotein.csv")
    save_category(ensg_novelccdprotein, "output/ensg_novelccdprotein.csv")

    # Save the gene names for each category
    save_category(names_ccdtranscript, "output/names_ccdtranscript.csv")
    save_category(names_nonccdtranscript, "output/names_nonccdtranscript.csv")
    save_category(names_ccdprotein, "output/names_ccdprotein.csv")
    save_category(names_nonccdprotein, "output/names_nonccdprotein.csv")
    save_category(names_ccdprotein_transcript_regulated, "output/names_ccdprotein_transcript_regulated.csv")
    save_category(names_ccdprotein_nontranscript_regulated, "output/names_ccdprotein_nontranscript_regulated.csv")
    save_category(names_genes_analyzed, "output/names_genes_analyzed.csv")
    save_category(names_ccd_regev_filtered, "output/names_ccd_regev_filtered.csv")
    save_category(names_genes_analyzed, "output/names_genes_analyzed.csv")
    save_category(names_ccd_filtered, "output/names_ccd_filtered.csv")
    
    return ((ensg_ccdtranscript, ensg_nonccdtranscript, ensg_ccdprotein, ensg_nonccdprotein, 
            ensg_ccdprotein_treg, ensg_ccdprotein_nontreg, 
            genes_analyzed, ccd_regev_filtered, ccd_filtered),
        (names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, 
            names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, 
            names_genes_analyzed, names_ccd_regev_filtered, names_ccd_filtered))
    