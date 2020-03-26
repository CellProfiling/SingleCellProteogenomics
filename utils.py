import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcolors
import matplotlib.patches as mpatches
import scanpy as sc
import os
import shutil
import scipy
import scipy.stats

def read_counts_and_phases(dd, count_or_rpkm, use_spike_ins, biotype_to_use):
    '''
    Read data into scanpy; Read phases and FACS intensities
    * dd: "All", "355", "356", "357"
    * count_or_rpkm: "Counts" or "Tpms"
    '''
    read_file = f"input/processed/scanpy/{count_or_rpkm}.csv" + (".ercc.csv" if use_spike_ins else "")
    if biotype_to_use != None and len(biotype_to_use) > 0:
        print(f"filtering for biotype: {biotype_to_use}")
        biotype_file = f"{read_file}.{biotype_to_use}.csv"
        if not os.path.exists(biotype_file):
            gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
            biotyped = gene_info[gene_info["biotype"] == biotype_to_use]["gene_id"]
            pd.read_csv(read_file)[biotyped ].to_csv(biotype_file, index=False)
        read_file = biotype_file

    adata = sc.read_csv(read_file)
    print(f"data shape: {adata.X.shape}")
    # adata.raw = adata

    phases = pd.read_csv("input/processed/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")
    
    # Assign phases and log intensities; require log intensity
    adata.obs["phase"] = np.array(phases["Stage"])
    adata.obs["Green530"] = np.array(phases["Green530"])
    adata.obs["Red585"] = np.array(phases["Red585"])
    adata = adata[pd.notnull(adata.obs["Green530"]) & pd.notnull(adata.obs["Red585"])] # removes dark mitotic cells
    
    # Read in fucci pseudotime from previous analysis
    if os.path.isfile("output/fucci_time.csv"):
        adata.obs["fucci_time"] = np.array(pd.read_csv("output/fucci_time.csv")["fucci_time"])

    # Get info about the genes
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", header=None, names=["name", "biotype", "description"], index_col=0)
    adata.var["name"] = gene_info["name"]
    adata.var["biotype"] = gene_info["biotype"]
    adata.var["description"] = gene_info["description"]

    return adata, phases


def qc_filtering(adata, do_log_normalize, do_remove_blob):
    '''QC and filtering; remove cluster of cells in senescent G0'''
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=100)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    if do_log_normalize: sc.pp.log1p(adata)
    louvain = np.load("input/processed/python/louvain.npy", allow_pickle=True)
    adata.obs["louvain"] = louvain
    if do_remove_blob: adata = adata[adata.obs["louvain"] != "5",:]
    phases = pd.read_csv("input/processed/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")
    phases_filt = phases[phases["Well_Plate"].isin(adata.obs_names)]
    phases_filt = phases_filt.reset_index(drop=True) # remove filtered indices
    print(f"data shape after filtering: {adata.X.shape}")
    return adata, phases_filt


def ccd_gene_lists(adata):
    '''Read in the published CCD genes / Diana's CCD / Non-CCD genes'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    ccd_regev=pd.read_csv("input/processed/manual/ccd_regev.txt")    
    ccd=pd.read_csv("output/picklestxt/ccd_compartment_ensg.txt")#"input/processed/manual/ccd_genes.txt")
    nonccd=pd.read_csv("output/picklestxt/nonccd_compartment_ensg.txt")#("input/processed/manual/nonccd_genes.txt")
    ccd_regev_filtered = list(gene_info[(gene_info["name"].isin(ccd_regev["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    ccd_filtered = list(gene_info[(gene_info["name"].isin(ccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    nonccd_filtered = list(gene_info[(gene_info["name"].isin(nonccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    return ccd_regev_filtered, ccd_filtered, nonccd_filtered

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

def np_save_overwriting(fn, arr):
    '''Helper function to always overwrite numpy pickles'''
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)