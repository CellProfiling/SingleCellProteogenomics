# -*- coding: utf-8 -*-
"""
Preparation of single-cell RNA sequencing data:
    - Loads FACS intensity data for FUCCI markers for each cell
    - Uses scanpy for loading, filtering, and making QC plots of RNA-Seq data across single cells
    - Evaluates global patterns of RNA expression with UMAP dimensionality reduction
    - Evaluates increase in reads and gene count over the cell cycle with increasing cell size

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import utils
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import seaborn as sbn
import scvelo as scv

def read_counts_and_phases(count_or_rpkm, use_spike_ins, biotype_to_use, u_plates, use_isoforms=False, load_velocities=False):
    '''
    Read data into scanpy; Read phases and FACS intensities
        - count_or_rpkm: Must be "Counts" or "Tpms"
    '''
    read_file = f"input/RNAData/{count_or_rpkm}{'_Isoforms' if use_isoforms else ''}.csv" + (".ercc.csv" if use_spike_ins else "")
    if biotype_to_use != None and len(biotype_to_use) > 0:
        print(f"filtering for biotype: {biotype_to_use}")
        biotype_file = f"{read_file}.{biotype_to_use}.csv"
        if not os.path.exists(biotype_file):
            gene_info = pd.read_csv(f"input/RNAData/IdsToNames{'_Isoforms' if use_isoforms else ''}.csv.gz", 
                                    index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
            biotyped = gene_info[gene_info["biotype"] == biotype_to_use]["gene_id"]
            pd.read_csv(read_file)[biotyped].to_csv(biotype_file, index=False)
        read_file = biotype_file

    adata = sc.read_csv(read_file)
    print(f"data shape: {adata.X.shape}")
    if load_velocities:
        adata.obs_names = pd.read_csv("input/RNAData/Tpms.obs_names.csv")["well_plate"]

    intensities, phases = [],[]
    for plate in u_plates:
        file = f"input/RNAData/180911_Fucci_single cell seq_ss2-18-{plate}_index sort export.csv"
        plateIntensities = pd.read_csv(file, skiprows=2)
        newColumns = list(plateIntensities.columns)
        newColumns[5] = "MeanGreen530"
        newColumns[6] = "MeanRed585"
        plateIntensities.columns = newColumns
        plateIntensities["Plate"] = [plate] * len(plateIntensities)
        plateIntensities["Well_Plate"] = [f"{w}_{plate}" for w in plateIntensities["Well"]]
        intensitiesSubFrame = plateIntensities[plateIntensities["Population"] == "All Events"]
        if len(intensities) == 0: intensities = intensitiesSubFrame
        else: intensities = intensities.append(intensitiesSubFrame, ignore_index=True)
        isPhaseRow = ~plateIntensities["Population"].isin(["All Events", "Cells", "Singlets"])
        phasesSubFrame = plateIntensities[isPhaseRow & (plateIntensities["% Total"] == "100.00%")]
        if len(phases) == 0: phases = phasesSubFrame
        else: phases = phases.append(phasesSubFrame, ignore_index=True)
    wp_idx = list(phases.columns).index("Well_Plate")
    pop_idx = list(phases.columns).index("Population")
    phases_lookup = dict([(row[1][wp_idx], row[1][pop_idx]) for row in phases.iterrows()])
            
    # Assign phases and log intensities; require log intensity
    intensities = intensities.sort_values(by="Well_Plate")
    adata.obs["Well_Plate"] = np.array(intensities["Well_Plate"])
    adata.obs["plate"] = np.array(intensities["Plate"])
    adata.obs["phase"] = np.array([phases_lookup[wp] if wp in phases_lookup else "N/A" for wp in intensities["Well_Plate"]])
    adata.obs["MeanGreen530"] = np.array(intensities["MeanGreen530"])
    adata.obs["MeanRed585"] = np.array(intensities["MeanRed585"])
    adata = adata[pd.notnull(adata.obs["MeanGreen530"]) & pd.notnull(adata.obs["MeanRed585"])] # removes 6 dark likely mitotic cells
    
    # Read in fucci pseudotime from previous analysis
    if os.path.isfile("output/fucci_time.csv"):
        adata.obs["fucci_time"] = np.array(pd.read_csv("output/fucci_time.csv")["fucci_time"])

    # Get info about the genes
    gene_info = pd.read_csv(f"input/RNAData/IdsToNames{'_Isoforms' if use_isoforms else ''}.csv.gz", 
                            header=None, names=["name", "biotype", "description"], index_col=0)
    adata.var["name"] = gene_info["name"]
    adata.var["biotype"] = gene_info["biotype"]
    adata.var["description"] = gene_info["description"]
    
    if load_velocities:
        ldata = scv.read("input/RNAData/a.loom", cache=True)
        ldata.obs_names = pd.read_csv("input/RNAData/a.obs_names.csv")["well_plate"]
        ldata.var["GeneName"] = ldata.var_names
        ldata.var_names = ldata.var["Accession"]
        adata = scv.utils.merge(adata, ldata, copy=True)

    return adata, phases

def zero_center_fucci(adata):
    '''Zero center and recenter after filtering'''
    plate = np.array(adata.obs["plate"])
    u_plate = np.unique(plate)
    well_plate = adata.obs["Well_Plate"]
    green_fucci, red_fucci = adata.obs["MeanGreen530"], adata.obs["MeanRed585"]
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
    adata.obs["Green530"] = log_green_fucci_zeroc_rescale # LogNormRescaleGreen530
    adata.obs["Red585"] = log_red_fucci_zeroc_rescale # LogNormRescaleRed585
    return adata
    
def qc_filtering(adata, do_log_normalize, do_remove_blob):
    '''QC and filtering; remove cluster of cells in senescent G0 from original clustering'''
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=100)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    if do_log_normalize: sc.pp.log1p(adata)
    louvain_file="input/RNAData/louvain_original.npy"
    louvain = np.load(louvain_file, allow_pickle=True)
    adata.obs["original_louvain_wp"] = louvain[:,0]
    adata.obs["original_louvain"] = louvain[:,1]
    if do_remove_blob: 
        removeThese = np.isin(adata.obs["Well_Plate"], adata.obs["original_louvain_wp"][adata.obs["original_louvain"] == "5"])
        adata = adata[~removeThese,:]
    phases_filt = np.array(adata.obs["phase"]) # remove filtered indices
    print(f"data shape after filtering: {adata.X.shape}")
    return adata, phases_filt

def is_ccd(adata, wp_ensg, ccd_comp, nonccd_comp, bioccd, ccd_regev_filtered):
    '''Return whether the genes in RNA-Seq analysis are 1) CCD by present protein analysis 2) non-CCD by present protein analysis, 3) curated published CCD'''
    ccdprotein = np.isin(adata.var_names, np.concatenate((wp_ensg[ccd_comp], bioccd)))
    nonccdprotein = np.isin(adata.var_names, wp_ensg[nonccd_comp]) & ~np.isin(adata.var_names, bioccd)
    regevccdgenes = np.isin(adata.var_names, ccd_regev_filtered)
    return ccdprotein, nonccdprotein, regevccdgenes

def general_plots(u_plates):
    '''Make plots to illustrate the results of the scRNA-Seq analysis'''
    valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
    adata, phases = read_counts_and_phases(valuetype, use_spikeins, biotype_to_use, u_plates)

    # QC plots before filtering
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save="AllCells.pdf")

    # Post filtering QC
    do_log_normalization = True
    do_remove_blob = False
    adata, phasesfilt = qc_filtering(adata, do_log_normalization, do_remove_blob)
    adata = zero_center_fucci(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, show=False, save="AllCells.pdf")

    # Louvain clustering and UMAP plots
    # Idea: based on the expression of all genes, do the cell cycle phases cluster together?
    # Execution: scanpy methods: UMAP statistics first, then make UMAP
    # Output: UMAP plots
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.louvain(adata)
    utils.np_save_overwriting("output/pickles/louvain.npy", adata.obs["louvain"])
    sc.tl.umap(adata)
    plt.rcParams['figure.figsize'] = (10, 10)
    sc.pl.umap(adata, color=["phase"], show=False, save="AllCellsSeqCenterPhase.pdf")

    # General display of RNA abundances in TPMs
    sbn.displot(np.concatenate(adata.X), color="tab:orange")
    plt.xlabel("TPM")
    plt.ylabel("Density")
    plt.savefig("figures/rna_abundance_density.pdf")
    plt.close()
    
def plot_markers_vs_reads(adata):
    utils.general_scatter(adata.obs["fucci_time"], adata.X[:,list(adata.var_names).index("ENSG00000112312")], 
                          "Fucci Pseudotime", "GMNN log10 TPM", "figures/GMNN_timeVsReadsScatter.png", False)
    utils.general_scatter(adata.obs["Green530"], adata.X[:,list(adata.var_names).index("ENSG00000112312")], 
                          "GMNN FUCCI Marker Intensity", "GMNN log10 TPM", "figures/GMNN_markerVsReadsScatter.png", False)
    utils.general_scatter(adata.obs["fucci_time"], adata.X[:,list(adata.var_names).index("ENSG00000167513")], 
                          "Fucci Pseudotime", "GMNN log10 TPM", "figures/CDT1_timeVsReadsScatter.png", False)
    utils.general_scatter(adata.obs["Red585"], adata.X[:,list(adata.var_names).index("ENSG00000167513")], 
                          "CDT1 FUCCI Marker Intensity", "CDT1 log10 TPM", "figures/CDT1_markerVsReadsScatter.png", False)

def plot_pca_for_batch_effect_analysis(adata, suffix):
    '''Make PCA plots to show batch effects if they exits'''
    sc.tl.pca(adata)
    sc.pl.pca(adata, color="plate", palette=['b','tab:orange','g','grey'], show=False, save=f"ByPlate_{suffix}.pdf")
    sc.pl.pca(adata, color="phase", palette=['b','tab:orange','g','grey'], show=False, save=f"ByPhase_{suffix}.pdf")
    pca_var = np.var(adata.obsm['X_pca'], axis=0)
    pd.DataFrame({"pca_var" : pca_var / sum(pca_var)}).to_csv(f"output/pca{suffix}_var.txt", sep="\t")
    
def analyze_noncycling_cells(u_plates):
    '''The raw UMAP shows a group of cells that appear sequestered from cycling; investigate those'''
    valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
    do_log_normalization = True
    do_remove_blob = False
    adata, phases = read_counts_and_phases(valuetype, use_spikeins, biotype_to_use, u_plates)
    adata, phasesfilt = qc_filtering(adata, do_log_normalization, do_remove_blob)
    adata = zero_center_fucci(adata)

    # Unsupervised clustering and generate the gene list for the uncycling cells, aka the unknown blob
    nneighbors = [5, 10, 15, 30, 100] # used nn=10 in the paper
    mindists = [0, 0.01, 0.05, 0.1, 0.5, 1] # used 0.5 (default) in the paper
    for nn in nneighbors:
        for md in mindists:
            sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=40)
            sc.tl.umap(adata, min_dist=md)
            sc.pl.umap(adata, color="louvain", show=False, save=f"louvain_clusters_before_nn{nn}_md{md}.pdf")
    sc.tl.rank_genes_groups(adata, groupby="louvain")
    p_blob=[a[5] for a in adata.uns["rank_genes_groups"]["pvals_adj"]]
    p_blob_sig = np.array(p_blob) < 0.01
    ensg_blob_sig=np.array([a[5] for a in adata.uns["rank_genes_groups"]["names"]])[p_blob_sig]
    np.savetxt("output/blob_genes.csv", ensg_blob_sig, fmt="%s", delimiter=",")

    # Remove the blob
    do_remove_blob = True
    adata, phases = read_counts_and_phases(valuetype, use_spikeins, biotype_to_use, u_plates)
    adata, phasesfilt = qc_filtering(adata, do_log_normalization, do_remove_blob)
    adata = zero_center_fucci(adata)
    for nn in nneighbors:
        for md in mindists:
            sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=40)
            sc.tl.umap(adata, min_dist=md)
            sc.pl.umap(adata, color="louvain", show=False, save=f"louvain_clusters_after_nn{nn}_md{md}.pdf")

def demonstrate_umap_cycle_without_ccd(adata):
    '''
    Idea: We should expect that the UMAP does not recognize a cycle when cycling transcripts are removed.
    Removing the CCD genes from the 93 curated genes or the 300-or-so CCD proteins identified in this work lead to UMAPs that are not cyclical. 
    '''
    ccd_regev_filtered, ccd_filtered, nonccd_filtered = utils.ccd_gene_lists(adata)
    adata_ccdregev = adata[:, [x for x in adata.var_names if x not in ccd_regev_filtered]]
    sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_ccdregev)
    sc.pl.umap(adata_ccdregev, color="fucci_time", show=False, save="AllCellsPhaseNonCcdCurated.pdf")

    adata_ccdregev = adata[:, [x for x in adata.var_names if x in nonccd_filtered]]
    sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_ccdregev)
    sc.pl.umap(adata_ccdregev, color="fucci_time", show=False, save="AllCellsPhaseNonCcd.pdf")

def readcount_and_genecount_over_pseudotime(u_plates):
    '''
    To demonstrate why we normalize read counts per cell, these plots show the increase in read count over the cell cycle as the cell grows.
    We also show the resulting increase in the number of genes detected.
    '''
    valuetype, use_spikeins, biotype_to_use = "Counts", False, "protein_coding"
    adata, phases = read_counts_and_phases(valuetype, use_spikeins, biotype_to_use, u_plates)
    expression_data = adata.X
    fucci_time_inds = np.argsort(adata.obs["fucci_time"])
    fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
    exp_sort = np.take(expression_data, fucci_time_inds, axis=0)

    # Total counts per cell, moving average
    exp = exp_sort.sum(axis=1)
    df = pd.DataFrame({"fucci_time" : fucci_time_sort, "total_counts" : exp})
    bin_size = 100
    plt.figure(figsize=(10,10))
    plt.plot(df["fucci_time"], 
            df["total_counts"].rolling(bin_size).mean(), 
            color="blue", 
            label=f"Moving Average by {bin_size} Cells")
    plt.xlabel("Fucci Pseudotime",size=36)
    plt.ylabel("Total RNA-Seq Counts",size=36)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig("figures/TotalCountsPseudotime.png")
    plt.close()

    # Total genes detected per cell, moving average
    gene_ct = np.count_nonzero(exp_sort, axis=1)
    df = pd.DataFrame({"fucci_time" : fucci_time_sort, "total_genes" : gene_ct})
    plt.figure(figsize=(10,10))
    plt.plot(df["fucci_time"], 
            df["total_genes"].rolling(bin_size).mean(), 
            color="blue", 
            label=f"Moving Average by {bin_size} Cells")
    plt.xlabel("Fucci Pseudotime",size=36)
    plt.ylabel("Total Genes Detected By RNA-Seq",size=36)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig("figures/TotalGenesPseudotime.png")
    plt.close()