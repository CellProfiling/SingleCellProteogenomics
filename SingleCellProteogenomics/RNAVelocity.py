#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 15:26:58 2021

@author: anthony.cesnik
"""

import pandas as pd
import numpy as np
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable
plt.rcParams['figure.figsize'] = (10, 10)

def analyze_rna_velocity(adata_withblob, mean_diff_from_rng, do_all_plotting):
    adata_blobless = adata_withblob[adata_withblob.obs["original_louvain"] != "5",:]
    
    # Make the chosen cutoff UMAP
    cutoff=0.08
    adata_withCutoff = adata_blobless[:,mean_diff_from_rng > cutoff]
    scv.pp.moments(adata_withCutoff, n_neighbors=20)
    scv.tl.velocity(adata_withCutoff, mode='stochastic')
    scv.tl.velocity_graph(adata_withCutoff)
    sc.tl.umap(adata_withCutoff)
    sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', show=False, save=f"{cutoff}CCD.pdf")
    scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', show=False, save=f"ScveloStreamUmap{cutoff}CCD.png", color_map="viridis", smooth=2, density=0.6)
    
    pd.DataFrame({
        "cell_well_plate" : adata_withCutoff.obs_names,
        "cell_cycle_hrs" : adata_withCutoff.obs["fucci_time"],
        "x_umap" : [x[0] for x in adata_withCutoff.obsm['velocity_umap']],
        "y_umap" : [x[1] for x in adata_withCutoff.obsm['velocity_umap']]}).to_csv(
            f"figures/velocityStreamUmap{cutoff}CCD.tsv", index=False, sep="\t")
    
    if not do_all_plotting:
        return
            
    # Profile different cutoffs for CCD and non-CCD, UMAP and tSNE to evaluate params and consistency
    for cutoff in (np.arange(20) + 1) / 100:
        adata_withCutoff = adata_blobless[:,mean_diff_from_rng <= cutoff]
        scv.pp.moments(adata_withCutoff, n_neighbors=20)
        scv.tl.velocity(adata_withCutoff, mode='stochastic')
        scv.tl.velocity_graph(adata_withCutoff)
        sc.tl.umap(adata_withCutoff)
        sc.tl.tsne(adata_withCutoff, n_pcs=2)
        sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', show=False, save=f"{cutoff}NonCCD.pdf")
        sc.pl.tsne(adata_withCutoff, color="fucci_time", cmap='viridis', show=False, save=f"{cutoff}NonCCD.pdf")
        scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', show=False, save=f"Umap{cutoff}NonCCD.png", color_map="viridis", smooth=2, density=0.6)
        scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='tsne', show=False, save=f"Tsne{cutoff}NonCCD.png", color_map="viridis", smooth=2, density=0.6)
        
        adata_withCutoff = adata_blobless[:,mean_diff_from_rng > cutoff]
        scv.pp.moments(adata_withCutoff, n_neighbors=20)
        scv.tl.velocity(adata_withCutoff, mode='stochastic')
        scv.tl.velocity_graph(adata_withCutoff)
        sc.tl.umap(adata_withCutoff)
        sc.tl.tsne(adata_withCutoff, n_pcs=2)
        sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', show=False, save=f"{cutoff}CCD.pdf")
        sc.pl.tsne(adata_withCutoff, color="fucci_time", cmap='viridis', show=False, save=f"{cutoff}CCD.pdf")
        scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', show=False, save=f"Umap{cutoff}CCD.png", color_map="viridis", smooth=2, density=0.6)
        scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='tsne', show=False, save=f"Tsne{cutoff}CCD.png", color_map="viridis", smooth=2, density=0.6)
        
    scv.pp.moments(adata_blobless, n_neighbors=20)
    scv.tl.velocity(adata_blobless, mode="stochastic")
    scv.tl.velocity_graph(adata_blobless)
    sc.tl.umap(adata_blobless)
    sc.tl.tsne(adata_blobless, n_pcs=2)
    scv.pl.velocity_embedding_stream(adata_blobless, color="fucci_time", basis='umap', show=False, save="UmapAll.png", color_map="viridis", smooth=2, density=0.6)
    scv.pl.velocity_embedding_stream(adata_blobless, color="fucci_time", basis='tsne', show=False, save="TsneAll.png", color_map="viridis", smooth=2, density=0.6)