# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:02:05 2020

Velocity calculations. Use the conda environment in velo.yaml.

Requires running the scripts 1_ through 3_ before running.

@author: antho
"""

import os
import pandas as pd
import numpy as np
import scvelo as scv
import scanpy as sc
import shutil
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable
plt.rcParams['figure.figsize'] = (10, 10)

adata = sc.read_csv("input/RNAData/Tpms.csv.protein_coding.csv")
adata.obs_names = pd.read_csv("input/RNAData/Tpms.obs_names.csv")["well_plate"]
phases = pd.read_csv("input/ProteinData/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")

# Assign phases and log intensities; require log intensity
adata.obs["phase"] = np.array(phases["Stage"])
adata.obs["Green530"] = np.array(phases["Green530"])
adata.obs["Red585"] = np.array(phases["Red585"])
adata = adata[pd.notnull(adata.obs["Green530"]) & pd.notnull(adata.obs["Red585"])] # removes dark mitotic cells
adata.obs["fucci_time"] = np.array(pd.read_csv("output/fucci_time.csv")["fucci_time"])

# Get info about the genes
gene_info = pd.read_csv("input/RNAData/IdsToNames.csv.gz", header=None, names=["name", "biotype", "description"], index_col=0)
adata.var["name"] = gene_info["name"]
adata.var["biotype"] = gene_info["biotype"]
adata.var["description"] = gene_info["description"]

ldata = scv.read("input/RNAData/a.loom", cache=True)
ldata.obs_names = pd.read_csv("input/RNAData/a.obs_names.csv")["well_plate"]
ldata.var["GeneName"] = ldata.var_names
ldata.var_names = ldata.var["Accession"]
adata = scv.utils.merge(adata, ldata, copy=True)

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=100)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

louvain = np.load("input/RNAData/louvain.npy", allow_pickle=True)
mean_diff_from_rng = np.load("output/pickles/mean_diff_from_rng.npy", allow_pickle=True)
adata.obs["louvain"] = louvain
adata_blobless = adata[adata.obs["louvain"] != "5",:]

# Make the chosen cutoff UMAP
cutoff=0.08
adata_withCutoff = adata_blobless[:,mean_diff_from_rng > cutoff]
scv.pp.moments(adata_withCutoff, n_neighbors=20)
scv.tl.velocity(adata_withCutoff, mode='stochastic')
scv.tl.velocity_graph(adata_withCutoff)
sc.tl.umap(adata_withCutoff)
sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', save=True)
shutil.move("figures/umap.pdf", f"figures/umap{cutoff}CCD.pdf")
scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', save="embedding_stream_umap.png", color_map="viridis", smooth=2, density=0.6)
shutil.move(f"figures/scvelo_embedding_stream_umap.png", f"figures/velocityStreamUmap{cutoff}CCD.png")

pd.DataFrame({
    "cell_well_plate" : adata_withCutoff.obs_names,
    "cell_cycle_hrs" : adata_withCutoff.obs["fucci_time"],
    f"x_umap" : [x[0] for x in adata_withCutoff.obsm['velocity_umap']],
    f"y_umap" : [x[1] for x in adata_withCutoff.obsm['velocity_umap']]}).to_csv(
        f"figures/velocityStreamUmap{cutoff}CCD.tsv", index=False, sep="\t")

# Profile different cutoffs for CCD and non-CCD, UMAP and tSNE to evaluate params and consistency
for cutoff in (np.arange(20) + 1) / 100:
    adata_withCutoff = adata_blobless[:,mean_diff_from_rng <= cutoff]
    scv.pp.moments(adata_withCutoff, n_neighbors=20)
    scv.tl.velocity(adata_withCutoff, mode='stochastic')
    scv.tl.velocity_graph(adata_withCutoff)
    sc.tl.umap(adata_withCutoff)
    sc.tl.tsne(adata_withCutoff, n_pcs=2)
    sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', save=True)
    sc.pl.tsne(adata_withCutoff, color="fucci_time", cmap='viridis', save=True)
    shutil.move("figures/umap.pdf", f"figures/umap{cutoff}NonCCD.pdf")
    shutil.move("figures/tsne.pdf", f"figures/tsne{cutoff}NonCCD.pdf")
    scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', save="embedding_stream_umap.png", color_map="viridis", smooth=2, density=0.6)
    scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='tsne', save="embedding_stream_tsne.png", color_map="viridis", smooth=2, density=0.6)
    shutil.move(f"figures/scvelo_embedding_stream_umap.png", f"figures/velocityStreamUmap{cutoff}NonCCD.png")
    shutil.move(f"figures/scvelo_embedding_stream_tsne.png", f"figures/velocityStreamTsne{cutoff}NonCCD.png")
    
    adata_withCutoff = adata_blobless[:,mean_diff_from_rng > cutoff]
    scv.pp.moments(adata_withCutoff, n_neighbors=20)
    scv.tl.velocity(adata_withCutoff, mode='stochastic')
    scv.tl.velocity_graph(adata_withCutoff)
    sc.tl.umap(adata_withCutoff)
    sc.tl.tsne(adata_withCutoff, n_pcs=2)
    sc.pl.umap(adata_withCutoff, color="fucci_time", cmap='viridis', save=True)
    sc.pl.tsne(adata_withCutoff, color="fucci_time", cmap='viridis', save=True)
    shutil.move("figures/umap.pdf", f"figures/umap{cutoff}CCD.pdf")
    shutil.move("figures/tsne.pdf", f"figures/tsne{cutoff}CCD.pdf")
    scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='umap', save="embedding_stream_umap.png", color_map="viridis", smooth=2, density=0.6)
    scv.pl.velocity_embedding_stream(adata_withCutoff, color="fucci_time", basis='tsne', save="embedding_stream_tsne.png", color_map="viridis", smooth=2, density=0.6)
    shutil.move(f"figures/scvelo_embedding_stream_umap.png", f"figures/velocityStreamUmap{cutoff}CCD.png")
    shutil.move(f"figures/scvelo_embedding_stream_tsne.png", f"figures/velocityStreamTsne{cutoff}CCD.png")
    
scv.pp.moments(adata_blobless, n_neighbors=20)
scv.tl.velocity(adata_blobless, mode="stochastic")
scv.tl.velocity_graph(adata_blobless)
sc.tl.umap(adata_blobless)
sc.tl.tsne(adata_blobless, n_pcs=2)
scv.pl.velocity_embedding_stream(adata_blobless, color="fucci_time", basis='umap', save="embedding_stream_umap.png", color_map="viridis", smooth=2, density=0.6)
scv.pl.velocity_embedding_stream(adata_blobless, color="fucci_time", basis='tsne', save="embedding_stream_tsne.png", color_map="viridis", smooth=2, density=0.6)
shutil.move(f"figures/scvelo_embedding_stream_umap.png", f"figures/velocityStreamUmapAll.png")
shutil.move(f"figures/scvelo_embedding_stream_tsne.png", f"figures/velocityStreamTsneAll.png")