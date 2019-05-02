#%% [markdown]
# # Cell Cycle scRNA-Seq Analysis
# We've shown in single cell imaging data that there is variability that is correlated to the cell cycle, as well as a majority of proteins that vary outside of the cell cycle, which might be due to metabolic processes or other sources of variation.
# 
# Here, we've collected single-cell RNA-Seq (scRNA-Seq) data for these cells.
# * How many cells were analyzed?
# * How many reads per cell?
# * How many genes show variability in expression at the RNA level?
# * How many of the genes that are temporally regulated over the cell cycle, using the fucci colors to build the cell cycle trajectory?
# * How many of the genes that show variability not correlated to the cell cycle?

#%%
# to render slides: open JupyterLab, hit plus, select terminal
# jupyter nbconvert C:\Users\antho\Dropbox\CodeLundbergGroup\TissueAtlasAnalysis.ipynb --to slides --post serve

get_ipython().system('pip install scanpy')
get_ipython().system('pip install pandas')
# get_ipython().system('pip install sklearn')
get_ipython().system('pip install ggplot')
get_ipython().system('conda install -y -c vtraag python-igraph ')

# Render our plots inline
get_ipython().run_line_magic('matplotlib', 'inline')

# some setup
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcolors
import scanpy as sc
import os
import shutil

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts355 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_355\\counts.tab", delimiter="\t", index_col=0)
    counts356 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_356\\counts.tab", delimiter="\t", index_col=0)
    counts357 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_357\\counts.tab", delimiter="\t", index_col=0)

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts355.columns += "_355"
    counts356.columns += "_356"
    counts357.columns += "_357"
    counts355.columns

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts = pd.concat([counts355,counts356,counts357], axis=1, sort=False)
    counts.T.sort_index().to_csv("AllCountsForScanpy.csv")
    counts355.T.sort_index().to_csv("355CountsForScanpy.csv")
    counts356.T.sort_index().to_csv("356CountsForScanpy.csv")
    counts357.T.sort_index().to_csv("357CountsForScanpy.csv")

#%% Read data into scanpy
dd = "All"
adata = sc.read_csv(dd + "CountsForScanpy.csv")
adata.var_names_make_unique()

#%%
phases = pd.read_csv("WellPlatePhases.csv").sort_values(by="Well_Plate")
phases_filt = phases[phases["Well_Plate"].isin(adata.obs_names)]
phases_filt = phases_filt.reset_index(drop=True) # remove filtered indices

#%%
adata.obs["phase"] = np.array(phases_filt["Stage"])
adata.obs["phase_ajc"] = np.array(phases_filt["StageAJC"])

#%%
sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=True)
shutil.move("figures/highest_expr_genes.pdf", f"figures/highest_expr_genes_{dd}Cells.pdf")

#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=True, save=True)
shutil.move("figures/filter_genes_dispersion.pdf", f"figures/filter_genes_dispersion{dd}Cells.pdf")

#%% UMAP statistics
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

#%% UMAP plot
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.umap(adata, color=["phase"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqCenterPhase.pdf")
sc.pl.umap(adata, color = ["phase_ajc"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhase.pdf")

#%% Read in the published CCD genes
ccd_regev=pd.read_csv("ccd_regev.txt")
ccd_regev_filtered = [gene for gene in ccd_regev["gene"] if gene in adata.var_names]
adata_ccdregev = adata[:, ccd_regev_filtered]
adata_ccdregev.var_names_make_unique()

#%% UMAP with just the Regev cell-cycle dependent genes (CCD) 
# https://science.sciencemag.org/content/352/6282/189
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="phase_ajc", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhaseCcdRegev.pdf")

#%% Heatmaps with published CCD genes again
sc.pl.heatmap(adata, ccd_regev_filtered, "phase", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdRegev.pdf")
sc.pl.heatmap(adata, ccd_regev_filtered, "phase_ajc", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdRegev.pdf")

#%%
ccd=pd.read_csv("ccd_genes.txt")
nonccd=pd.read_csv("nonccd_genes.txt")
ccd_filtered = [gene for gene in ccd["gene"] if gene in adata.var_names]
nonccd_filtered = [gene for gene in nonccd["gene"] if gene in adata.var_names]

#%% Heatmaps with Diana's first few CCD genes
sc.pl.heatmap(adata, ccd_filtered[:50], "phase", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectCcdDiana.pdf")
sc.pl.heatmap(adata, ccd_filtered[:50], "phase_ajc", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectCcdDiana.pdf")

# Heatmaps with Diana's CCD genes
sc.pl.heatmap(adata, ccd_filtered, "phase", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdDiana.pdf")
sc.pl.heatmap(adata, ccd_filtered, "phase_ajc", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdDiana.pdf")

#%% Heatmaps with Diana's first few Non-CCD genes
sc.pl.heatmap(adata, nonccd_filtered[:50], "phase", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectNonCcdDiana.pdf")
sc.pl.heatmap(adata, nonccd_filtered[:50], "phase_ajc", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectNonCcdDiana.pdf")

# Heatmaps with Diana's Non-CCD genes
sc.pl.heatmap(adata, nonccd_filtered, "phase", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseNonCcdDiana.pdf")
sc.pl.heatmap(adata, nonccd_filtered, "phase_ajc", show=True, save=True)
shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseNonCcdDiana.pdf")

#%% Louvian clustering
sc.tl.louvain(adata)

#%%
