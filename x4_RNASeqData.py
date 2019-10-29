# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:19:25 2019

@author: antho
"""

#%% [markdown]
# # Summary: Cell Cycle scRNA-Seq Analysis
# We've shown in single cell imaging data that there is variability that is correlated to the cell cycle, as well as a majority of proteins that vary outside of the cell cycle, which might be due to metabolic processes or other sources of variation.
# 
# Here, we've collected single-cell RNA-Seq (scRNA-Seq) data for these cells.
# * How many cells were analyzed?
# * How many reads per cell?
# * How many genes show variability in expression at the RNA level?
# * How many of the genes that are temporally regulated over the cell cycle, using the fucci colors to build the cell cycle trajectory?
# * How many of the genes that show variability not correlated to the cell cycle?

#%% 

#%% Heatmaps with published CCD genes again
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, ccd_regev_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdRegev.pdf")
# sc.pl.heatmap(adata, ccd_regev_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdRegev.pdf")

#%% Heatmaps with Diana's first few CCD genes
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, ccd_filtered[:50], "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectCcdDiana.pdf")
# sc.pl.heatmap(adata, ccd_filtered[:50], "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectCcdDiana.pdf")

# # Heatmaps with Diana's CCD genes
# sc.pl.heatmap(adata, ccd_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdDiana.pdf")
# sc.pl.heatmap(adata, ccd_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdDiana.pdf")

#%% Heatmaps with Diana's first few Non-CCD genes
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, nonccd_filtered[:50], "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectNonCcdDiana.pdf")
# sc.pl.heatmap(adata, nonccd_filtered[:50], "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectNonCcdDiana.pdf")

# # Heatmaps with Diana's Non-CCD genes
# sc.pl.heatmap(adata, nonccd_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseNonCcdDiana.pdf")
# sc.pl.heatmap(adata, nonccd_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseNonCcdDiana.pdf")

#%% Louvian clustering -- this doesn't work; they use C++ packages that don't install easily.
# I checked and I can use Seurat for this in R
# sc.tl.louvain(adata)

#%% partition-based graph abstraction (PAGA), and this doesn't work because of louvain either
# sc.tl.paga(adata)
# sc.pl.paga_compare(adata, edges=True, threshold=0.05)

#%% pseudotime reconstruction (FUCCI is better; this one is kind of arbitrary)
# adata.uns["iroot"] = 0
# sc.tl.dpt(adata, n_branchings=0)
# sc.pl.umap(adata, color='dpt_pseudotime', show=True, save=True)
# shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPredictedPseudotime.pdf")
# sc.pl.diffmap(adata, color='dpt_pseudotime', projection="3d", show=True, save=True)
# shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsPredictedPseudotime3d.pdf")