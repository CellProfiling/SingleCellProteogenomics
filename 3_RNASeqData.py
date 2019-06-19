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

#%% Imports
from imports import *

#%% Read data into scanpy; Read phases and FACS intensities
# Use RPKMs for comparing genes
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
counts_or_rpkms = "Rpkms"
adata, phases_filt = read_counts_and_phases(dd, counts_or_rpkms)

#%% QC and filtering
do_log_normalization = True
sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=True)
shutil.move("figures/highest_expr_genes.pdf", f"figures/highest_expr_genes_{dd}Cells.pdf")
qc_filtering(adata, do_log_normalization)

# Post filtering QC
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=True, save=True)
shutil.move("figures/filter_genes_dispersion.pdf", f"figures/filter_genes_dispersion{dd}Cells.pdf")

# Read in the published CCD genes / Diana's CCD / Non-CCD genes
# filter for genes that weren't filtered in QC
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% UMAP plots
# Idea: based on the expression of all genes, do the cell cycle phases cluster together?
# Execution: scanpy methods
# Output: UMAP plots

# UMAP statistics first
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.umap(adata, color=["phase"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqCenterPhase.pdf")
sc.pl.umap(adata, color = ["phase_ajc"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhase.pdf")


#%% UMAP with just the Regev cell-cycle dependent genes (CCD) 
# Idea: Do the UMAPs from before look similar to UMAPs using only published CCD genes 
# (https://science.sciencemag.org/content/352/6282/189)?
# Execution: scanpy
# Output: UMAPs
adata_ccdregev = adata[:, ccd_regev_filtered]
adata_ccdregev.var_names_make_unique()
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="phase_ajc", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhaseCcdRegev.pdf")

#%% There is a nodule in the UMAP; are there any highly expressing genes in there?
# Execution: subset the cells in that corner of the UMAP
# Output: Genes that are differentially expressed in this quadrant compared to the rest in G1

#%% Expression UMAP, uncomment to run again
# Idea: For each gene in the lists, overlay expression on the UMAP
# Execution: scanpy, using the non-log normalized counts because we're not comparing genes
# Output: lots of UMAPs

dd = "All"
counts_or_rpkms = "Counts"
do_log_normalization = False
adata, phases_filt = read_counts_and_phases(dd, counts_or_rpkms)
qc_filtering(adata, do_log_normalization)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

# UMAP statistics first
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

def plot_expression_umap(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        adata.obs[gene] = normalized_exp_data
        sc.pl.umap(adata, color=gene, show=False, save=True)
        shutil.move("figures/umap.pdf", f"{outfolder}/{gene}.pdf")
        plt.close()
        adata.obs.drop(gene, 1)

# plot_expression_umap(ccd_regev_filtered, "figures/RegevGeneExpressionUmap")
# plot_expression_umap(ccd_filtered, "figures/DianaCcdGeneExpressionUmap")
# plot_expression_umap(nonccd_filtered, "figures/DianaNonCcdGeneExpressionUmap")

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

