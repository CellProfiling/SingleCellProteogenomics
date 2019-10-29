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
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists

#%% Read data into scanpy; Read phases and FACS intensities
# Use RPKMs for comparing genes
dd = "All"
counts_or_rpkms = "Tpms"
use_spike_ins = False
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, use_spike_ins)

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


#%% UMAP with just the Regev cell-cycle dependent genes (CCD) 
# Idea: Do the UMAPs from before look similar to UMAPs using only published CCD genes 
# (https://science.sciencemag.org/content/352/6282/189)?
# Execution: scanpy
# Output: UMAPs
adata_ccdregev = adata[:, ccd_regev_filtered]
adata_ccdregev.var_names_make_unique()
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="phase", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPhaseCcdRegev.pdf")

#%% There is a nodule in the UMAP; are there any highly expressing genes in there?
# Execution: subset the cells in that corner of the UMAP
# Output: Genes that are differentially expressed in this quadrant compared to the rest in G1

#%% Expression UMAP, uncomment to run again
# Idea: For each gene in the lists, overlay expression on the UMAP
# Execution: scanpy, using the non-log normalized counts because we're not comparing genes
# Output: lots of UMAPs

dd = "All"
counts_or_rpkms = "Tpms"
do_log_normalization = True
adata, phases_filt = read_counts_and_phases(dd, counts_or_rpkms, use_spike_ins)
qc_filtering(adata, do_log_normalization)

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

#%%
