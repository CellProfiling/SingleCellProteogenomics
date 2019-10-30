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

#%% Senescent G0 investigation start; archiving this because the louvain clustering helped a bunch to investigate this. 
# There is a nodule in the UMAP; are there any highly expressing genes in there?
# Execution: subset the cells in that corner of the UMAP
# Output: Genes that are differentially expressed in this quadrant compared to the rest in G1

#% Subset based on phase and separate out the weird nodule
# Idea: There's a nodule in the data that is interesting; what makes it unique?
# Exec: Use rank gene groups
# Outp: Volcano plots and genes that are significant for the different groups (>2 FC and <1e-9 pval_adj )
nodulecells = adata.obsm.X_umap[:,0] < -4 
phasessss = np.array(phases_filt["Stage"])
phasessss[nodulecells] = "G1-not"
adata.obs["phase+"] = phasessss
sc.tl.rank_genes_groups(adata, "phase+", groups=["G1", "G1-not", "G2M", "S-ph"], n_genes=len(adata.X[0,:]))
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.diffmap(adata, color='phase+', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsFucciPhasePlus.pdf")
sc.pl.umap(adata, color=["phase+"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqFucciPhasePlus.pdf")
# adata.uns["rank_genes_groups"]


groups = ["G1", "G1-not", "G2M", "S-ph"]
sig_gene_lists = []
for idx, ggg in enumerate(groups):
    plt.scatter([x[idx] for x in adata.uns["rank_genes_groups"]["logfoldchanges"]], 
        -np.log10([x[idx] for x in adata.uns["rank_genes_groups"]["pvals_adj"]]))
    plt.title(ggg)
    plt.xlabel("Log Fold Changes")
    plt.ylabel("P-Value, BH Adj")
    plt.savefig(f"figures/RankGenes{ggg}.png")
    plt.close()

    highly_expressed = np.array([x[idx] for x in adata.uns["rank_genes_groups"]["logfoldchanges"]]) > 1
    sig = np.array([x[idx] for x in adata.uns["rank_genes_groups"]["pvals_adj"]]) < 1e-9
    sig_genes = [g[idx] for g in adata.uns["rank_genes_groups"]["names"][(highly_expressed) & (sig)]]
    sig_gene_foldchange = [g[idx] for g in adata.uns["rank_genes_groups"]["logfoldchanges"][(highly_expressed) & (sig)]]
    sig_gene_pvaladj = [g[idx] for g in adata.uns["rank_genes_groups"]["pvals_adj"][(highly_expressed) & (sig)]]
    np.savetxt(f"output/{ggg}SignificantGenes.tsv", np.column_stack((sig_genes, sig_gene_foldchange, sig_gene_pvaladj)), fmt="%s", delimiter="\t")
    sig_gene_lists.append(sig_genes)


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

for idx, ggg in enumerate(groups):
    plot_expression_umap(sig_gene_lists[idx], f"figures/{ggg}GeneExpressionUmap")
