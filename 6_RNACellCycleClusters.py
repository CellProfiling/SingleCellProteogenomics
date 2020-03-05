#%% [markdown]
# # Cell Cycle scRNA-Seq Analysis
# We've shown in single cell imaging data that there is variability that is correlated to the cell cycle, as well as a majority of proteins that vary outside of the cell cycle, 
# which might be due to metabolic processes or other sources of variation.
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
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

#%% Read RNA-Seq data and illustrate the data
# Idea: Read in RNA-Seq data and do filtering
# Execution: Read data into scanpy; Read phases and FACS intensities; Use TPMs for comparing genes
# Output: QC, filtering, and UMAP plots (phased, expression, sanity check with known genes, fucci pseudotime)
dd = "All"
counts_or_rpkms = "Tpms"
use_spike_ins = False
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, use_spike_ins, "protein_coding")

# QC plots before filtering
sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=True)
shutil.move("figures/highest_expr_genes.pdf", f"figures/highest_expr_genes_{dd}Cells.pdf")

# Post filtering QC
do_log_normalization = True
do_remove_blob = False
adata, phasesfilt = qc_filtering(adata, do_log_normalization, do_remove_blob)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=True, save=True)
shutil.move("figures/filter_genes_dispersion.pdf", f"figures/filter_genes_dispersion{dd}Cells.pdf")

# Unsupervised clustering and gene list for the blob
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color="louvain", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap_louvain_clusters_before.pdf")
sc.tl.rank_genes_groups(adata, groupby="louvain")
p_blob=[a[5] for a in adata.uns["rank_genes_groups"]["pvals_adj"]]
p_blob_sig = np.array(p_blob) < 0.01
ensg_blob_sig=np.array([a[5] for a in adata.uns["rank_genes_groups"]["names"]])[p_blob_sig]
np.savetxt("output/blob_genes.csv", ensg_blob_sig, fmt="%s", delimiter=",")

# Remove the blob
do_remove_blob = True
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, use_spike_ins, "protein_coding")
adata, phasesfilt = qc_filtering(adata, do_log_normalization, do_remove_blob)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color="louvain", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap_louvain_clusters_after.pdf")

# Read in the published CCD genes / Diana's CCD / Non-CCD genes
# filter for genes that weren't filtered in QC
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

# UMAP plots
# Idea: based on the expression of all genes, do the cell cycle phases cluster together?
# Execution: scanpy methods: UMAP statistics first, then make UMAP
# Output: UMAP plots
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.umap(adata, color=["phase"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqCenterPhase.pdf")

# Displaying relative expression on the UMAP (using unlogged expression values)
# Output: a folder of UMAPs for each gene in several sets
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

# UMAP with just the Regev cell-cycle dependent genes (CCD) 
# Idea: Do the UMAPs from before look similar to UMAPs using only published CCD genes (https://science.sciencemag.org/content/352/6282/189)?
# Execution: scanpy
# Output: UMAPs
adata_ccdregev = adata[:, ccd_regev_filtered]
adata_ccdregev.var_names_make_unique()
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="phase", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPhaseCcdRegev.pdf")

# fucci pseudotime
# Idea: display pseudotime on the UMAP created from the gene expression
sc.tl.diffmap(adata)
sc.pl.diffmap(adata, color='fucci_time', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsFucciPseudotime3d.pdf")
sc.pl.umap(adata, color=["fucci_time"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqFucciPseudotime.pdf")



#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
counts_or_rpkms = "Tpms"
do_log_normalization = True
use_spike_ins = False
do_remove_blob = True
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, use_spike_ins, "protein_coding")
adata, phases_filt = qc_filtering(adata, do_log_normalization, do_remove_blob)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)


#%% Expression vs Pseudotime, uncomment to run again
# Idea: plot the expression of genes against the pseudotime calculated using the fucci FACS intensities
# Exec: scatters (expression values; moving avg of expression values)
# Outp: scatters
expression_data = adata.X
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T

def plot_expression_pseudotime(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        nexp = normalized_exp_data[:,list(adata.var_names).index(gene)]
        plt.scatter(adata.obs["fucci_time"], nexp)
        plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
        plt.ylabel("RNA-Seq Expression, Normalized By Cell",size=36,fontname='Arial')
        plt.title(gene,size=36,fontname='Arial')
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

# plot_expression_pseudotime(ccd_regev_filtered, "figures/RegevGeneProfiles")
# plot_expression_pseudotime(ccd_filtered, "figures/DianaCcdGeneProfiles")
# plot_expression_pseudotime(nonccd_filtered, "figures/DianaNonCcdGeneProfiles")
# for idx, ggg in enumerate(groups):
#     plot_expression_pseudotime(sig_gene_lists[idx], f"figures/{ggg}GeneProfiles")

def plot_expression_facs(genelist, pppp, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        nexp = normalized_exp_data[:,list(adata.var_names).index(gene)]
        plt.scatter(pppp["Green530"], pppp["Red585"], nexp)
        plt.tight_layout()
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel("Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
        plt.title(gene,size=20,fontname='Arial')
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

phasesFilt = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585)] # stage may be null

# plot_expression_facs(ccd_regev_filtered, phasesFilt, "figures/RegevGeneFucci")
# plot_expression_facs(ccd_filtered, phasesFilt, "figures/DianaCcdGeneFucci")
# plot_expression_facs(nonccd_filtered, phasesFilt, "figures/DianaNonCcdGeneFucci")
# for idx, ggg in enumerate(groups):
#     plot_expression_facs(sig_gene_lists[idx], phasesFilt, f"figures/{ggg}GeneFucci")


# MOVING AVERAGE EXPRESSION PLOTS
# Moving avg Expression vs Pseudotime, uncomment to run again
def moving_average(a, n):
    '''A formula for the moving average of an array (a) with window size (n)'''
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
cell_time_sort = pd.DataFrame({"fucci_time" : fucci_time_sort, "cell" : np.take(adata.obs_names, fucci_time_inds)})
cell_time_sort.to_csv("output/CellPseudotimes.csv")

def plot_expression_avg_pseudotime(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        nexp = norm_exp_sort[:,list(adata.var_names).index(gene)]
        df = pd.DataFrame({"fucci_time" : fucci_time_sort, gene : nexp})
        # plt.scatter(df["fucci_time"], df[gene], label="Normalized Expression")
        bin_size = 100
        plt.plot(df["fucci_time"], 
            df[gene].rolling(bin_size).mean(), 
            color="blue", 
            label=f"Moving Average by {bin_size} Cells")
        plt.fill_between(df["fucci_time"], 
            df[gene].rolling(bin_size).quantile(0.10),
            df[gene].rolling(bin_size).quantile(0.90), 
            color="lightsteelblue", 
            label="10th & 90th Percentiles")
        # plt.plot(df["fucci_time"], color="orange", label="Normalized Expression, 10th Percentile")
        # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() + 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
        # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() - 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
        plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
        plt.ylabel("RNA-Seq Counts, Normalized By Cell",size=36,fontname='Arial')
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.title(gene,size=36,fontname='Arial')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

plot_expression_avg_pseudotime(["ENSG00000104833"], "figures/OtherGeneProfileAvgs")
# plot_expression_avg_pseudotime(ccd_regev_filtered, "figures/RegevGeneProfileAvgs")
# plot_expression_avg_pseudotime(ccd_filtered, "figures/DianaCcdGeneProfileAvgs")
# plot_expression_avg_pseudotime(nonccd_filtered, "figures/DianaNonCcdGeneProfileAvgs")
# for idx, ggg in enumerate(groups):
#     plot_expression_avg_pseudotime(sig_gene_lists[idx], f"figures/{ggg}GeneProfileAvgs")


#%%
# Idea: Are there any differences in raw counts and such from the different parts
#       of the cell cycle?
# Execution: Based on the phase, before and after filtering,
#       summarize the 1) total count of reads per cell,
#       2) total number of genes detected per cell
# Output: total counts per cell; total genes per cell before and after filtering
#       plot moving avg total counts per cell on pseudotime
#       plot moving avg toal # genes per cell on pseudotime

dd = "All"
counts_or_rpkms = "Counts"
do_log_normalization = False
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, False, "protein_coding")
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
plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
plt.ylabel("Total RNA-Seq Counts",size=36,fontname='Arial')
plt.xticks(size=14)
plt.yticks(size=14)
# plt.title("Total Counts",size=36,fontname='Arial')
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(f"figures/TotalCountsPseudotime.png")
plt.show()
plt.close()

# Total genes detected per cell, moving average
gene_ct = np.count_nonzero(exp_sort, axis=1)
df = pd.DataFrame({"fucci_time" : fucci_time_sort, "total_genes" : gene_ct})
plt.figure(figsize=(10,10))
plt.plot(df["fucci_time"], 
        df["total_genes"].rolling(bin_size).mean(), 
        color="blue", 
        label=f"Moving Average by {bin_size} Cells")
plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
plt.ylabel("Total Genes Detected By RNA-Seq",size=36,fontname='Arial')
plt.xticks(size=14)
plt.yticks(size=14)
# plt.title("Total Genes ",size=36,fontname='Arial')
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(f"figures/TotalGenesPseudotime.png")
plt.show()
plt.close()

#%% Idea: What does the UMAP look like without the transcript CCD genes?
# Execution: filter out the transcript CCD genes from the counts table
# Output: UMAP for the RNA-Seq data
adata_ccdregev = adata[:, [x for x in adata.var_names if x not in ccd_regev_filtered]]
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="fucci_time", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPhaseNotCcdRegev.pdf")

adata_ccdregev = adata[:, [x for x in adata.var_names if x not in ccd_filtered]]
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="fucci_time", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPhaseNotCcdDiana.pdf")

adata_ccdregev = adata[:, [x for x in adata.var_names if x not in np.concatenate((ccd_filtered, ccd_regev_filtered, nonccd_filtered))]]
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="fucci_time", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPhaseNotListed.pdf")

#%% [markdown]
# Conclusion:
# There is information about the cell cycle in the genes that aren't cell cycle dependent


