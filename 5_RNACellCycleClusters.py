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

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
counts_or_rpkms = "Tpms"
do_log_normalization = True
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, False, "protein_coding")
phases_filt = qc_filtering(adata, do_log_normalization)
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

#%% Make the peak RNA heatmap
# Idea: generate a heatmap similar to the one for protein data
# Execution: use the list of significantly differentially expressed genes to pluck the normalized intensities for those genes
#      and order them by the time of peak expression
# output: heatmap

def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')

# imports
dd = "All"
counts_or_rpkms = "Counts"
do_log_normalization = True
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, False, "protein_coding")
expression_data = adata.X # log normalized
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T # divide for cell
# normalized_exp_data = (normalized_exp_data / np.max(normalized_exp_data, axis=1)[:,None]) # divide for gene
fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)

