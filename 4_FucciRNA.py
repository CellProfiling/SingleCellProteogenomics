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

#%% Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
counts_or_rpkms = "Tpms"
do_log_normalization = True
adata, phases = read_counts_and_phases(dd, counts_or_rpkms, False, "protein_coding")
phases_filt = qc_filtering(adata, do_log_normalization)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% fucci pseudotime
# Idea: display pseudotime on the UMAP created from the gene expression
# Exec: Scanpy
# Outp: umap visualization
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.diffmap(adata)
umap_coords = sc.tl.umap(adata)

dd = "All"
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.diffmap(adata, color='fucci_time', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsFucciPseudotime3d.pdf")
sc.pl.umap(adata, color=["fucci_time"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqFucciPseudotime.pdf")

#%% Subset based on phase and separate out the weird nodule
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
moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, 100)
max_exp_inds = np.argmax(moving_averages, 0)
moving_avg_peak_inds = np.argsort(max_exp_inds)
moving_avg_sort = np.take(moving_averages, moving_avg_peak_inds, axis=1)
gene_id_sort = np.take(np.array(adata.var))

# take the moving average and sort by max time


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

#%%

# take the moving average and sort by max time
moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, 100)
max_exp_inds = np.argmax(moving_averages, 0)
moving_avg_peak_inds = np.argsort(max_exp_inds)
moving_avg_sort = np.take(moving_averages, moving_avg_peak_inds, axis=1)
gene_id_sort = np.take(np.array(adata.var_names), moving_avg_peak_inds)
plt.figure(figsize=(10,10))
allccd_transcript_regulated = np.array(pd.read_csv("output/allccd_transcript_regulated.csv")["gene"])
plt.imshow(moving_avg_sort.T[np.isin(gene_id_sort, allccd_transcript_regulated)], interpolation="nearest")
plt.tight_layout()
plt.legend()
# plt.savefig(os.path.join("figures",'RNA_heatmap.pdf'), transparent=True)
plt.show()
plt.close()

#%%

#create localization matrix
#sort the coefficients based on the max value point
loc_mat = np.zeros([len(model_free_list),len(loc_set)])
for i,response in enumerate(model_free_list):
    curr_well = response[0]
    plate = int(curr_well.split("_")[-1])
    print(curr_well)
    print(plate_well_to_ab[plate][curr_well])
    curr_locs = plate_well_to_ab[plate][curr_well][3]
    loc_mat[i,curr_locs] = 1

model_free_list.sort(key=operator.itemgetter(3))
prob_interact = np.zeros([len(model_free_list),len(model_free_list)])
sum_rows = []
prot_list = []
sorted_gene_array = []
sorted_ensg_array = []
sorted_maxpol_array = []
sorted_expression_table_ccd = []
sorted_expression_table_all = []
for i,response in enumerate(model_free_list):
    sorted_gene_array.append(response[4]) # this is the expression values
    sorted_maxpol_array.append(response[3])
    sorted_ensg_array.append(response[6])
    for i, exp in enumerate(response[4]):
        sorted_expression_table_ccd.append([response[6], (i + 1) / len(response[4]), exp])
    curr_well = response[0]
    #set all other rows that share a loc to 1
    curr_locs = loc_mat[i,:]
    match_rows = np.greater(np.sum(loc_mat*curr_locs,axis=1),0)*response[3]
    sum_rows.append(np.sum(match_rows))
    prob_interact[i,:] = match_rows
    prot_list.append(response[1])

for i, response in enumerate(model_free_list_all):
    for i, exp in enumerate(response[4]):
        sorted_expression_table_all.append([response[6], (i + 1) / len(response[4]), exp])   

fig, ax = plt.subplots(figsize=(10, 10))
sc = ax.imshow(sorted_gene_array, interpolation='nearest')
np.savetxt("output/temporal_protein_expression_ccd.tsv", npv.array(sorted_expression_table_ccd), fmt="%s", delimiter="\t")
np.savetxt("output/temporal_protein_expression_all.tsv", np.array(sorted_expression_table_all), fmt="%s", delimiter="\t")

#Arange the x ticks
xtick_labels = [str(np.around(x * TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)]
xtick_labels = xtick_labels#+['G1/S','S/G2']
xphase_labels = ['G1/S','S/G2']
my_xticks = np.arange(-.5, 20, 2)
num_ticks = 20

phase_trans = np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])

ax.set_xticks(my_xticks,minor=True)
ax.set_xticklabels(xtick_labels,minor=True)
ax.set_xticks(phase_trans, minor=False)
ax.set_xticklabels(xphase_labels, minor=False)
ax.tick_params(length=12)
plt.xticks(size=12,fontname='Arial')
plt.yticks(size=10,fontname='Arial')

#Do the y ticks
ytick_locs = [prot_list.index(p) for p in HIGHLIGHTS]
# ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
ax.set_yticks(ytick_locs,minor=False)
ax.set_yticklabels(HIGHLIGHTS,minor=False)
# ax.set_yticks(ytick_locs_minor,minor=True)
# ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)
ax.tick_params(direction='out', length=12, width=2, colors='k', 
                axis='x',which='major')
ax.tick_params(direction='out', length=56, width=1, colors='k', 
                axis='y',which='major')
ax.set_aspect('auto')
plt.xlabel('Division Cycle, hrs',size=20,fontname='Arial')
plt.ylabel('Gene',size=20,fontname='Arial')

divider1 = make_axes_locatable(ax)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(sc,cax = cax1)
cbar.set_label('Relative expression',fontname='Arial',size=20)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(os.path.join("output_devin",'sorted_heatmap21_sw30_take4.pdf'), transparent=True)
plt.show()