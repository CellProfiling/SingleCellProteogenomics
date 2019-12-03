#%% Imports
from imports import *
import numpy as np
from scipy import stats
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms" # so that the gene-specific results scales match for cross-gene comparisons
print("reading scRNA-Seq data")
biotype_to_use="protein_coding"
# biotype_to_use="lincRNA" # there are only 10 genes marked lincRNA, so need to analyze isoform data, I think
adata, phases = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=False, biotype_to_use=biotype_to_use)
adata, phasesfilt = qc_filtering(adata, do_log_normalize=False, do_remove_blob=True)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Cell cycle regulated ANOVA
# Idea: separate the expression of each gene into G1/S/G2 phases, and use ANOVA to test if there's variation
# Execution: scipy.stats.f_oneway and multiple testing correction
# Output: table of all genes and ANOVA test results

# expression_data = np.exp(adata.X) - 1
# normalized_exp_data = expression_data / np.max(expression_data)
expression_data = adata.X
normalized_exp_data = expression_data / np.max(expression_data)
stages = np.array(adata.obs["phase"])
g1_exp = np.take(normalized_exp_data, np.nonzero(stages == "G1")[0], axis=0)
s_exp = np.take(normalized_exp_data, np.nonzero(stages == "S-ph")[0], axis=0)
g2_exp = np.take(normalized_exp_data, np.nonzero(stages == "G2M")[0], axis=0)
tests_fp = [scipy.stats.f_oneway(g1_exp[:,geneidx], s_exp[:,geneidx], g2_exp[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
pvals = [p for (F, p) in tests_fp]

# benjimini-hochberg multiple testing correction
# source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
alpha = 0.01
def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection'''
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)

pvals_sortind = np.argsort(pvals)
pvals_sorted = np.take(pvals, pvals_sortind)
ecdffactor = _ecdf(pvals_sorted)
reject = pvals_sorted <= ecdffactor*alpha
reject = pvals_sorted <= ecdffactor*alpha
if reject.any():
    rejectmax = max(np.nonzero(reject)[0])
    reject[:rejectmax] = True
pvals_corrected_raw = pvals_sorted / ecdffactor
pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
del pvals_corrected_raw
pvals_corrected[pvals_corrected>1] = 1
pvals_corrected_ = np.empty_like(pvals_corrected)

# deal with sorting
pvals_corrected_[pvals_sortind] = pvals_corrected
del pvals_corrected
reject_ = np.empty_like(reject)
reject_[pvals_sortind] = reject

# BenjiHoch is actually pretty liberal for this dataset. What about bonferroni?
# bonferroni MTC
alphaBonf = alpha / float(len(pvals))
rejectBonf = pvals_sorted <= alphaBonf
pvals_correctedBonf = pvals_sorted * float(len(pvals))
pvals_correctedBonf_unsorted = np.empty_like(pvals_correctedBonf)
pvals_correctedBonf_unsorted[pvals_sortind] = pvals_correctedBonf
rejectBonf_unsorted = np.empty_like(rejectBonf)
rejectBonf_unsorted[pvals_sortind] = rejectBonf

anova_tests = pd.DataFrame(
    {"gene" : adata.var_names,
    "pvalue" : pvals,
    "pvaladj_BH" : pvals_corrected_,
    "reject_BH" : reject_,
    "pvaladj_B" : pvals_correctedBonf_unsorted,
    "reject_B" : rejectBonf_unsorted})
anova_tests.to_csv(f"output/transcript_regulation{biotype_to_use}.csv")

#%% [markdown]
## Summary of ANOVA analysis
#
# I looked for transcript regulation across the cell cycle  on all 17,199 genes that made it through filtering 
# (using ANOVA with Benjimini-Hochberg multiple testing correction).  
# Of those, about a third have evidence of cell cycle regulation at the transcript level (at 1% FDR). 
# Of the gene set from the Regev study (single cell transcriptomics), all but one (99%) of the genes 
# have evidence of transcript regulation over the cell cycle. An example of one of the surprises is above, 
# where geminin is significant (rejecting the hypothesis that all stages have the same expression level), 
# even while this effect is subtle in the pseudotime analysis. It also seems to make sense that 
# G1 and S phases have higher expression, while G2M has lower expression of geminin.
#
# In Diana's 464 cell-cycle dependent genes, 197 had evidence of cell-cycle regulation at the 
# transcript level (48% of the 410 genes that weren't filtered in scRNA-Seq preprocessing). 
# Strangely, of the 890 non cell-cycle dependent genes, including 811 that made it through 
# scRNA-Seq preprocessing, 362 genes (45%) showed evidence for cell cycle dependence at the 
# transcript level.
#
# I tried a Bonferroni correction (alpha=0.01), which is a bit more conservative, 
# and 97% of the Regev genes, 27% of Diana's CCD proteins, and 20% of Diana's 
# non-CCD proteins showed cell cycle regulation for transcripts.

#%% Expression boxplots
# Idea: Use boxplots to display the variance and overlay ANOVA results for significance
# Execution: Separate G1/S/G2 and generate boxplots
# Output: Boxplots
def format_p(p):
    '''3 decimal places, scientific notation'''
    return '{:0.3e}'.format(p)

def plot_expression_boxplots(genelist, outanova, outfolder):
    notfound = []
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        if gene not in adata.var_names:
            notfound.append(gene)
            continue
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        expression_stage = pd.DataFrame({gene : normalized_exp_data, "stage" : adata.obs["phase"]})
        exp_stage_filt = expression_stage[expression_stage.stage != "nan"].reset_index(drop=True)
        boxplot = exp_stage_filt.boxplot(gene, by="stage", figsize=(12, 8))
        boxplot.set_xlabel("Cell Cycle Stage", size=36,fontname='Arial')
        boxplot.set_ylabel("Normalized Transcript Expression", size=36,fontname='Arial')
        boxplot.tick_params(axis="both", which="major", labelsize=16)
        qqq = anova_tests[anova_tests.gene == gene].iloc[0]["pvaladj"]
        rej = anova_tests[anova_tests.gene == gene].iloc[0]["reject"]
        boxplot.set_title(f"{gene} (Q_bh = {format_p(qqq)})",size=36,fontname='Arial')
        boxplot.get_figure().savefig(f"{outfolder}/{gene}_boxplot.png")
        plt.close()
    pd.DataFrame({"gene" : notfound}).to_csv(f"{outfolder}/GenesFilteredInScRNASeqAnalysis.csv")

# plot_expression_boxplots(ccd_regev_filtered, "ccd_regev_filtered", "figures/RegevGeneBoxplots")
# plot_expression_boxplots(ccd_filtered, "ccd_filtered", "figures/DianaCcdGeneBoxplots")
# plot_expression_boxplots(nonccd_filtered, "nonccd_filtered", "figures/DianaNonCcdGeneBoxplots")

# Output the stage in which there is peak expression for each gene in the known set
# expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
# normalized_exp_data = expression_data / np.max(expression_data)
# expression_stage = pd.DataFrame({gene : normalized_exp_data, "stage" : adata.obs["phase"]})
# exp_stage_filt = expression_stage[expression_stage.stage != "nan"].reset_index(drop=True)


#%% Plotting variances of gene expression
# Idea: Using the spread as a measure of variability, is there higher variability for CCD genes (likely due to cyclicity)?
# Execution: mean and stdevs
# Output: Scatter of stdev vs mean expression; histogram of stdevs
means = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
variances = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
mean_regev = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_regev_filtered]
variances_regev = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_regev_filtered]
mean_dianaccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
variances_dianaccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
mean_diananonccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]
variances_diananonccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x = means, y = variances, label="All Genes")
ax1.scatter(x = mean_regev, y = variances_regev, label="Regev CCD Genes")
ax1.scatter(x = mean_dianaccd, y = variances_dianaccd, label="Fucci CCD Genes")
ax1.scatter(x = mean_diananonccd, y = variances_diananonccd, label="Fucci Non-CCD Genes")
plt.legend(loc="upper left")
plt.xlabel("Mean Expression")
plt.ylabel("Stdev Expression") 
plt.tight_layout()
plt.savefig(f"figures/stdev_expression{biotype_to_use}.png")
plt.show()
plt.close()

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

fig = plt.figure()
ax1 = fig.add_subplot(111)
bins=np.histogram(np.hstack((variances, variances_regev, variances_dianaccd, variances_diananonccd)), bins=40)[1] #get the bin edges
ax1.hist(variances_regev, bins=bins, weights=weights(variances_regev), 
    label="Regev CCD Genes")
ax1.hist(variances_dianaccd, bins=bins, weights=weights(variances_dianaccd),
    label="Fucci CCD Genes")
ax1.hist(variances_diananonccd, bins=bins, weights=weights(variances_diananonccd), 
    label="Fucci Non-CCD Genes")
ax1.hist(variances, bins=bins, weights=weights(variances), 
    label="All Genes")
plt.legend(loc="upper right")
plt.xlabel("Stdev Expression")
plt.ylabel("Count, Normalized to 1")
plt.tight_layout()
plt.savefig(f"figures/stdev_expression_hist_{biotype_to_use}.png")
plt.show()
plt.close()
 

#%% Redo inputs:
# Log normalize and use RPKMS because we're comparing genes in the next cell
count_or_rpkm = "Tpms"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, False, biotype_to_use)
adata, phasesfilt = qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)

adata_spikeins, phases_spikeins = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=True, biotype_to_use="")
sc.pp.filter_genes(adata_spikeins, min_cells=100)
print(f"data shape after filtering: {adata_spikeins.X.shape}")

#%% Moving average percent variance
# Idea: What percentage of variance in transcript expression is due to the cell cycle?
# Execution: Calculate variance due to cell cycle as residuals from a moving average
# Output: Scatters of the percent and total variances for each gene 
def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')

expression_data = adata.X # log normalized
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, 100)
cell_cycle_variance = np.apply_along_axis(np.var, 0, moving_averages)
total_variance = np.apply_along_axis(np.var, 0, norm_exp_sort)
total_cv = np.apply_along_axis(stats.variation, 0, norm_exp_sort)
percent_ccd_variance = cell_cycle_variance / total_variance
avg_expression = np.apply_along_axis(np.median, 0, norm_exp_sort)

def plot_variances(total_var, percent_var, expression_color, title, file_tag):
    '''Plots percent variances from cell line against total variance'''
    plt.figure(figsize=(10,10))
    plt.scatter(total_var, percent_var, c=expression_color)
    plt.xlabel("Total Variance", size=36,fontname='Arial')
    plt.ylabel("Percent Variance Due to Cell Cycle",size=36,fontname='Arial')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.title(title, size=36, fontname='Arial')
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 24
    cbar.ax.set_ylabel("Median Log-Normalized RNA-Seq Counts", rotation=270,size=24,fontname='Arial')
    plt.savefig(f"figures/VariancePlot{file_tag}_{biotype_to_use}.png")
    plt.show()
    plt.close()

regevccdgenes = np.isin(adata.var_names, ccd_regev_filtered)
dianaccdgenes = np.isin(adata.var_names, ccd_filtered)
diananonccdgenes = np.isin(adata.var_names, nonccd_filtered)
plot_variances(total_variance, percent_ccd_variance, avg_expression, "All Genes", "All")
plot_variances(total_variance[regevccdgenes], percent_ccd_variance[regevccdgenes], avg_expression[regevccdgenes], "Regev CCD Genes", "RegevCcd")
plot_variances(total_variance[dianaccdgenes], percent_ccd_variance[dianaccdgenes], avg_expression[dianaccdgenes], "Diana CCD Genes", "DianaCcd")
plot_variances(total_variance[diananonccdgenes], percent_ccd_variance[diananonccdgenes], avg_expression[diananonccdgenes], "Diana Non-CCD Genes", "DianaNonCCD")

# Displaying the ANOVA results on the percent variances
def plot_variances_tf(total_var, percent_var, expression_color, title, file_tag, percent_var_cutoff, total_var_cutoff):
    '''Plots percent variances from cell line against total variance'''
    # plt.figure(figsize=(10,10))
    plt.scatter(total_var, percent_var, c=expression_color, cmap="bwr_r")
    plt.axhline(percent_var_cutoff)
    plt.axvline(total_var_cutoff)
    plt.xlabel("Total Variance", size=24,fontname='Arial')
    plt.ylabel("Percent Variance Due to Cell Cycle",size=24,fontname='Arial')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlim(-0.005, 0.1)
    plt.ylim(-0.05, 0.8)
    plt.title(title, size=24, fontname='Arial')
    plt.legend()
    # plt.savefig(f"figures/VarianceSignificancePlot{file_tag}.png")
    # plt.show()
    # plt.close()

# plot_variances_tf(total_variance, percent_ccd_variance, rejectBonf_unsorted, "All Genes", "All")
# plot_variances_tf(total_variance[regevccdgenes], percent_ccd_variance[regevccdgenes], rejectBonf_unsorted[regevccdgenes], "Regev CCD Genes", "RegevCcd")
# plot_variances_tf(total_variance[dianaccdgenes], percent_ccd_variance[dianaccdgenes], rejectBonf_unsorted[dianaccdgenes], "Diana CCD Genes", "DianaCcd")
# plot_variances_tf(total_variance[diananonccdgenes], percent_ccd_variance[diananonccdgenes], rejectBonf_unsorted[diananonccdgenes], "Diana Non-CCD Genes", "DianaNonCCD")

# What are these low percent variance ones that are significant?
print("What are these low percent variance ones that are significant?")
low_variance_signif = (percent_ccd_variance < 0.03) & dianaccdgenes & rejectBonf_unsorted
print(np.array(adata.var_names)[low_variance_signif])

# Choose a cutoff for total percent variance based on the spike-in genes
expression_data_spike = adata_spikeins.X # log normalized
normalized_exp_data_spike = (expression_data_spike.T / np.max(expression_data_spike, axis=0)[:,None]).T
fucci_time_inds_spike = np.argsort(adata_spikeins.obs["fucci_time"])
# fucci_time_sort_spike = np.take(np.array(adata_spikeins.obs["fucci_time"]), fucci_time_inds_spike)
norm_exp_sort_spike = np.take(normalized_exp_data_spike, fucci_time_inds_spike, axis=0)
moving_averages_spike = np.apply_along_axis(mvavg, 0, norm_exp_sort_spike, 100)
cell_cycle_variance_spike = np.apply_along_axis(np.var, 0, moving_averages_spike)
total_variance_spike = np.apply_along_axis(np.var, 0, norm_exp_sort_spike)
total_cv_spike = np.apply_along_axis(stats.variation, 0, norm_exp_sort_spike)
percent_ccd_variance_spike = cell_cycle_variance_spike / total_variance_spike
# avg_expression_spike = np.apply_along_axis(np.median, 0, norm_exp_sort_spike)
print(f"mean +/- stdev of spike-in variance: {np.mean(total_variance_spike)} +/- {np.std(total_variance_spike)}")
print(f"median of spike-in variance: {np.median(total_variance_spike)}")
print(f"mean +/- stdev of spike-in variance explained by cell cycle: {np.mean(percent_ccd_variance_spike)} +/- {np.std(percent_ccd_variance_spike)}")
print(f"mean of spike-in variance explained by cell cycle: {np.median(percent_ccd_variance_spike)}")

# Okay, we're going to take the ones that passed ANOVA and have > 10% variance explained
total_var_cutoff = np.mean(total_variance_spike) + 1 * np.std(total_variance_spike)
percent_var_cutoff = np.mean(percent_ccd_variance_spike) + 1 * np.std(percent_ccd_variance_spike)
print(f"using total variance cutoff {total_var_cutoff} and percent CCD variance cutoff {percent_var_cutoff}")
plt.figure(figsize=(15,15))
plt.subplot(221)
plot_variances_tf(total_variance, percent_ccd_variance, rejectBonf_unsorted & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff), "All Genes", "All", percent_var_cutoff, total_var_cutoff)
plt.subplot(222)
plot_variances_tf(total_variance[regevccdgenes], percent_ccd_variance[regevccdgenes], (rejectBonf_unsorted[regevccdgenes]) & (percent_ccd_variance[regevccdgenes] > percent_var_cutoff) & (total_variance[regevccdgenes] > total_var_cutoff), "Regev CCD Genes", "RegevCcd", percent_var_cutoff, total_var_cutoff)
plt.subplot(223)
plot_variances_tf(total_variance[dianaccdgenes], percent_ccd_variance[dianaccdgenes], (rejectBonf_unsorted[dianaccdgenes]) & (percent_ccd_variance[dianaccdgenes] > percent_var_cutoff) & (total_variance[dianaccdgenes] > total_var_cutoff), "Diana CCD Genes", "DianaCcd", percent_var_cutoff, total_var_cutoff)
plt.subplot(224)
plot_variances_tf(total_variance[diananonccdgenes], percent_ccd_variance[diananonccdgenes], (rejectBonf_unsorted[diananonccdgenes]) & (percent_ccd_variance[diananonccdgenes] > percent_var_cutoff) & (total_variance[diananonccdgenes] > total_var_cutoff), "Diana Non-CCD Genes", "DianaNonCCD", percent_var_cutoff, total_var_cutoff)
plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}.png")
plt.show()

# Let's output those
gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
gene_ids = list(gene_info["gene_id"])
gene_names = list(gene_info["name"])
gene_id_name = dict([(gene_ids[idxx], gene_names[idxx]) for idxx in range(len(gene_info))])
percent_variance_tests = pd.DataFrame(
    {"gene" : adata.var_names, 
    "name" : [gene_id_name[x] if x in gene_id_name else "" for x in adata.var_names],
    "pvalue" : pvals, 
    "pvaladj_BH" : pvals_corrected_, 
    "reject_BH" : reject_, 
    "pvaladj_B" : pvals_correctedBonf_unsorted, 
    "reject_B" : rejectBonf_unsorted,
    "percent_variance" : percent_ccd_variance,
    "total_variance" : total_variance,
    "significant" : (rejectBonf_unsorted) & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff),
    "regev_ccd" : regevccdgenes,
    "diana_ccd" : dianaccdgenes,
    "diana_nonccd" : diananonccdgenes})
percent_variance_tests.to_csv(f"output/transcript_regulation{biotype_to_use}.csv")

# And output the CCD genes that aren't in diana's CCD gene set
percent_variance_tests[percent_variance_tests.significant & ~percent_variance_tests.diana_ccd].to_csv("output/transcript_regulation_significant_ccd_notindianasset.csv")

# And keep track of the ccd genes with and without transcript regulation
def np_save_overwriting(fn, arr):
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)

ccd_transcript_regulated = rejectBonf_unsorted & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff)
dianaccd_transcript_regulated = dianaccdgenes & rejectBonf_unsorted & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff)
dianaccd_nontranscript_regulated = dianaccdgenes & ~dianaccd_transcript_regulated
ccd_transcript_regulated_names = np.array(adata.var_names)[ccd_transcript_regulated]
dianaccd_transcript_regulated_names = np.array(adata.var_names)[dianaccd_transcript_regulated]
dianaccd_nontranscript_regulated_names = np.array(adata.var_names)[dianaccd_nontranscript_regulated]
np_save_overwriting("output/ccd_transcript_regulated.npy", ccd_transcript_regulated)
np_save_overwriting("output/dianaccd_transcript_regulated.npy", dianaccd_transcript_regulated)
np_save_overwriting("output/dianaccd_nontranscript_regulated.npy", dianaccd_nontranscript_regulated)
pd.DataFrame({"gene" : ccd_transcript_regulated_names}).to_csv("output/allccd_transcript_regulated.csv")
pd.DataFrame({"gene" : dianaccd_transcript_regulated_names}).to_csv("output/ccd_transcript_regulated.csv")
pd.DataFrame({"gene" : dianaccd_nontranscript_regulated_names}).to_csv("output/ccd_nontranscript_regulated.csv")
pd.DataFrame({"gene" : adata.var_names}).to_csv("output/gene_names.csv")


#%%
