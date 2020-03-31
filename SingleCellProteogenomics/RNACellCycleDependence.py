# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 19:31:42 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, FucciCellCycle, FucciPseudotime, RNADataPreparation
from SingleCellProteogenomics import mvavg

WINDOW = 100
PERMUTATIONS = 10000
MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM = 0.08
fucci = FucciCellCycle()

def boxplot_result(g1, s, g2, outfolder, ensg):
    '''Make boxplots of RNA expression by phase'''
    if not os.path.exists(f"{outfolder}_png"): os.mkdir(f"{outfolder}_png")
    if not os.path.exists(f"{outfolder}_pdf"): os.mkdir(f"{outfolder}_pdf")
    mmmm = np.concatenate((g1, s, g2))
    cccc = (["G1"] * len(g1))
    cccc.extend(["G1/S"] * len(s))
    cccc.extend(["G2"] * len(g2))
    boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=True)
    boxplot.set_xlabel("", size=36,fontname='Arial')
    boxplot.set_ylabel("Normalized Log10(TPM)", size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    plt.title("")
    plt.savefig(f"{outfolder}_png/GaussianClusteringProtein_{ensg}.png")
    plt.savefig(f"{outfolder}_pdf/GaussianClusteringProtein_{ensg}.pdf")
    plt.close()

def plot_expression_umap(genelist, outfolder):
    '''
    Displaying relative expression on the UMAP (using unlogged expression values)
    Output: a folder of UMAPs for each gene in a list
    '''
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        adata.obs[gene] = normalized_exp_data
        sc.pl.umap(adata, color=gene, show=False, save=True)
        shutil.move("figures/umap.pdf", f"{outfolder}/{gene}.pdf")
        plt.close()
        adata.obs.drop(gene, 1)

def plot_expression_facs(genelist, normalized_exp_data, phases, var_names, outfolder):
    '''
    Displaying relative expression on the FUCCI plot (log-log intensity of FUCCI markers)
    Output: a folder of FUCCI plots for each gene in a list
    '''
    phases_validInt = phases[pd.notnull(phases.Green530) & pd.notnull(phases.Red585)]
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        nexp = normalized_exp_data[:,list(var_names).index(gene)]
        plt.scatter(phases_validInt["Green530"], phases_validInt["Red585"], nexp)
        plt.tight_layout()
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel("Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
        plt.title(gene,size=20,fontname='Arial')
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

def analyze_ccd_variation_rna(adata):
    expression_data = adata.X # log normalized
    normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
    fucci_time_inds = np.argsort(adata.obs["fucci_time"])
    norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
    moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, WINDOW)
    mvavg_xvals = mvavg(adata.obs["fucci_time"][fucci_time_inds], WINDOW)
    cell_cycle_variance = np.var(moving_averages, 0)
    total_variance = np.var(norm_exp_sort, 0)
    total_gini = np.apply_along_axis(utils.gini, 0, norm_exp_sort)
    percent_ccd_variance = cell_cycle_variance / total_variance
    avg_expression = np.median(norm_exp_sort, 0)

    # randomize and calculate the mean difference in percent variances from random
    percent_ccd_variance_rng, mean_diff_from_rng = [],[]
    perms = np.asarray([np.random.permutation(len(adata.obs)) for nnn in np.arange(PERMUTATIONS)])
    if not os.path.exists("output/pickles/percent_ccd_variance_rng.npy") or not os.path.exists("output/pickles/mean_diff_from_rng.npy"):
        norm_exp_sort_perm = np.asarray([np.take(normalized_exp_data, perm, axis=0) for perm in perms])
        moving_averages_perm = np.apply_along_axis(mvavg, 1, norm_exp_sort_perm, WINDOW)
        percent_ccd_variance_rng = np.var(moving_averages_perm, axis=1) / np.var(norm_exp_sort_perm, axis=1)
        for iii, perm in enumerate(perms):
            if iii % 50 == 0: print(f"permutation {iii}")
            norm_exp_sort_perm = np.take(normalized_exp_data, perm, axis=0)
            moving_averages_perm = np.apply_along_axis(mvavg, 0, norm_exp_sort_perm, WINDOW)
            percent_ccd_variance_rng.append(
                    np.var(moving_averages_perm, axis=0) / np.var(norm_exp_sort_perm, axis=0))
        utils.np_save_overwriting("output/pickles/percent_ccd_variance_rng.npy", percent_ccd_variance_rng)
        utils.np_save_overwriting("output/pickles/mean_diff_from_rng.npy", mean_diff_from_rng)
    else: 
        percent_ccd_variance_rng = np.load("output/pickles/percent_ccd_variance_rng.npy", allow_pickle=True)
        mean_diff_from_rng = np.load("output/pickles/mean_diff_from_rng.npy", allow_pickle=True)
    percent_ccd_variance_rng = np.asarray(percent_ccd_variance_rng)
    mean_diff_from_rng = np.mean((percent_ccd_variance - percent_ccd_variance_rng).T, 1)

    # Statistical testing based on randomization analysis
    alpha_ccd = 0.01
    pass_meandiff = mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM
    ccd_var_comp_rng_wilcoxp = np.apply_along_axis(scipy.stats.wilcoxon, 1, (percent_ccd_variance - percent_ccd_variance_rng).T, None, "wilcox", False, "greater").T[1].T
    eq_percvar_adj, pass_eq_percvar_adj = bonf(alpha_ccd, ccd_var_comp_rng_wilcoxp)
    gtpass_eq_percvar_adj = pass_eq_percvar_adj & (percent_ccd_variance > np.median(percent_ccd_variance_rng, axis=0))
    
    return percent_ccd_variance, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms

def figures_ccd_analysis_rna(adata, percent_ccd_variance, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, wp_ensg, ccd_comp, ccd_regev_filtered):
    '''Print some figures of merit for the RNA CCD analysis'''
    print(f"{sum(pass_meandiff)}: number of transcripts that pass CCD by mean percvar diff from random > {MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM}")
    print(f"{sum(np.isin(adata.var_names[mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM], wp_ensg[ccd_comp]))}: number of CCD transcripts that are also CCD proteins")
    print(f"{sum(np.isin(adata.var_names[mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM], ccd_regev_filtered)) / len(ccd_regev_filtered)}: percent of Regev CCD transcripts that are called CCD with this analysis")

    plt.figure(figsize=(10,10))
    plt.scatter(percent_ccd_variance, mean_diff_from_rng, c=pass_meandiff, cmap="bwr_r")
    plt.hlines(MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM, np.min(percent_ccd_variance), np.max(percent_ccd_variance), color="gray")
    plt.xlabel("Percent Variance Explained by Cell Cycle")
    plt.ylabel("Mean Difference from Random")
    plt.savefig("figures/MedianDiffFromRandom_RNA.png")
    plt.savefig("figures/MedianDiffFromRandom_RNA.pdf")
    plt.show()
    plt.close()

    eq_percvar_adj_nextafter = np.nextafter(eq_percvar_adj, eq_percvar_adj + 1)
    plt.figure(figsize=(10,10))
    plt.scatter(mean_diff_from_rng, -np.log10(eq_percvar_adj_nextafter), c=pass_meandiff, cmap="bwr_r")
    plt.vlines(MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM, np.min(-np.log10(eq_percvar_adj_nextafter)), np.max(-np.log10(eq_percvar_adj_nextafter)), color="gray")
    plt.xlabel("Mean Difference from Random")
    plt.ylabel("-log10 adj p-value from randomization")
    plt.savefig("figures/MedianDiffFromRandomVolcano_RNA.png")
    plt.savefig("figures/MedianDiffFromRandomVolcano_RNA.pdf")
    plt.show()
    plt.close()

def mvavg_plots_pergene(adata, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals):
    for iii, ensg in enumerate(adata.var_names):
        plt.close('all')
        if iii % 500 == 0: print(f"well {iii} of {len(adata.var_names)}")  
        windows = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(len(adata.obs["fucci_time"][fucci_time_inds]) - WINDOW + 1)])
        MovingAverages.temporal_mov_avg_rna(adata.obs["fucci_time"][fucci_time_inds], norm_exp_sort[:,iii], mvavg_xvals, moving_averages[:,iii], windows, f"figures/RNAPseudotimes", f"{ensg}")

def ccd_analysis_of_spikeins(adata_spikeins, perms):
    '''ERCC spikeins were used as an internal control. We can use them to get an idea of the noise for this analysis.'''
    expression_data_spike = adata_spikeins.X # log normalized
    normalized_exp_data_spike = (expression_data_spike.T / np.max(expression_data_spike, axis=0)[:,None]).T
    fucci_time_inds_spike = np.argsort(adata_spikeins.obs["fucci_time"])
    # fucci_time_sort_spike = np.take(np.array(adata_spikeins.obs["fucci_time"]), fucci_time_inds_spike)
    norm_exp_sort_spike = np.take(normalized_exp_data_spike, fucci_time_inds_spike, axis=0)
    moving_averages_spike = np.apply_along_axis(mvavg, 0, norm_exp_sort_spike, 100)
    cell_cycle_variance_spike = np.apply_along_axis(np.var, 0, moving_averages_spike)
    total_variance_spike = np.apply_along_axis(np.var, 0, norm_exp_sort_spike)
    total_cv_spike = np.apply_along_axis(scipy.stats.variation, 0, norm_exp_sort_spike)
    percent_ccd_variance_spike = cell_cycle_variance_spike / total_variance_spike
    # avg_expression_spike = np.apply_along_axis(np.median, 0, norm_exp_sort_spike)
    print(f"mean +/- stdev of spike-in variance: {np.mean(total_variance_spike)} +/- {np.std(total_variance_spike)}")
    print(f"median of spike-in variance: {np.median(total_variance_spike)}")
    print(f"mean +/- stdev of spike-in variance explained by cell cycle: {np.mean(percent_ccd_variance_spike)} +/- {np.std(percent_ccd_variance_spike)}")
    print(f"median of spike-in variance explained by cell cycle: {np.median(percent_ccd_variance_spike)}")

    percent_ccd_variance_rng_spike = []
    for iii, perm in enumerate(perms):
        if iii % 50 == 0: print(f"permutation {iii}")
        norm_exp_sort_perm_spike = np.take(normalized_exp_data_spike, perm, axis=0)
        moving_averages_perm_spike = np.apply_along_axis(mvavg, 0, norm_exp_sort_perm_spike, WINDOW)
        percent_ccd_variance_rng_spike.append(
                np.var(moving_averages_perm_spike, axis=0) / np.var(norm_exp_sort_perm_spike, axis=0))
    percent_ccd_variance_rng_spike = np.asarray(percent_ccd_variance_rng_spike)
    mean_diff_from_rng_spike = np.mean((percent_ccd_variance - percent_ccd_variance_rng).T, 1)
    print(f"mean +/- stdev of spike-in mean difference from random percent variance: {np.mean(mean_diff_from_rng_spike)} +/- {np.std(mean_diff_from_rng_spike)}")