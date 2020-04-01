# -*- coding: utf-8 -*-
"""
Methods for assessing cell cycle dependence of RNA abundances in single cells.
-  Percent variance attributed to the cell cycle was calculated using the (variance of moving average / total variance)
-  Randomization analysis was used to determine statistical significance of high percent variances due to the cell cycle

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, MovingAverages, FucciCellCycle, FucciPseudotime, RNADataPreparation

WINDOW = 100 # Number of points for moving average window for protein analysis
PERMUTATIONS = 10000 # Number of permutations used for randomization analysis
MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM = 0.08 # Cutoff used for percent additional variance explained by the cell cycle than random

fucci = FucciCellCycle.FucciCellCycle() # Object representing FUCCI cell cycle phase durations

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

def plot_expression_umap(adata, genelist, outfolder):
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

def plot_expression_boxplots(adata, genelist, bulk_phase_tests, outfolder):
    '''Generate plots of RNA expression by phase in a boxplots'''
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
        qqq = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["pvaladj"]
        rej = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["reject"]
        boxplot.set_title(f"{gene} (Q_bh = {format_p(qqq)})",size=36,fontname='Arial')
        boxplot.get_figure().savefig(f"{outfolder}/{gene}_boxplot.png")
        plt.close()
    pd.DataFrame({"gene" : notfound}).to_csv(f"{outfolder}/GenesFilteredInScRNASeqAnalysis.csv")

def panel_variances_tf(total_var, percent_var, expression_color, title, file_tag):
    '''Displaying the bulk_phase results on the percent variances: plot percent variances from cell line against total variance'''
    plt.scatter(total_var, percent_var, c=expression_color, cmap="bwr_r")
    plt.xlabel("Gini of Gene Expression", size=24,fontname='Arial')
    plt.ylabel("Percent Variance Due to Cell Cycle",size=24,fontname='Arial')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlim(-0.005, 1)
    plt.ylim(-0.05, 0.8)
    plt.title(title, size=24, fontname='Arial')
    plt.legend()

def plot_overall_and_ccd_variances(adata, biotype_to_use, total_gini, percent_ccd_variance, pass_meandiff, ccdprotein, nonccdprotein, regevccdgenes):
    '''Plots displaying CCD variance and overall variances'''
    
    # Paneled plot for variances of several groups
    plt.figure(figsize=(15,15))
    plt.subplot(221)
    panel_variances_tf(total_gini, percent_ccd_variance, pass_meandiff, "All Genes", "All")
    plt.subplot(222)
    panel_variances_tf(total_gini[regevccdgenes], percent_ccd_variance[regevccdgenes], pass_meandiff[regevccdgenes], "Regev CCD Genes", "RegevCcd")
    plt.subplot(223)
    panel_variances_tf(total_gini[ccdprotein], percent_ccd_variance[ccdprotein], pass_meandiff[ccdprotein], "CCD Proteins", "DianaCcd")
    plt.subplot(224)
    panel_variances_tf(total_gini[nonccdprotein], percent_ccd_variance[nonccdprotein], pass_meandiff[nonccdprotein], "Non-CCD Proteins", "DianaNonCCD")
    plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}.pdf")
    plt.show()
    plt.close()

    # A second plot with just the variance for all
    panel_variances_tf(total_gini, percent_ccd_variance, pass_meandiff, "All Genes", "All")
    plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}_allgenes.pdf")
    plt.show()
    plt.close()

def analyze_ccd_variation_by_phase_rna(adata, normalized_exp_data, biotype_to_use):
    stages = np.array(adata.obs["phase"])
    g1_exp = np.take(normalized_exp_data, np.nonzero(stages == "G1")[0], axis=0)
    s_exp = np.take(normalized_exp_data, np.nonzero(stages == "S-ph")[0], axis=0)
    g2_exp = np.take(normalized_exp_data, np.nonzero(stages == "G2M")[0], axis=0)
    tests_fp = [scipy.stats.kruskal(g1_exp[:,geneidx], s_exp[:,geneidx], g2_exp[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
    pvals = [p for (F, p) in tests_fp]
    pvals_corrected_BH, reject_BH = utils.benji_hoch(0.01, pvals)
    pvals_correctedBonf, rejectBonf = utils.bonf(0.01, pvals)
    bulk_phase_tests = pd.DataFrame({
        "gene" : adata.var_names,
        "pvalue" : pvals,
        "pvaladj_BH" : pvals_corrected_BH,
        "reject_BH" : reject_BH,
        "pvaladj_B" : pvals_correctedBonf,
        "reject_B" : rejectBonf
        })
    bulk_phase_tests.to_csv(f"output/phase_clustered_transcript_CCD_analysis_{biotype_to_use}.csv")
    return bulk_phase_tests

def analyze_ccd_variation_by_mvavg_rna(adata, wp_ensg, ccd_comp, bioccd, adata_nonccdproteins, adata_regevccdgenes, biotype_to_use):
    expression_data = adata.X # log normalized
    normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
    fucci_time_inds = np.argsort(adata.obs["fucci_time"])
    norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
    moving_averages = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort, WINDOW)
    mvavg_xvals = MovingAverages.mvavg(adata.obs["fucci_time"][fucci_time_inds], WINDOW)
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
        moving_averages_perm = np.apply_along_axis(MovingAverages.mvavg, 1, norm_exp_sort_perm, WINDOW)
        percent_ccd_variance_rng = np.var(moving_averages_perm, axis=1) / np.var(norm_exp_sort_perm, axis=1)
        for iii, perm in enumerate(perms):
            if iii % 50 == 0: print(f"permutation {iii}")
            norm_exp_sort_perm = np.take(normalized_exp_data, perm, axis=0)
            moving_averages_perm = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort_perm, WINDOW)
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

    ccdprotein = np.isin(adata.var_names, np.concatenate((wp_ensg[ccd_comp], bioccd)))
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    gene_ids = list(gene_info["gene_id"])
    gene_names = list(gene_info["name"])
    gene_id_name = dict([(gene_ids[idxx], gene_names[idxx]) for idxx in range(len(gene_info))])
    ccdstring = np.array(["No                 "] * len(ccdprotein))
    ccdstring[np.isin(adata.var_names, wp_ensg[ccd_comp])] = "Pseudotime"
    ccdstring[np.isin(adata.var_names, bioccd)] = "Mitotic"
    ccdstring[np.isin(adata.var_names, wp_ensg[ccd_comp]) & np.isin(adata.var_names, bioccd)] = "Pseudotime&Mitotic"
    percent_variance_tests = pd.DataFrame(
        {"gene" : adata.var_names, 
        "name" : [gene_id_name[x] if x in gene_id_name else "" for x in adata.var_names],
        "ccd_transcript" : pass_meandiff, 
        "regev_ccd" : adata_regevccdgenes,
        "ccd_protein" : ccdstring,
        "nonccd_protein" : adata_nonccdproteins,
        "mean_diff_from_rng":mean_diff_from_rng,
        "-log10 CCD FDR":-np.log10(eq_percvar_adj)})
    percent_variance_tests.to_csv(f"output/transcript_regulation{biotype_to_use}.csv", index=False)

    # And keep track of the ccd genes with and without transcript regulation
    ccdtranscript = pass_meandiff
    ccdprotein_transcript_regulated = ccdprotein & pass_meandiff
    ccdprotein_nontranscript_regulated = ccdprotein & ~pass_meandiff
    ccdtranscript_names = np.array(adata.var_names)[ccdtranscript]
    proteinccd_transcript_regulated_names = np.array(adata.var_names)[ccdprotein_transcript_regulated]
    proteinccd_nontranscript_regulated_names = np.array(adata.var_names)[ccdprotein_nontranscript_regulated]
    utils.np_save_overwriting("output/pickles/ccdprotein.npy", ccdprotein) # pseudotime/mitotic ccd, might not have all the proteins, since this only has proteins not filtered in RNA-Seq analysis
    utils.np_save_overwriting("output/pickles/ccdtranscript.npy", ccdtranscript)
    utils.np_save_overwriting("output/pickles/ccdprotein_transcript_regulated.npy", ccdprotein_transcript_regulated)
    utils.np_save_overwriting("output/pickles/ccdprotein_nontranscript_regulated.npy", ccdprotein_nontranscript_regulated)
    pd.DataFrame({"gene" : ccdtranscript_names}).to_csv("output/all_ccdtranscript_names.csv")
    pd.DataFrame({"gene" : proteinccd_transcript_regulated_names}).to_csv("output/proteinccd_transcript_regulated_names.csv")
    pd.DataFrame({"gene" : proteinccd_nontranscript_regulated_names}).to_csv("output/proteinccd_nontranscript_regulated_names.csv")
    pd.DataFrame({"gene" : adata.var_names}).to_csv("output/gene_names.csv")
    
    # make folders
    folder = f"figures/RNAPseudotimes"
    ccdtransccdprotfolder = f"figures/RNA_CCDTranscriptCCDProtein"
    ccdtransnonccdprotfolder = f"figures/RNA_CCDTranscriptNonCCDProtein"
    nonccdtransccdprotfolder = f"figures/RNA_NonCCDTranscriptCCDProtein"
    nonccdfolder = f"figures/RNA_NonCCD"
    for f in [ccdtransccdprotfolder,ccdtransnonccdprotfolder,nonccdtransccdprotfolder,nonccdfolder]:
        if not os.path.exists(f): os.mkdir(f)
    # CCD transcript & not CCD protein
    for ensg in adata.var_names[ccdtranscript & ~ccdprotein]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdtransnonccdprotfolder, ensg +'_mvavg.pdf'))
    # CCD transcript & CCD Protein
    for ensg in adata.var_names[ccdprotein_transcript_regulated]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(ccdtransccdprotfolder, ensg +'_mvavg.pdf'))
    # Not CCD transcript & CCD Protein
    for ensg in adata.var_names[ccdprotein_nontranscript_regulated]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(nonccdtransccdprotfolder, ensg +'_mvavg.pdf'))
    # Non-CCD 
    for ensg in adata.var_names[~ccdtranscript & ~ccdprotein]:
        shutil.copy(os.path.join(folder, ensg+'_mvavg.pdf'), os.path.join(nonccdfolder, ensg+'_mvavg.pdf'))

    # Figures of merit
    with open("output/figuresofmerit.txt", "a") as file:
        fom = "--- RNA pseudotime\n\n"
        fom += f"We identified {sum(ccdtranscript)} genes of {len(ccdtranscript)} protein-coding genes analyzed ({100 * sum(ccdtranscript) / len(ccdtranscript)}%) to have variance in expression levels correlated to cell cycle progression" + "\n\n"
        fom += f"We can attribute only {100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}% of proteomic cell cycle regulation to transcriptomic cycling with single-cell RNA sequencing" + "\n\n"
        fom += f"This includes {100 * sum(np.isin(adata.var_names[mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM], adata.var_names[adata_regevccdgenes])) / sum(adata_regevccdgenes)}% of known CCD transcripts. Of these, {sum(ccdprotein_transcript_regulated)} were also cell cycle dependent proteins ({100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}%). Of the {sum(ccdprotein)} CCD proteins, {sum(ccdprotein_nontranscript_regulated)} did not have CCD transcripts, including DUSP18 (Figure 2E). There were {sum(ccdtranscript & adata_nonccdproteins)} CCD transcripts that were Non-CCD as proteins." + "\n\n"
        fom += f"" + "\n\n"
        print(fom)
        file.write(fom)
    
    return percent_ccd_variance, total_gini, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms

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
    moving_averages_spike = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort_spike, 100)
    cell_cycle_variance_spike = np.apply_along_axis(np.var, 0, moving_averages_spike)
    total_variance_spike = np.apply_along_axis(np.var, 0, norm_exp_sort_spike)
    total_cv_spike = np.apply_along_axis(scipy.stats.variation, 0, norm_exp_sort_spike)
    percent_ccd_variance_spike = cell_cycle_variance_spike / total_variance_spike
    # avg_expression_spike = np.apply_along_axis(np.median, 0, norm_exp_sort_spike)
    print("Percent variance of spike-in:")
    print(f"mean +/- stdev of spike-in variance explained by cell cycle: {np.mean(percent_ccd_variance_spike)} +/- {np.std(percent_ccd_variance_spike)}")
    print(f"median of spike-in variance explained by cell cycle: {np.median(percent_ccd_variance_spike)}")

    percent_ccd_variance_rng_spike = []
    for iii, perm in enumerate(perms):
        if iii % 1000 == 0: print(f"permutation {iii}")
        norm_exp_sort_perm_spike = np.take(normalized_exp_data_spike, perm, axis=0)
        moving_averages_perm_spike = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort_perm_spike, WINDOW)
        percent_ccd_variance_rng_spike.append(
                np.var(moving_averages_perm_spike, axis=0) / np.var(norm_exp_sort_perm_spike, axis=0))
    percent_ccd_variance_rng_spike = np.asarray(percent_ccd_variance_rng_spike)
    mean_diff_from_rng_spike = np.mean((percent_ccd_variance_spike - percent_ccd_variance_rng_spike).T, 1)
    print("Percent additional variance CCD than random of spike-in")
    print(f"mean +/- stdev of spike-in mean additional percent variance from random: {np.mean(mean_diff_from_rng_spike)} +/- {np.std(mean_diff_from_rng_spike)}")
    print(f"median of spike-in addtional variance explained by cell cycle than random: {np.median(mean_diff_from_rng_spike)}")

    utils.general_boxplot((percent_ccd_variance_spike, mean_diff_from_rng_spike), ("Percent Variance\nCCD Spike-In", "Percent Additional\nCCD Variance Spike-In"), "", "Percent Variance CCD", "", True, "figures/RNASpikeinVarianceBoxplot.png")