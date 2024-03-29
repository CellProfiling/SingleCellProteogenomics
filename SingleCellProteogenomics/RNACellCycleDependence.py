# -*- coding: utf-8 -*-
"""
Methods for assessing cell cycle dependence of RNA abundances in single cells.
-  Percent variance attributed to the cell cycle was calculated using the (variance of moving average / total variance)
-  Randomization analysis was used to determine statistical significance of high percent variances due to the cell cycle

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import (utils, 
                                      MovingAverages, 
                                      FucciCellCycle, 
                                      RNADataPreparation)
from sklearn.linear_model import MultiTaskLassoCV
from sklearn.impute import KNNImputer
import warnings, os, pickle, shutil
import numpy as np
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt
import scanpy as sc
import scipy

np.random.seed(0) # Get the same results each time
WINDOW = 100 # Number of points for moving average window for protein analysis
PERMUTATIONS = 10000 # Number of permutations used for randomization analysis
PERMUTATIONS_ISOFORMS = 100
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
    boxplot.set_xlabel("", size=36)
    boxplot.set_ylabel("Normalized Log10(TPM)", size=18)
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

def plot_expression_facs(adata, genelist, normalized_exp_data, outfolder):
    '''
    Displaying relative expression on the FUCCI plot (log-log intensity of FUCCI markers)
    Output: a folder of FUCCI plots for each gene in a list
    '''
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    var_name_list= list(adata.var_names)
    for gene in genelist:
        nexp = normalized_exp_data[:,var_name_list.index(gene)]
        plt.scatter(adata.obs["Green530"], adata.obs["Red585"], c=nexp)
        plt.tight_layout()
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel("Normalized RNA-Seq Counts", rotation=270,size=16)
        plt.title(gene,size=20)
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
        boxplot.set_xlabel("Cell Cycle Stage", size=36)
        boxplot.set_ylabel("Normalized Transcript Expression", size=36)
        boxplot.tick_params(axis="both", which="major", labelsize=16)
        qqq = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["pvaladj"]
        rej = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["reject"]
        boxplot.set_title(f"{gene} (Q_bh = {utils.format_p(qqq)})",size=36)
        boxplot.get_figure().savefig(f"{outfolder}/{gene}_boxplot.png")
        plt.close()
    pd.DataFrame({"gene" : notfound}).to_csv(f"{outfolder}/GenesFilteredInScRNASeqAnalysis.csv")

def panel_variances_tf(total_var, percent_var, expression_color, title, file_tag):
    '''Displaying the bulk_phase results on the percent variances: plot percent variances from cell line against total variance'''
    plt.scatter(total_var, percent_var, c=expression_color, cmap="bwr_r")
    plt.xlabel("Gini of Gene Expression", size=24)
    plt.ylabel("Percent Variance Due to Cell Cycle",size=24)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlim(-0.005, 1)
    plt.ylim(-0.05, 0.8)
    plt.title(title, size=24)
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
    plt.close()

    # A second plot with just the variance for all
    panel_variances_tf(total_gini, percent_ccd_variance, pass_meandiff, "All Genes", "All")
    plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}_allgenes.pdf")
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

def analyze_ccd_variation_by_mvavg_rna(adata, wp_ensg, ccd_comp, bioccd, adata_nonccdproteins, adata_regevccdgenes, 
               biotype_to_use, use_isoforms=False, make_mvavg_plots_isoforms=False):
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
    perms = np.asarray([np.random.permutation(len(adata.obs)) for nnn in np.arange(PERMUTATIONS if not use_isoforms else PERMUTATIONS_ISOFORMS)])
    picklePath = f"output/pickles/percent_ccd_variance_rng{'' if not use_isoforms else 'Isoforms'}.npy"
    meandiffPath = f"output/pickles/mean_diff_from_rng{'' if not use_isoforms else 'Isoforms'}.npy"
    if not os.path.exists(picklePath):
        # norm_exp_sort_perm = np.asarray([np.take(normalized_exp_data, perm, axis=0) for perm in perms])
        # moving_averages_perm = np.apply_along_axis(MovingAverages.mvavg, 1, norm_exp_sort_perm, WINDOW)
        # percent_ccd_variance_rng = np.var(moving_averages_perm, axis=1) / np.var(norm_exp_sort_perm, axis=1)
        for iii, perm in enumerate(perms):
            if iii % 50 == 0: print(f"permutation {iii}")
            norm_exp_sort_perm = np.take(normalized_exp_data, perm, axis=0)
            moving_averages_perm = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort_perm, WINDOW)
            percent_ccd_variance_rng.append(
                    np.var(moving_averages_perm, axis=0) / np.var(norm_exp_sort_perm, axis=0))
        utils.np_save_overwriting(picklePath, percent_ccd_variance_rng)
    else: 
        percent_ccd_variance_rng = np.load(picklePath, allow_pickle=True)
    percent_ccd_variance_rng = np.asarray(percent_ccd_variance_rng)
    mean_diff_from_rng = np.mean((percent_ccd_variance - percent_ccd_variance_rng).T, 1)
    utils.np_save_overwriting(meandiffPath, mean_diff_from_rng)

    # Statistical testing based on randomization analysis
    alpha_ccd = 0.01
    pass_meandiff = mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM
    ccd_var_comp_rng_wilcoxp = np.apply_along_axis(scipy.stats.wilcoxon, 1, (percent_ccd_variance - percent_ccd_variance_rng).T, None, "wilcox", False, "greater").T[1].T
    eq_percvar_adj, pass_eq_percvar_adj = utils.bonf(alpha_ccd, ccd_var_comp_rng_wilcoxp)
    gtpass_eq_percvar_adj = pass_eq_percvar_adj & (percent_ccd_variance > np.median(percent_ccd_variance_rng, axis=0))

    ccdprotein = np.isin(adata.var_names, np.concatenate((wp_ensg[ccd_comp], bioccd)))
    gene_info = pd.read_csv("input/RNAData/IdsToNames.csv.gz", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
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
    percent_variance_tests.to_csv(f"output/transcript_regulation{biotype_to_use}{'' if not use_isoforms else 'Isoforms'}.csv", index=False)

    # And keep track of the ccd genes with and without transcript regulation
    ccdtranscript = pass_meandiff
    ccdprotein_transcript_regulated = ccdprotein & pass_meandiff
    ccdprotein_nontranscript_regulated = ccdprotein & ~pass_meandiff
    ccdtranscript_names = np.array(adata.var_names)[ccdtranscript]
    proteinccd_transcript_regulated_names = np.array(adata.var_names)[ccdprotein_transcript_regulated]
    proteinccd_nontranscript_regulated_names = np.array(adata.var_names)[ccdprotein_nontranscript_regulated]
    utils.np_save_overwriting(f"output/pickles/ccdprotein{'' if not use_isoforms else 'Isoforms'}.npy", ccdprotein) # pseudotime/mitotic ccd, might not have all the proteins, since this only has proteins not filtered in RNA-Seq analysis
    utils.np_save_overwriting(f"output/pickles/ccdtranscript{'' if not use_isoforms else 'Isoforms'}.npy", ccdtranscript)
    utils.np_save_overwriting(f"output/pickles/ccdprotein_transcript_regulated{'' if not use_isoforms else 'Isoforms'}.npy", ccdprotein_transcript_regulated)
    utils.np_save_overwriting(f"output/pickles/ccdprotein_nontranscript_regulated{'' if not use_isoforms else 'Isoforms'}.npy", ccdprotein_nontranscript_regulated)
    pd.DataFrame({"gene" : ccdtranscript_names}).to_csv(f"output/all_ccdtranscript_names{'' if not use_isoforms else 'Isoforms'}.csv")
    pd.DataFrame({"gene" : proteinccd_transcript_regulated_names}).to_csv(f"output/proteinccd_transcript_regulated_names{'' if not use_isoforms else 'Isoforms'}.csv")
    pd.DataFrame({"gene" : proteinccd_nontranscript_regulated_names}).to_csv(f"output/proteinccd_nontranscript_regulated_names{'' if not use_isoforms else 'Isoforms'}.csv")
    pd.DataFrame({"gene" : adata.var_names}).to_csv(f"output/gene_names{'' if not use_isoforms else 'Isoforms'}.csv")
    
    # make folders
    mvpercs = [] if use_isoforms and not make_mvavg_plots_isoforms else mvavg_plots_pergene(adata, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, use_isoforms)
    if not use_isoforms or make_mvavg_plots_isoforms:
        folder = f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNAPseudotimes{'' if not use_isoforms else 'Isoforms'}"
        ccdtransccdprotfolder = f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNA_CCDTranscriptCCDProtein{'' if not use_isoforms else 'Isoforms'}"
        ccdtransnonccdprotfolder = f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNA_CCDTranscriptNonCCDProtein{'' if not use_isoforms else 'Isoforms'}"
        nonccdtransccdprotfolder = f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNA_NonCCDTranscriptCCDProtein{'' if not use_isoforms else 'Isoforms'}"
        nonccdfolder = f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNA_NonCCD"
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
        fom += f"We identified {sum(ccdtranscript)} {'genes' if use_isoforms else 'transcript isoforms'} of {len(ccdtranscript)} protein-coding {'genes' if use_isoforms else 'transcript isoforms'} analyzed ({100 * sum(ccdtranscript) / len(ccdtranscript)}%) to have variance in expression levels correlated to cell cycle progression" + "\n\n"
        if not use_isoforms:
            fom += f"We can attribute only {100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}% of proteomic cell cycle regulation to transcriptomic cycling with single-cell RNA sequencing" + "\n\n"
            fom += f"This includes {100 * sum(np.isin(adata.var_names[mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM], adata.var_names[adata_regevccdgenes])) / sum(adata_regevccdgenes)}% of known CCD transcripts. Of these, {sum(ccdprotein_transcript_regulated)} were also cell cycle dependent proteins ({100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}%). Of the {sum(ccdprotein)} CCD proteins, {sum(ccdprotein_nontranscript_regulated)} did not have CCD transcripts, including DUSP18 (Figure 2E). There were {sum(ccdtranscript & adata_nonccdproteins)} CCD transcripts that were Non-CCD as proteins." + "\n\n"
        fom += f"" + "\n\n"
        print(fom)
        file.write(fom)
    
    return percent_ccd_variance, total_gini, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms, ccdtranscript, ccdprotein, mvpercs

def compare_genes_to_isoforms(adata, ccdprotein, ccdtranscript, adata_nonccdprotein, adata_isoform, ccdtranscript_isoform):
    '''Check out the isoform results at the gene level'''
    gene_varnames = list(adata.var_names)
    isoform_varnames = list(adata_isoform.var_names)
    isoformToGene = pd.read_csv("input/RNAData/IsoformToGene.csv.gz", index_col=False, header=None, names=["transcript_id", "gene_id"])
    isoformIdList = list(isoformToGene["transcript_id"])
    isoform_varnames_geneids = np.array([isoformToGene["gene_id"][isoformIdList.index(t)] for t in isoform_varnames])
    ccdIsoformWithCcdGene = ccdtranscript_isoform[np.isin(isoform_varnames_geneids, gene_varnames)] & np.array([ccdtranscript[gene_varnames.index(gene_id)] for gene_id in isoform_varnames_geneids if gene_id in gene_varnames])
    numIsoformsPerGene = isoformToGene.groupby("gene_id")["gene_id"].value_counts()
    perGene_geneIds = np.array([x[0] for x in numIsoformsPerGene.index])
    useGene = np.isin(perGene_geneIds, gene_varnames)
    numIsoformsPerGene = np.array(numIsoformsPerGene[useGene])
    ccdIsoformsPerGene = np.array([sum(ccdtranscript_isoform[isoform_varnames_geneids == gene_id]) for gene_id in perGene_geneIds[useGene]])
    ccdIsoformsPerCcdProtein = np.array([sum(ccdtranscript_isoform[isoform_varnames_geneids == gene_id]) for gene_id in adata.var_names[ccdprotein]])
    ccdIsoformsPerNonCcdProtein = np.array([sum(ccdtranscript_isoform[isoform_varnames_geneids == gene_id]) for gene_id in adata.var_names[adata_nonccdprotein]])
    ccdAndNonCcdIsoformsPerGene = np.array([numIsoformsPerGene[ii] != ccdIsoformsPerGene[ii] and ccdIsoformsPerGene[ii] > 0 for ii in np.arange(len(numIsoformsPerGene))])
    print(f"{sum(ccdtranscript_isoform)} CCD transcripts, of which {sum(ccdIsoformWithCcdGene)} ({sum(ccdIsoformWithCcdGene) / sum(ccdtranscript_isoform) * 100}%) correspond to genes that were also found to be CCD.")
    print(f"of the {sum(ccdtranscript)} CCD genes, {sum(ccdAndNonCcdIsoformsPerGene)} were found to have both CCD and non-CCD transcript isoforms.")
    print(f"for the {sum(ccdIsoformsPerGene > 1)} genes with multiple CCD transcripts, the time of peak expression...")
    print(f"{sum(ccdIsoformsPerCcdProtein > 0) / sum(ccdprotein)}% of CCD proteins had at least one CCD transcript isoform")
    print(f"{sum(ccdIsoformsPerNonCcdProtein > 0) / sum(adata_nonccdprotein)}% of non-CCD proteins had at least one CCD transcript isoform")
    
def analyze_isoforms(adata, ccdtranscript, wp_ensg, ccd_comp, nonccd_comp, u_plates, make_mvavg_plots_isoforms=False):
    '''Analyze the isoform-level results over the cell cycle'''
    # Read in the data & QC analysis
    valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
    adata_isoform, phases_isoform = RNADataPreparation.read_counts_and_phases(valuetype, use_spikeins, biotype_to_use, u_plates, use_isoforms=True)
    # RNADataPreparation.plot_pca_for_batch_effect_analysis(adata_isoform, "BeforeRemovingNoncycling_Isoform")
    adata_isoform, phasesfilt_isoform = RNADataPreparation.qc_filtering(adata_isoform, do_log_normalize=True, do_remove_blob=True)
    adata_isoform = RNADataPreparation.zero_center_fucci(adata_isoform)
    # RNADataPreparation.plot_pca_for_batch_effect_analysis(adata_isoform, "AfterRemovingNoncycling_Isoform")
    # FucciPseudotime.pseudotime_umap(adata_isoform, isIsoform=True)
   
    # Cell cycle analysis    
    bioccd = np.genfromtxt("input/ProteinData/BiologicallyDefinedCCD.txt", dtype='str') # from mitotic structures in the protein work
    ccd_regev_filtered_isoform, ccd_filtered_isoform, nonccd_filtered_isoform = utils.ccd_gene_lists(adata_isoform)
    adata_ccdprotein_isoform, adata_nonccdprotein_isoform, adata_regevccdgenes_isoform = RNADataPreparation.is_ccd(adata_isoform, wp_ensg, ccd_comp, nonccd_comp, bioccd, ccd_regev_filtered_isoform)
    rna_ccd_analysis_results = analyze_ccd_variation_by_mvavg_rna(adata_isoform, wp_ensg, ccd_comp, bioccd, adata_nonccdprotein_isoform, adata_regevccdgenes_isoform, biotype_to_use, use_isoforms=True, make_mvavg_plots_isoforms=make_mvavg_plots_isoforms)
    percent_ccd_variance_isoform, total_gini_isoform, mean_diff_from_rng_isoform, pass_meandiff_isoform, eq_percvar_adj_isoform, fucci_time_inds_isoform, norm_exp_sort_isoform, moving_averages_isoform, mvavg_xvals_isoform, perms_isoform, ccdtranscript_isoform, ccdprotein_isoform, mvpercs_isoform = rna_ccd_analysis_results    
    return adata_isoform, ccdtranscript_isoform

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
    plt.close()

    eq_percvar_adj_nextafter = np.nextafter(eq_percvar_adj, eq_percvar_adj + 1)
    plt.figure(figsize=(10,10))
    plt.scatter(mean_diff_from_rng, -np.log10(eq_percvar_adj_nextafter), c=pass_meandiff, cmap="bwr_r")
    plt.vlines(MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM, np.min(-np.log10(eq_percvar_adj_nextafter)), np.max(-np.log10(eq_percvar_adj_nextafter)), color="gray")
    plt.xlabel("Mean Difference from Random")
    plt.ylabel("-log10 adj p-value from randomization")
    plt.savefig("figures/MedianDiffFromRandomVolcano_RNA.png")
    plt.savefig("figures/MedianDiffFromRandomVolcano_RNA.pdf")
    plt.close()

def mvavg_plots_pergene(adata, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, use_isoforms=False):
    '''Generate moving average plots for all the genes'''
    mvpercs = []
    examples = ["ENSG00000011426", "ENSG00000006747", "ENSG00000072756", "ENSG00000102908", 
                "ENSG00000258947", "ENSG00000091651", "ENSG00000169740", "ENSG00000105173", 
                "ENSG00000162999", "ENSG00000134057", "ENSG00000178999", "ENSG00000156970", 
                "ENSG00000167065", "ENSG00000132768", "ENSG00000138801", "ENSG00000156239", 
                "ENSG00000019144", "ENSG00000151702", "ENSG00000123607", "ENSG00000173599", "ENSG00000109814"]
    for iii, ensg in enumerate(adata.var_names):
        plt.close('all')
        if iii % 500 == 0: print(f"well {iii} of {len(adata.var_names)}")  
        windows = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(len(adata.obs["fucci_time"][fucci_time_inds]) - WINDOW + 1)])
        mvperc = MovingAverages.mvpercentiles(norm_exp_sort[:,iii][windows])
        mvpercs.append(mvperc if not use_isoforms else [])
        MovingAverages.temporal_mov_avg_rna(adata.obs["fucci_time"][fucci_time_inds], norm_exp_sort[:,iii], mvavg_xvals, moving_averages[:,iii], mvperc, f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNAPseudotimes{'' if not use_isoforms else 'Isoforms'}", f"{ensg}")
        if ensg in examples:
            print(ensg)
            for plate in np.unique(adata.obs["plate"]):
                isFromPlate = adata.obs["plate"] == plate
                windowPlate = int(WINDOW / len(adata.obs["plate"]) * sum(isFromPlate))
                mvavgPlate = MovingAverages.mvavg(norm_exp_sort[isFromPlate, iii], windowPlate)
                mvavgXPlate = MovingAverages.mvavg(adata.obs["fucci_time"][fucci_time_inds][isFromPlate], windowPlate)
                windowsPlate = np.asarray([np.arange(start, start + windowPlate) for start in np.arange(sum(isFromPlate) - windowPlate + 1)])
                mvpercPlate = MovingAverages.mvpercentiles(norm_exp_sort[isFromPlate, iii][windowsPlate])
                MovingAverages.temporal_mov_avg_rna(adata.obs["fucci_time"][fucci_time_inds][isFromPlate], norm_exp_sort[adata.obs["plate"] == plate,iii], mvavgXPlate, mvavgPlate, mvpercPlate, f"{'f:/CellCycle/' if use_isoforms else ''}figures/RNAPseudotimesExamples", f"{ensg}{plate}")
    return np.array(mvpercs)

def plot_umap_ccd_cutoffs(adata, mean_diff_from_rng):
    '''Make UMAPs with CCD transcripts eliminated at various cutoffs between meandiff > 0.01 and meandiff > 0.1'''
    if not os.path.exists("figures/RNACcdUmapCutoffs"): 
        os.mkdir("figures/RNACcdUmapCutoffs")
    for cutoff in (np.arange(10) + 1) / 100:
        adata_withoutCCd = adata[:,mean_diff_from_rng < cutoff]
        sc.pp.neighbors(adata_withoutCCd, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata_withoutCCd)
        sc.pl.umap(adata_withoutCCd, color="fucci_time", show=False, save=True)
        shutil.move("figures/umap.pdf", f"figures/RNACcdUmapCutoffs/umap{cutoff}cutoff.pdf")

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
 
def compare_to_lasso_analysis(adata, ccdtranscript):
    '''Perform a comparison of pseudotime analysis to LASSO analysis for finding CCD proteins'''
    prevPlotSize = plt.rcParams['figure.figsize']
    plt.rcParams['figure.figsize'] = (6, 5)

    print("ANALYZING SC-RNA-SEQ WITH LASSO")
    warnings.filterwarnings("ignore")
    fucci_rna_data = [(adata.obs["Red585"][ii], adata.obs["Green530"][ii]) for ii in np.arange(len(adata.obs))]
    imputer = KNNImputer(missing_values=0)
    expression = imputer.fit_transform(adata.X)
    fucci_rna_path = "output/pickles/fucci_rna_imputed_lasso.pkl"
    if os.path.exists(fucci_rna_path):
        fucci_rna = np.load(open(fucci_rna_path, 'rb'), allow_pickle=True)
    else:
        fucci_rna = MultiTaskLassoCV()
        fucci_rna.fit(expression, fucci_rna_data)
        pickle.dump(fucci_rna, open(fucci_rna_path, 'wb'))
    nz_coef = np.sum(fucci_rna.coef_, axis=0) != 0
    print(f"{sum(nz_coef)}: number of nonzero lasso coefficients")
    print(f"{adata.var_names[nz_coef]}: genes with nonzero lasso coeff")
    print(f"{sum(ccdtranscript[nz_coef]) / sum(nz_coef)}: % nonzero lasso found as CCD transcripts")
    print(f"{np.sum(fucci_rna.coef_, axis=0)[nz_coef]}: coefficients for nonzero lasso coeff")
    
    # Generate UMAP for CCD and nonCCD for the LASSO model
    adataCCd = adata[:,nz_coef]
    sc.pp.neighbors(adataCCd, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adataCCd)
    sc.pl.umap(adataCCd, color="fucci_time", show=False, save=True)
    shutil.move("figures/umap.pdf", "figures/umapRNALassoCCD.pdf")
    adataNonCCd = adata[:,~nz_coef]
    sc.pp.neighbors(adataNonCCd, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adataNonCCd)
    sc.pl.umap(adataNonCCd, color="fucci_time", show=False, save=True)
    shutil.move("figures/umap.pdf", "figures/umapRNALassoNonCCD.pdf")
    plt.rcParams['figure.figsize'] = prevPlotSize
    warnings.filterwarnings("default")

def analyze_cnv_calls(adata, ccdtranscript):
    '''Take results from cnvkit calls to analyze effects of copy number variation'''
    cnsresults = pd.read_csv("input/RNAData/CnsCallSummary.tsv", sep="\t")
    cnsresults_gene = cnsresults["gene"]
    cnsresults_allgenes = np.concatenate([g.split(',') for g in cnsresults_gene])
    genenamedict = utils.getGeneNameDict()
    adata_names = np.array(utils.ccd_gene_names_gapped(adata.var_names[ccdtranscript], genenamedict))
    adata_ccd_isInCns = adata[np.isin(adata.obs["Well_Plate"], cnsresults.columns), 
                              np.arange(len(ccdtranscript))[ccdtranscript][np.isin(adata_names, cnsresults_allgenes)]]
    adata_ccd_isInCns_names = utils.ccd_gene_names_gapped(adata_ccd_isInCns.var_names, genenamedict)
    cnsresultIdx = np.array([[n in genelist for genelist in cnsresults_gene] for n in adata_ccd_isInCns_names])
    geneInJustOneList = np.array([sum(x) == 1 for x in cnsresultIdx])
    adata_ccd_isInCns_inJustOneList = adata_ccd_isInCns[:, geneInJustOneList]
    adata_ccd_isInCns_inJustOneList_names = utils.ccd_gene_names_gapped(adata_ccd_isInCns_inJustOneList.var_names, genenamedict)
    cnsresultIdx_inJustOneList = cnsresultIdx[geneInJustOneList]
    cnsResultsCellData = np.array(cnsresults)[:, np.isin(cnsresults.columns, adata_ccd_isInCns_inJustOneList.obs["Well_Plate"])]
    
    # evaluate consistency of CNVs
    heatmap = np.zeros(cnsResultsCellData.T.shape)
    heatmap[cnsResultsCellData.T == -5] = -1
    heatmap[(cnsResultsCellData.T > -5) & (cnsResultsCellData.T < 1)] = 0
    heatmap[cnsResultsCellData.T == 1] = 1
    heatmap[cnsResultsCellData.T == 2] = 2
    heatmap[cnsResultsCellData.T > 2] = 3
    clustergrid = sbn.clustermap(heatmap[:,:-8], col_cluster=False)
    plt.savefig("figures/CnvConsistency.pdf")
    plt.close()
    
    # heatmaps for phases
    adata_idx = np.array([list(adata.obs["Well_Plate"]).index(wp) for wp in cnsresults.columns[np.isin(cnsresults.columns, 
                                               adata_ccd_isInCns_inJustOneList.obs["Well_Plate"])]])
    sbn.heatmap([adata_ccd_isInCns.obs["phase"][np.asarray(clustergrid.dendrogram_row.reordered_ind)] == "G1",
                 adata_ccd_isInCns.obs["phase"][np.asarray(clustergrid.dendrogram_row.reordered_ind)] == "S-ph",
                 adata_ccd_isInCns.obs["phase"][np.asarray(clustergrid.dendrogram_row.reordered_ind)] == "G2M"],
                yticklabels=["G1", "S", "G2"])
    plt.savefig("figures/CnvConsistencyPhases.pdf")
    plt.close()
    
    # is there enrichment for phase in the highly amplified genes?
    # print(adata_ccd_isInCns.obs["phase"][clustergrid.dendrogram_row.reordered_ind[:100]].value_counts())
    
    # yes, so is there correlation?
    x = adata_ccd_isInCns.obs["fucci_time"]
    y = np.mean(cnsResultsCellData, axis=0)
    linearModel = scipy.stats.linregress(np.asarray(x).astype(float), np.asarray(y).astype(float))
    plt.scatter(x * fucci.TOT_LEN, y)
    plt.scatter(x * fucci.TOT_LEN, linearModel.intercept + x * linearModel.slope)
    plt.xlabel("Cell Division Time, hrs")
    plt.ylabel("Mean CNV of All Chromosome Arms")
    plt.savefig("figures/CnvCorrelation.pdf")
    plt.close()
    
    print(f"{linearModel[3]}: p-value for nonzero slope by two-sided t test")
    residualLinearModel = scipy.stats.linregress(np.asarray(x).astype(float), np.asarray(y - (linearModel.intercept + x * linearModel.slope)).astype(float))
    residualNormality = scipy.stats.normaltest(np.asarray(y - (linearModel.intercept + x * linearModel.slope)))
    print(f"{residualLinearModel[3]}: p-value for nonzero slope of residuals by two-sided t-test")
    print(f"{residualNormality[1]}: p-value for normality of residuals")
    
    # what if we only look at one phase? G1 before doubling? for all genes?
    adata_names = np.array(utils.ccd_gene_names_gapped(adata.var_names, genenamedict))
    adata_ccd_isInCns = adata[np.isin(adata.obs["Well_Plate"], cnsresults.columns) & (adata.obs["phase"] == "G1"), np.arange(len(adata_names))[np.isin(adata_names, cnsresults_allgenes)]]
    adata_ccd_isInCns_names = utils.ccd_gene_names_gapped(adata_ccd_isInCns.var_names, genenamedict)
    cnsresultIdx = np.array([[n in genelist for genelist in cnsresults_gene] for n in adata_ccd_isInCns_names])
    geneInJustOneList = np.array([sum(x) == 1 for x in cnsresultIdx])
    adata_ccd_isInCns_inJustOneList = adata_ccd_isInCns[:, geneInJustOneList]
    adata_ccd_isInCns_inJustOneList_names = utils.ccd_gene_names_gapped(adata_ccd_isInCns_inJustOneList.var_names, genenamedict)
    cnsresultIdx_inJustOneList = cnsresultIdx[geneInJustOneList]
    cnsResultsCellData = np.array(cnsresults)[:, np.isin(cnsresults.columns, adata_ccd_isInCns_inJustOneList.obs["Well_Plate"])]
    cnvAmplified, cnvPvalOneSided = [],[]
    cnvDeleted, cnvPvalOneSidedDeleted = [],[]
    amplifiedTpmsAll, neutralTpmsAll, deletionTpmsAll = [],[],[]
    for ii, tpm in enumerate(adata_ccd_isInCns.X.T[geneInJustOneList]):
        cnv = np.concatenate(cnsResultsCellData[cnsresultIdx_inJustOneList[ii],:])
        missingData = cnv == -5
        amplified, amplifiedTpms = cnv[~missingData & (cnv > 1)], tpm[~missingData & (cnv > 1)]
        neutral, neutralTpms = cnv[~missingData & (cnv == 1)], tpm[~missingData & (cnv == 1)]
        deletion, deletionTpms = cnv[~missingData & (cnv < 1)], tpm[~missingData & (cnv < 1)]
        cnvAmplified.append(np.median(amplifiedTpms) > np.median(tpm[~missingData]))
        cnvPvalOneSided.append(scipy.stats.kruskal(amplifiedTpms, neutralTpms)[1] * 2)
        cnvDeleted.append(np.median(deletionTpms) < np.median(tpm[~missingData]))
        cnvPvalOneSidedDeleted.append(scipy.stats.kruskal(deletionTpms, neutralTpms)[1] * 2)
        amplifiedTpmsAll.extend(amplifiedTpms)
        neutralTpmsAll.extend(neutralTpms)
        deletionTpmsAll.extend(deletionTpms)
    cnvAmplified = np.asarray(cnvAmplified)
    cnvTestPvals_BH, cnvTestPvals_rejectBH = utils.benji_hoch(0.01, cnvPvalOneSided)
    cnvTestPvalsDel_BH, cnvTestPvalsDel_rejectBH = utils.benji_hoch(0.01, cnvPvalOneSidedDeleted)
    print(f"{sum(cnvAmplified & cnvTestPvals_rejectBH)}: number of novel CCD with significantly higher expression with amplified CNVs than neutral")
    print(f"{sum(cnvDeleted & cnvTestPvalsDel_rejectBH)}: number of novel CCD with significantly higher expression with amplified CNVs than neutral")
    utils.general_boxplot([amplifiedTpmsAll, neutralTpmsAll, deletionTpmsAll], 
                          ["amplified", "neutral", "deletion"], "", "logTPMs", "", False, "figures/CNVStateBoxplot.pdf")
    print(f"Of {len(cnvAmplified)} genes:")
    print(f"{scipy.stats.kruskal(amplifiedTpmsAll, neutralTpmsAll, deletionTpmsAll)[1]}: kruskal two sided pval that there's a difference between the three")
    print(f"{scipy.stats.kruskal(amplifiedTpmsAll, neutralTpmsAll)[1]}: kruskal two sided pval that there's a difference between amplified/neutral")


def cell_data_string(adata, idx):
    '''Make a string representing the data for each cell'''
    cell_data = [adata.obs['fucci_time'][idx], 
                 adata.obs['phase'][idx], 
                 adata.obs['Red585'][idx],
                 adata.obs['Green530'][idx]]
    return "|".join([str(xx) for xx in cell_data])
    
def make_plotting_dataframe(adata, ccdtranscript, wp_ensg, bioccd, norm_exp_sort, mvavgs_x, moving_averages, mvpercs):
    '''Make a single table for HPA website figures on RNA pseudotime, boxplots, and fucci plots'''
    hasProteinData = np.isin(adata.var_names, np.concatenate((wp_ensg, bioccd)))
    filterccd = ccdtranscript | hasProteinData
    pd.DataFrame({
        "ENSG" : adata.var_names[filterccd],
        "CCD" : ccdtranscript[filterccd],
        "cell_pseudotime" : [",".join([str(xx) for xx in adata.obs['fucci_time']]) for ii in np.arange(sum(filterccd))],
        "cell_intensity" :  [",".join([str(yyy) for yyy in yy]) for ii, yy in enumerate(norm_exp_sort.T[filterccd])],
        "cell_fred" : [",".join([str(xx) for xx in adata.obs['Red585']]) for ii in np.arange(sum(filterccd))],
        "cell_fgreen" : [",".join([str(xx) for xx in adata.obs['Green530']]) for ii in np.arange(sum(filterccd))],
        "mvavg_x" : [",".join([str(xx) for xx in mvavgs_x]) for ii in np.arange(sum(filterccd))],
        "mvavg_y" : [",".join([str(yyy) for yyy in yy]) for yy in moving_averages.T[filterccd]],
        "mvavgs_10p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[filterccd,0,:]],
        "mvavgs_90p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[filterccd,-1,:]],
        "mvavgs_25p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[filterccd,1,:]],
        "mvavgs_75p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[filterccd,-2,:]],
        "phase" : [",".join([str(xx) for xx in adata.obs['phase']]) for ii in np.arange(sum(filterccd))]
        }).to_csv("output/RNAPseudotimePlotting.csv.gz", index=False, sep="\t")
    pd.DataFrame({
        "ENSG" : adata.var_names,
        "CCD" : ccdtranscript,
        "cell_pseudotime" : [",".join([str(xx) for xx in adata.obs['fucci_time']]) for ii in np.arange(len(adata.var_names))],
        "cell_intensity" :  [",".join([str(yyy) for yyy in yy]) for ii, yy in enumerate(norm_exp_sort.T)],
        "cell_fred" : [",".join([str(xx) for xx in adata.obs['Red585']]) for ii in np.arange(len(adata.var_names))],
        "cell_fgreen" : [",".join([str(xx) for xx in adata.obs['Green530']]) for ii in np.arange(len(adata.var_names))],
        "mvavg_x" : [",".join([str(xx) for xx in mvavgs_x]) for ii in np.arange(len(adata.var_names))],
        "mvavg_y" : [",".join([str(yyy) for yyy in yy]) for yy in moving_averages.T],
        "mvavgs_10p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[:,0,:]],
        "mvavgs_90p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[:,-1,:]],
        "mvavgs_25p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[:,1,:]],
        "mvavgs_75p" : [",".join([str(yyy) for yyy in yy]) for yy in mvpercs[:,-2,:]],
        "phase" : [",".join([str(xx) for xx in adata.obs['phase']]) for ii in np.arange(len(adata.var_names))]
        }).to_csv("output/RNAPseudotimePlotting_Unfiltered.csv.gz", index=False, sep="\t")