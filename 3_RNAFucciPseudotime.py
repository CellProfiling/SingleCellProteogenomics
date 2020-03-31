#%% imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, stretch_time, FucciCellCycle, FucciPseudotime, RNADataPreparation, RNACellCycleDependence
from scipy.optimize import least_squares
import decimal
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)

#%% Convert FACS intensities for FUCCI markers to pseudotime using the same polar coordinate methods as for protein
# Idea: Use the polar coordinate pseudotime calculations to calculate the pseudotime for each cell
# Execution: Adapt Devin's code for the cells sorted for RNA-Seq
# Output: Make log-log fucci intensity plots for the cells analyzed by RNA-Seq; Plot of all fucci pseudotimes; table of pseudotimes for each cell
adata, phases_filt = RNADataPreparation.read_counts_and_phases("All", "Counts", False, "protein_coding") # no qc, yet
FucciPseudotime.pseudotime_rna(adata, phases_filt)

#%% Single cell RNA-Seq data preparation and general analysis
RNADataPreparation.general_plots()
RNADataPreparation.analyze_noncycling_cells()

#%% Idea: Similar to mock-bulk analysis for proteins, we can evaluate each gene bundled by phase across cells
# Execution: Make boxplots of RNA expression by phase
# Output: boxplots for each gene
plate, valuetype, use_spikeins, biotype_to_use = "All", "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spikeins, biotype_to_use)
adata, phasesfilt = RNADataPreparation.qc_filtering(adata, do_log_normalize=True, do_remove_blob=True)
g1, s, g2 = adata.obs["phase"] == "G1", adata.obs["phase"] == "S-ph", adata.obs["phase"] == "G2M"
do_make_boxplots = False
if do_make_boxplots:
    for iii, ensg in enumerate(adata.var_names):
        maxtpm = np.max(np.concatenate((adata.X[g1,iii], adata.X[s,iii], adata.X[g2,iii])))
        RNACellCycleDependence.boxplot_result(adata.X[g1,iii] / maxtpm, adata.X[s,iii] / maxtpm, adata.X[g2,iii] / maxtpm, "figures/RNABoxplotByPhase", ensg)

#%% Idea: Display general RNA expression patterns in single cells using UMAP dimensionality reduction, and display with FUCCI pseudotime overlayed
FucciPseudotime.pseudotime_umap(adata) # Generate a UMAP with the pseudotime overlayed

# We can also show that the cycle pattern disappears when the curated CCD genes or CCD proteins are removed,
# demonstrating that there's valuable information about cycling in these datasets
RNADataPreparation.demonstrate_loss_of_umap_cycle(adata)

# Read in the currated CCD genes / CCD proteins from the present work / Non-CCD genes from the present work; filter for genes that weren't filtered in QC of RNA-Seq
ccd_regev_filtered, ccd_filtered, nonccd_filtered = RNADataPreparation.ccd_gene_lists(adata)

# Generate plots with expression of genes overlayed
do_make_gene_expression_plots = False
expression_data = adata.X
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
if do_make_gene_expression_plots:
    # UMAPs with RNA expression overlayed
    RNACellCycleDependence.plot_expression_umap(ccd_regev_filtered, "figures/RegevGeneExpressionUmap") # curated CCD genes from a different scRNA-Seq analysis
    RNACellCycleDependence.plot_expression_umap(ccd_filtered, "figures/CcdGeneExpressionUmap") # CCD proteins from the present immunofluorescense work
    RNACellCycleDependence.plot_expression_umap(nonccd_filtered, "figures/NonCcdGeneExpressionUmap") # non-CCD proteins from the present immunofluorescense work

    # Log-log FUCCI plot with RNA expression overlayed
    RNACellCycleDependence.plot_expression_facs(ccd_regev_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/RegevGeneFucci")
    RNACellCycleDependence.plot_expression_facs(ccd_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/CcdGeneFucci")
    RNACellCycleDependence.plot_expression_facs(nonccd_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/NonCcdGeneFucci")


#%% Cluster the expression into phases and analyze it that way
stages = np.array(adata.obs["phase"])
g1_exp = np.take(normalized_exp_data, np.nonzero(stages == "G1")[0], axis=0)
s_exp = np.take(normalized_exp_data, np.nonzero(stages == "S-ph")[0], axis=0)
g2_exp = np.take(normalized_exp_data, np.nonzero(stages == "G2M")[0], axis=0)
tests_fp = [scipy.stats.kruskal(g1_exp[:,geneidx], s_exp[:,geneidx], g2_exp[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
pvals = [p for (F, p) in tests_fp]
pvals_corrected_BH, reject_BH = utils.benji_hoch(0.01, pvals)
pvals_correctedBonf, rejectBonf = utils.bonf(0.01, pvals)
bulk_phase_tests = pd.DataFrame(
    {"gene" : adata.var_names,
    "pvalue" : pvals,
    "pvaladj_BH" : pvals_corrected_BH,
    "reject_BH" : reject_BH,
    "pvaladj_B" : pvals_correctedBonf,
    "reject_B" : rejectBonf})
bulk_phase_tests.to_csv(f"output/transcript_regulation{biotype_to_use}.csv")

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
        qqq = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["pvaladj"]
        rej = bulk_phase_tests[bulk_phase_tests.gene == gene].iloc[0]["reject"]
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

#% Plotting variances of gene expression
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
 
#%% Moving average calculations and randomization analysis for RNA
rna_ccd_analysis_results = RNACellCycleDependence.analyze_ccd_variation_rna(adata)
percent_ccd_variance, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms = rna_ccd_analysis_results
RNACellCycleDependence.figures_ccd_analysis_rna(adata, percent_ccd_variance, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, wp_ensg, ccd_comp, ccd_regev_filtered)
RNACellCycleDependence.mvavg_plots_pergene(adata, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals)

#%% Moving average calculations and randomization analysis for the spike-in internal controls
adata_spikeins, phases_spikeins = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spike_ins=True, biotype_to_use="")
sc.pp.filter_genes(adata_spikeins, min_cells=100)
print(f"data shape after filtering: {adata_spikeins.X.shape}")

RNACellCycleDependence.ccd_analysis_of_spikeins(adata_spikeins, perms)


#%% Variances and outputs
def plot_variances(total_var, percent_var, expression_color, title, file_tag):
    '''Plots percent variances from cell line against total variance'''
    plt.figure(figsize=(10,10))
    plt.scatter(total_var, percent_var, c=pass_meandiff, cmap="bwr_r")
    plt.xlabel("Gini of Gene Expression", size=36,fontname='Arial')
    plt.ylabel("Percent Variance Due to Cell Cycle",size=36,fontname='Arial')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.title(title, size=36, fontname='Arial')
#    cbar = plt.colorbar()
#    cbar.ax.get_yaxis().labelpad = 24
#    cbar.ax.set_ylabel("Median Log-Normalized RNA-Seq Counts", rotation=270,size=24,fontname='Arial')
    plt.savefig(f"figures/VariancePlot{file_tag}_{biotype_to_use}.png")
    plt.show()
    plt.close()

regevccdgenes = np.isin(adata.var_names, ccd_regev_filtered)
dianaccdgenes = np.isin(adata.var_names, ccd_filtered)
diananonccdgenes = np.isin(adata.var_names, nonccd_filtered)
plot_variances(total_gini, percent_ccd_variance, avg_expression, "All Genes", "All")
#plot_variances(total_variance[regevccdgenes], percent_ccd_variance[regevccdgenes], avg_expression[regevccdgenes], "Regev CCD Genes", "RegevCcd")
#plot_variances(total_variance[dianaccdgenes], percent_ccd_variance[dianaccdgenes], avg_expression[dianaccdgenes], "Diana CCD Genes", "DianaCcd")
#plot_variances(total_variance[diananonccdgenes], percent_ccd_variance[diananonccdgenes], avg_expression[diananonccdgenes], "Diana Non-CCD Genes", "DianaNonCCD")

# Displaying the bulk_phase results on the percent variances
def plot_variances_tf(total_var, percent_var, expression_color, title, file_tag, percent_var_cutoff, total_var_cutoff):
    '''Plots percent variances from cell line against total variance'''
    # plt.figure(figsize=(10,10))
    plt.scatter(total_var, percent_var, c=expression_color, cmap="bwr_r")
#    plt.axhline(percent_var_cutoff)
#    plt.axvline(total_var_cutoff)
    plt.xlabel("Gini of Gene Expression", size=24,fontname='Arial')
    plt.ylabel("Percent Variance Due to Cell Cycle",size=24,fontname='Arial')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlim(-0.005, 1)
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

# Okay, we're going to take the ones that passed bulk_phase and have > 10% variance explained
total_var_cutoff = np.mean(total_variance_spike) + 1 * np.std(total_variance_spike)
percent_var_cutoff = np.mean(percent_ccd_variance_spike) + 1 * np.std(percent_ccd_variance_spike)
print(f"using total variance cutoff {total_var_cutoff} and percent CCD variance cutoff {percent_var_cutoff}")
plt.figure(figsize=(15,15))
plt.subplot(221)
plot_variances_tf(total_gini, percent_ccd_variance, pass_meandiff, "All Genes", "All", percent_var_cutoff, total_var_cutoff)
plt.subplot(222)
plot_variances_tf(total_gini[regevccdgenes], percent_ccd_variance[regevccdgenes], pass_meandiff[regevccdgenes], "Regev CCD Genes", "RegevCcd", percent_var_cutoff, total_var_cutoff)
plt.subplot(223)
ccdprotein = np.isin(adata.var_names, np.concatenate((wp_ensg[ccd_comp], bioccd)))
plot_variances_tf(total_gini[ccdprotein], percent_ccd_variance[ccdprotein], pass_meandiff[ccdprotein], "CCD Proteins", "DianaCcd", percent_var_cutoff, total_var_cutoff)
plt.subplot(224)
nonccdprotein = np.isin(adata.var_names, wp_ensg[nonccd_comp]) & ~np.isin(adata.var_names, bioccd)
plot_variances_tf(total_gini[nonccdprotein], percent_ccd_variance[nonccdprotein], pass_meandiff[nonccdprotein], "Non-CCD Proteins", "DianaNonCCD", percent_var_cutoff, total_var_cutoff)
plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}.pdf")
plt.show()
plt.close()

plot_variances_tf(total_gini, percent_ccd_variance, pass_meandiff, "All Genes", "All", percent_var_cutoff, total_var_cutoff)
plt.savefig(f"figures/VarianceSignificancePlots{biotype_to_use}_allgenes.pdf")
plt.show()
plt.close()

# Let's output those
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
    "regev_ccd" : regevccdgenes,
    "ccd_protein" : ccdstring,
    "nonccd_protein" : nonccdprotein,
    "mean_diff_from_rng":mean_diff_from_rng,
    "-log10 CCD FDR":-np.log10(pvals_corrected_)})
percent_variance_tests.to_csv(f"output/transcript_regulation{biotype_to_use}.csv", index=False)

# And output the CCD genes that aren't in diana's CCD gene set
# percent_variance_tests[percent_variance_tests.significant & ~percent_variance_tests.ccd_protein].to_csv("output/transcript_regulation_significant_ccd_notindianasset.csv")

# And keep track of the ccd genes with and without transcript regulation
ccdtranscript = pass_meandiff
ccdprotein_transcript_regulated = ccdprotein & pass_meandiff
ccdprotein_nontranscript_regulated = ccdprotein & ~pass_meandiff
ccdtranscript_names = np.array(adata.var_names)[ccdtranscript]
proteinccd_transcript_regulated_names = np.array(adata.var_names)[ccdprotein_transcript_regulated]
proteinccd_nontranscript_regulated_names = np.array(adata.var_names)[ccdprotein_nontranscript_regulated]
np_save_overwriting("output/pickles/ccdprotein.npy", ccdprotein) # pseudotime/mitotic ccd, might not have all the proteins, since this only has proteins not filtered in RNA-Seq analysis
np_save_overwriting("output/pickles/ccdtranscript.npy", ccdtranscript)
np_save_overwriting("output/pickles/ccdprotein_transcript_regulated.npy", ccdprotein_transcript_regulated)
np_save_overwriting("output/pickles/ccdprotein_nontranscript_regulated.npy", ccdprotein_nontranscript_regulated)
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

#% Figures of merit
with open("output/figuresofmerit.txt", "a") as file:
    fom = "--- RNA pseudotime\n\n"
    fom += f"We identified {sum(ccdtranscript)} genes of {len(ccdtranscript)} protein-coding genes analyzed ({100 * sum(ccdtranscript) / len(ccdtranscript)}%) to have variance in expression levels correlated to cell cycle progression" + "\n\n"
    fom += f"We can attribute only {100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}% of proteomic cell cycle regulation to transcriptomic cycling with single-cell RNA sequencing" + "\n\n"
    fom += f"This includes {100 * sum(np.isin(adata.var_names[mean_diff_from_rng > MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM], ccd_regev_filtered)) / len(ccd_regev_filtered)}% of known CCD transcripts. Of these, {sum(ccdprotein_transcript_regulated)} were also cell cycle dependent proteins ({100 * sum(ccdprotein_transcript_regulated) / sum(ccdprotein)}%). Of the {sum(ccdprotein)} CCD proteins, {sum(ccdprotein_nontranscript_regulated)} did not have CCD transcripts, including DUSP18 (Figure 2E). There were {sum(ccdtranscript & nonccdprotein)} CCD transcripts that were Non-CCD as proteins." + "\n\n"
    fom += f"" + "\n\n"
    print(fom)
    file.write(fom)
