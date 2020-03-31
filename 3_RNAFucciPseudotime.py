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
adata_regevccdgenes = np.isin(adata.var_names, ccd_regev_filtered)

# Generate plots with expression of genes overlayed
do_make_gene_expression_plots = False
expression_data = adata.X
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
if do_make_gene_expression_plots:
    # UMAPs with RNA expression overlayed
    RNACellCycleDependence.plot_expression_umap(adata, ccd_regev_filtered, "figures/RegevGeneExpressionUmap") # curated CCD genes from a different scRNA-Seq analysis
    RNACellCycleDependence.plot_expression_umap(adata, ccd_filtered, "figures/CcdGeneExpressionUmap") # CCD proteins from the present immunofluorescense work
    RNACellCycleDependence.plot_expression_umap(adata, nonccd_filtered, "figures/NonCcdGeneExpressionUmap") # non-CCD proteins from the present immunofluorescense work

    # Log-log FUCCI plot with RNA expression overlayed
    RNACellCycleDependence.plot_expression_facs(ccd_regev_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/RegevGeneFucci")
    RNACellCycleDependence.plot_expression_facs(ccd_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/CcdGeneFucci")
    RNACellCycleDependence.plot_expression_facs(nonccd_filtered, normalized_exp_data, phasesfilt, adata.var_names, "figures/NonCcdGeneFucci")

#%% Cluster the expression into phases and analyze it that way
bulk_phase_tests = RNACellCycleDependence.analyze_ccd_variation_by_phase_rna(adata, normalized_exp_data, adata_regevccdgenes, biotype_to_use)
if do_make_gene_expression_plots:
     # Remove?
    RNACellCycleDependence.plot_expression_boxplots(adata, ccd_regev_filtered, bulk_phase_tests, "figures/RegevGeneBoxplots")
    RNACellCycleDependence.plot_expression_boxplots(adata, ccd_filtered, bulk_phase_tests, "figures/DianaCcdGeneBoxplots")
    RNACellCycleDependence.plot_expression_boxplots(adata, nonccd_filtered, bulk_phase_tests, "figures/DianaNonCcdGeneBoxplots")

#%% Moving average calculations and randomization analysis for RNA
adata_ccdprotein, adata_nonccdprotein, adata_regevccdgenes = RNADataPreparation.is_ccd(adata, wp_ensg, ccd_comp, nonccd_comp, bioccd, ccd_regev_filtered)

rna_ccd_analysis_results = RNACellCycleDependence.analyze_ccd_variation_by_mvavg_rna(adata, wp_ensg, ccd_comp, bioccd, adata_nonccdprotein, adata_regevccdgenes, biotype_to_use)
percent_ccd_variance, total_gini, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals, perms = rna_ccd_analysis_results

RNACellCycleDependence.figures_ccd_analysis_rna(adata, percent_ccd_variance, mean_diff_from_rng, pass_meandiff, eq_percvar_adj, wp_ensg, ccd_comp, ccd_regev_filtered)
RNACellCycleDependence.mvavg_plots_pergene(adata, fucci_time_inds, norm_exp_sort, moving_averages, mvavg_xvals)
RNACellCycleDependence.plot_overall_and_ccd_variances(adata, biotype_to_use, total_gini, percent_ccd_variance, pass_meandiff, adata_ccdprotein, adata_nonccdprotein, adata_regevccdgenes)

# (Remove?) What are these low percent variance ones that are significant?
print("What are these low percent variance ones that are significant?")
low_variance_signif = (percent_ccd_variance < 0.03) & ccd_comp & pass_meandiff
print(np.array(adata.var_names)[low_variance_signif])


#%% Moving average calculations and randomization analysis for the spike-in internal controls
adata_spikeins, phases_spikeins = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spike_ins=True, biotype_to_use="")
sc.pp.filter_genes(adata_spikeins, min_cells=100)
print(f"data shape after filtering: {adata_spikeins.X.shape}")

RNACellCycleDependence.ccd_analysis_of_spikeins(adata_spikeins, perms)