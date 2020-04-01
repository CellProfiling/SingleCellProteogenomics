#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, Loaders, RNADataPreparation, PTMAnalysis, ProteinStabilityAnalysis
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

#%% Import the genes names we're analyzing
# Read in RNA-Seq data again and the CCD gene lists
plate, valuetype, use_spikeins, biotype_to_use = "All", "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spikeins, biotype_to_use)
adata, phasesfilt = RNADataPreparation.qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)

import_dict = Loaders.load_ptm_and_stability(adata)
wp_max_pol = import_dict["wp_max_pol"]
ccdtranscript, nonccdtranscript = import_dict["ccdtranscript_names"], import_dict["nonccdtranscript_names"]
ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated = import_dict["ccdprotein_transcript_regulated_names"], import_dict["ccdprotein_nontranscript_regulated_names"]
genes_analyzed = import_dict["genes_analyzed_names"]
ccd_regev_filtered, ccd_filtered = import_dict["ccd_regev_filtered_names"], import_dict["ccd_filtered_names"] 
ccdprotein, nonccdprotein = import_dict["ccdprotein_names"], import_dict["nonccdprotein_names"]
utils.save_gene_names_by_category(adata) # We're using gene names for this analysis instead of ENSG gene identifiers, so save those gene names

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in MetaMorpheus results
genemodsBulk = PTMAnalysis.analyze_ptms("input/raw/U2OSBulkAllProteinGroups.tsv", ccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, genes_analyzed)
genemodsPhospho = PTMAnalysis.analyze_ptms("input/raw/U2OSPhosphoAllProteinGroups.tsv", ccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, genes_analyzed)

# Analyze the modification occupancy for each PTM site
dfBulk, occdfBulk = PTMAnalysis.process_genemods(genemodsBulk)
dfPhospho, occdfPhospho = PTMAnalysis.process_genemods(genemodsPhospho)

# Compare the PTM regulation between transcript and non-transcript regulated groups
PTMAnalysis.compare_ptm_regulation(dfBulk, occdfBulk, "Bulk", genes_analyzed, ccd_regev_filtered, ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, ccdprotein, nonccdprotein)
PTMAnalysis.compare_ptm_regulation(dfPhospho, occdfPhospho, "Phospho", genes_analyzed, ccd_regev_filtered, ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, ccdprotein, nonccdprotein)

#%% Look at time of peak expression for phosphosites on proteins that are CCD
# Conclusion: there's not much going on in terms of correlation to the peak expression of the protein
hngcWithGaps = list(geneIdToHngc_withgaps(wp_ensg))
ccd_modctsdf, ccd_modcts, ccd_avgocc, ccd_modctsoccdf, ccd_modocc = count_mods(ccdprotein, dfPhospho, occdfPhospho)
maxpol_per_gene = np.array([wp_max_pol[hngcWithGaps.index(gene)] if gene in hngcWithGaps else -1 for gene in ccd_modctsoccdf["gene"]])

plt.scatter(maxpol_per_gene[maxpol_per_gene >= 0], ccd_modocc[maxpol_per_gene >= 0])
plt.close()

phase = np.array(["g1" if pol * fucci.TOT_LEN < fucci.G1_LEN else "g1s" if pol * fucci.TOT_LEN >= fucci.G1_LEN and pol * fucci.TOT_LEN < fucci.G1_LEN + fucci.G1_S_TRANS else "g2" for pol in maxpol_per_gene if pol >= 0])
g1 = ccd_modocc[maxpol_per_gene >= 0][phase == "g1"]
g1s = ccd_modocc[maxpol_per_gene >= 0][phase == "g1s"]
g2 = ccd_modocc[maxpol_per_gene >= 0][phase == "g2"]

utils.boxplot_with_stripplot((g1, g1s, g2), ("G1", "G1S", "G2"),
        "", "Occupancy per Modified Residue", "", False, f"figures/ModsOccupancyBoxplotSelected{analysisTitle}.pdf",alpha=0.2, size=7, jitter=0.25, ylim=(-0.01,1.01))


#%% Finessing
# Idea: What are these modifications that have high occupancies?
# Exec: pandas
# Output: table of modifications and the counts at different cutoffs
print(all_modctsoccdf[all_modctsoccdf.occupancy == 1]["modification"].value_counts())
print(ccd_rt_modctsoccdf[ccd_rt_modctsoccdf.occupancy == 1]["modification"].value_counts())
print(ccd_at_modctsoccdf[ccd_at_modctsoccdf.occupancy == 1]["modification"].value_counts())
print(ccd_t_modctsoccdf[ccd_t_modctsoccdf.occupancy == 1]["modification"].value_counts())
print(ccd_n_modctsoccdf[ccd_n_modctsoccdf.occupancy == 1]["modification"].value_counts())

#%%
print(all_modctsoccdf[all_modctsoccdf.occupancy >= 0.5]["modification"].value_counts())
print(ccd_rt_modctsoccdf[ccd_rt_modctsoccdf.occupancy >= 0.5]["modification"].value_counts())
print(ccd_at_modctsoccdf[ccd_at_modctsoccdf.occupancy >= 0.5]["modification"].value_counts())
print(ccd_t_modctsoccdf[ccd_t_modctsoccdf.occupancy >= 0.5]["modification"].value_counts())
print(ccd_n_modctsoccdf[ccd_n_modctsoccdf.occupancy >= 0.5]["modification"].value_counts())

#%%
print(all_modctsoccdf[all_modctsoccdf.occupancy >= 0.5]["modification"].value_counts() / all_modctsoccdf["modification"].value_counts())
print(ccd_rt_modctsoccdf[ccd_rt_modctsoccdf.occupancy >= 0.5]["modification"].value_counts() / ccd_rt_modctsoccdf["modification"].value_counts())
print(ccd_at_modctsoccdf[ccd_at_modctsoccdf.occupancy >= 0.5]["modification"].value_counts() / ccd_at_modctsoccdf["modification"].value_counts())
print(ccd_t_modctsoccdf[ccd_t_modctsoccdf.occupancy >= 0.5]["modification"].value_counts() / ccd_t_modctsoccdf["modification"].value_counts())
print(ccd_n_modctsoccdf[ccd_n_modctsoccdf.occupancy >= 0.5]["modification"].value_counts() / ccd_n_modctsoccdf["modification"].value_counts())

#%%
print(all_modctsoccdf[all_modctsoccdf.occupancy >= 0]["modification"].value_counts())
print(ccd_rt_modctsoccdf[ccd_rt_modctsoccdf.occupancy >= 0]["modification"].value_counts())
print(ccd_at_modctsoccdf[ccd_at_modctsoccdf.occupancy >= 0]["modification"].value_counts())
print(ccd_t_modctsoccdf[ccd_t_modctsoccdf.occupancy >= 0]["modification"].value_counts())
print(ccd_n_modctsoccdf[ccd_n_modctsoccdf.occupancy >= 0]["modification"].value_counts())

#%% Perform melting point analysis
ProteinStabilityAnalysis.melting_point_analysis(ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, nonccdprotein)