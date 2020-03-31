#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, Loaders, RNADataPreparation, PTMAnalysis
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

#%% Import the genes names we're analyzing
# Read in RNA-Seq data again and the CCD gene lists
plate, valuetype, use_spikeins, biotype_to_use = "All", "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spikeins, biotype_to_use)
adata, phasesfilt = RNADataPreparation.qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)
import_dict = Loaders.load_ptm_and_stability(adata)

wp_max_pol = import_dict["wp_max_pol"]
ccdtranscript, nonccdtranscript = import_dict["ccdtranscript_names"], import_dict["nonccdtranscript_names"]
ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated = import_dict["ccdprotein_transcript_regulated"], import_dict["ccdprotein_nontranscript_regulated"]
genes_analyzed = import_dict["genes_analyzed_names"]
ccd_regev_filtered, ccd_filtered = import_dict["ccd_regev_filtered_names"], import_dict["ccd_filtered_names"] 
ccdprotein, nonccdprotein = import_dict["ccdprotein_names"], import_dict["nonccdprotein_names"]

utils.save_gene_names_by_category(adata)

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in MetaMorpheus results
genemodsBulk = PTMAnalysis.analyze_ptms("input/raw/U2OSBulkAllProteinGroups.tsv")
genemodsPhospho = PTMAnalysis.analyze_ptms("input/raw/U2OSPhosphoAllProteinGroups.tsv")

# Analyze the modification occupancy for each PTM site
dfBulk, occdfBulk = PTMAnalysis.process_genemods(genemodsBulk)
dfPhospho, occdfPhospho = PTMAnalysis.process_genemods(genemodsPhospho)

# Compare the PTM regulation between transcript and non-transcript regulated groups
PTMAnalysis.compare_ptm_regulation(dfBulk, occdfBulk, "Bulk")
PTMAnalysis.compare_ptm_regulation(dfPhospho, occdfPhospho, "Phospho")

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


#%% [markdown]
# # Bulk U2OS analysis results (modification counts)
#
# It turns out the significance went away when I filtered out an artifact I missed (iron adducts),
# and then after checking all transcript-regulated CCD genes, rather than just the 
# set in Diana's study, it all went away.
# 
# On the bright side, there are some really interesting modifications even in the 
# bulk sequencing: plenty of phosphorylations and even a nitrosylation on BUB1B that's
# not in UniProt.
#
# This also doesn't work with the metals anymore, nor does it work with
# just phosphorylations.
#
# When I drop the non-CCD proteins in there, there is




#%% [markdown]
# Some thoughts on occupancy. It looks like there's a significant difference in occupancies between
# the different protein sets. There are a lot of ones (fully occupied), and I'm sure that's biased by
# there less labile modifications like n-acetylation. Even so, when I look at the modifications
# at high occupancy, I'm surprised to find that the non-CCD set is pretty rich with a few phosphos,
# a hydroxylation, and a citrullination.
#
# Well, I could do the count of mods with occupancy > 0.5. It looks like the ratio of those mod sites
# is pretty high.
# 
# But I wonder if I should take this all with a grain of salt. The number of mod sites I'm detecting
# is pretty low for the sub categories. Yet, there are still some interesting phosphorylations in
# there, and the N-acetylations don't dominate as much as I expected.

#%%
# Idea: Do a modification class based occupancy comparison
# Exec: Group the modifications by type; calculate the differences in occupancy between categories
# Output: Table of modification categories and statistical results comparing the categories

