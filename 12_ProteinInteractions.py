#%% Imports
from imports import *
from Bio import SeqIO
import sys
import re
import math
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names

#%% Import the genes names we're analyzing
allccd_transcript_regulated = np.array(pd.read_csv("output/allccd_transcript_regulated.csv")["gene"])
dianaccd_transcript_regulated = np.array(pd.read_csv("output/ccd_transcript_regulated.csv")["gene"])
dianaccd_nontranscript_regulated = np.array(pd.read_csv("output/ccd_nontranscript_regulated.csv")["gene"])
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])
gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
ccd_regev=list(gene_info[gene_info["name"].isin(pd.read_csv("input/ccd_regev.txt")["gene"])]["gene_id"])
ccd=list(gene_info[gene_info["name"].isin(pd.read_csv("input/ccd_genes.txt")["gene"])]["gene_id"])
nonccd=list(gene_info[gene_info["name"].isin(pd.read_csv("input/nonccd_genes.txt")["gene"])]["gene_id"])

# Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
count_or_rpkm = "Tpms" # so that the results match for cross-gene comparisons
adata, phases_filt = read_counts_and_phases(dd, count_or_rpkm, False, "")
qc_filtering(adata, False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

allccd_transcript_regulated_names = set(ccd_gene_names(allccd_transcript_regulated))
dianaccd_transcript_regulated_names = set(ccd_gene_names(dianaccd_transcript_regulated))
dianaccd_nontranscript_regulated_names = set(ccd_gene_names(dianaccd_nontranscript_regulated))
diananonccd_names = set(ccd_gene_names(nonccd_filtered))
genes_analyzed_names = set(ccd_gene_names(genes_analyzed))
ccd_regev_filtered_names = set(ccd_gene_names(ccd_regev_filtered))
ccd_filtered_names = set(ccd_gene_names(ccd_filtered))

#%% Import the protein-protein interactions
interactionHuriUnion = 'C:\\Users\\antho\\Box\ProjectData\\ProteinInteractions\\HuRI\\HI-Union.tsv'
interactions = {}
with open(interactionHuriUnion) as file:
    for line in file:
        if not line.startswith("ENSG"): continue
        linesplit = line.split('\t')
        ensg1 = linesplit[0]
        ensg2 = linesplit[1]
        if ensg1 in interactions: interactions[ensg1].add(ensg2)
        else: interactions[ensg1] = set([ensg2])
interactionCounts = dict([(ensg, len(interactions[ensg])) for ensg in interactions])

#%% 
# Idea: compare the distributions of counts of interactions for the different sets of proteins
# Exec: create numpy arrays of the counts for the genes in the set
# Output: histograms
all_interaction_counts = list(interactionCounts.values())
ccd_interaction_counts = [interactionCounts[ensg] for ensg in ccd if ensg in interactionCounts]
ccd_regev_interaction_counts = [interactionCounts[ensg] for ensg in ccd_regev if ensg in interactionCounts]
nonccd_interaction_counts = [interactionCounts[ensg] for ensg in nonccd if ensg in interactionCounts]
ccd_a_interaction_counts = [interactionCounts[ensg] for ensg in allccd_transcript_regulated if ensg in interactionCounts]
ccd_t_interaction_counts = [interactionCounts[ensg] for ensg in dianaccd_transcript_regulated if ensg in interactionCounts]
ccd_n_interaction_counts = [interactionCounts[ensg] for ensg in dianaccd_nontranscript_regulated if ensg in interactionCounts]


# Generate a histogram of effective mod counts with bins normalized to 1'''
def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

# bar histogram (modcts)
bins=np.histogram(np.hstack((all_interaction_counts, ccd_interaction_counts, ccd_regev_interaction_counts, nonccd_interaction_counts)), bins=40)[1] #get the bin edges
plt.hist(all_interaction_counts, bins=bins, weights=weights(all_interaction_counts), histtype="bar", alpha=0.4, label="All Proteins")
plt.hist(ccd_interaction_counts, bins=bins, weights=weights(ccd_interaction_counts), histtype="bar", alpha=0.4, label="Protein Reg. CCD")
plt.hist(ccd_regev_interaction_counts, bins=bins, weights=weights(ccd_regev_interaction_counts), histtype="bar", alpha=0.5, label="Regev Transcript Reg. CCD")
plt.hist(nonccd_interaction_counts, bins=bins, weights=weights(nonccd_interaction_counts), histtype="bar", alpha=0.5, label="Protein Reg. Non-CCD")
plt.legend(loc="upper left")
plt.xlabel("Interactions per Protein")
plt.savefig("figures/ProteinInteractionsHistogram1.png")
plt.show()
plt.close()

bins=np.histogram(np.hstack((all_interaction_counts, ccd_a_interaction_counts, ccd_t_interaction_counts, ccd_n_interaction_counts)), bins=40)[1] #get the bin edges
plt.hist(all_interaction_counts, bins=bins, weights=weights(all_interaction_counts), histtype="bar", alpha=0.4, label="All Proteins")
plt.hist(ccd_a_interaction_counts, bins=bins, weights=weights(ccd_a_interaction_counts), histtype="bar", alpha=0.4, label="All Transcript Reg. CCD")
plt.hist(ccd_t_interaction_counts, bins=bins, weights=weights(ccd_t_interaction_counts), histtype="bar", alpha=0.5, label="Diana's Transcript Reg. CCD")
plt.hist(ccd_n_interaction_counts, bins=bins, weights=weights(ccd_n_interaction_counts), histtype="bar", alpha=0.5, label="Diana's Non-Transcript Reg. CCD")
plt.legend(loc="upper left")
plt.xlabel("Interactions per Protein")
plt.savefig("figures/ProteinInteractionsHistogram2.png")
plt.show()
plt.close()

mmmm = np.concatenate((all_interaction_counts, ccd_interaction_counts, ccd_regev_interaction_counts, nonccd_interaction_counts, ccd_a_interaction_counts, ccd_t_interaction_counts, ccd_n_interaction_counts))
cccc = (["All\nProteins"] * len(all_interaction_counts))
cccc.extend(["Diana's\nCCD"] * len(ccd_interaction_counts))
cccc.extend(["Regev\nTranscript\nReg. CCD"] * len(ccd_regev_interaction_counts))
cccc.extend(["Diana's\nNon-CCD"] * len(nonccd_interaction_counts))
cccc.extend(["All\nTranscript\nCCD"] * len(ccd_a_interaction_counts))
cccc.extend(["Diana's\nTranscript\nCCD"] * len(ccd_t_interaction_counts))
cccc.extend(["Diana's\nNon-Transcript\nCCD"] * len(ccd_n_interaction_counts))
moddf = pd.DataFrame({"time": mmmm, "category" : cccc})
boxplot = moddf.boxplot("time", by="category", figsize=(12, 8), showfliers=True)
boxplot.set_xlabel("", size=36,fontname='Arial')
boxplot.set_ylabel("Interactions per Protein", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/ProteinInteractionsBoxplot.png")
plt.show()
plt.close()

#%%
# Idea: compare the distrubutions of multilocalizing and single localizing proteins
# Exec: pandas
# Output: histogram and boxplot
compartment_counts = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\ProteinLocalizations\\CompartmentCountsThul2017.csv")
multilocalizing = compartment_counts[compartment_counts["compartment_count"] > 1]["ensg"]
singlelocalizing = compartment_counts[compartment_counts["compartment_count"] == 1]["ensg"]
multiloc_counts = [interactionCounts[ensg] for ensg in multilocalizing if ensg in interactionCounts]
singleloc_counts = [interactionCounts[ensg] for ensg in singlelocalizing if ensg in interactionCounts]

# bar histogram (modcts)
bins=np.histogram(np.hstack((multiloc_counts, singleloc_counts)), bins=40)[1] #get the bin edges
plt.hist(multiloc_counts, bins=bins, weights=weights(multiloc_counts), histtype="bar", alpha=0.5, label="Multilocalizing Proteins")
plt.hist(singleloc_counts, bins=bins, weights=weights(singleloc_counts), histtype="bar", alpha=0.5, label="Singlelocalizing Proteins")
plt.legend(loc="upper left")
plt.xlabel("Interactions per Protein")
plt.savefig("figures/ProteinInteractionsHistogramMultiSingleLocalizing.png")
plt.show()
plt.close()

def interaction_multiloc_boxplot(outliers):
    mmmm = np.concatenate((multiloc_counts, singleloc_counts))
    cccc = (["Multilocalizing\nProteins"] * len(multiloc_counts))
    cccc.extend(["Singlelocalizing\nProteins"] * len(singleloc_counts))
    moddf = pd.DataFrame({"time": mmmm, "category" : cccc})
    boxplot = moddf.boxplot("time", by="category", figsize=(12, 8), showfliers=outliers)
    boxplot.set_xlabel("", size=36,fontname='Arial')
    boxplot.set_ylabel("Interactions per Protein", size=36,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=16)
    plt.title("")
    plt.savefig(f"figures/ProteinInteractionsBoxplotOutliers{outliers}.png")
    plt.show()
    plt.close()

interaction_multiloc_boxplot(True)
interaction_multiloc_boxplot(False)

#%%
