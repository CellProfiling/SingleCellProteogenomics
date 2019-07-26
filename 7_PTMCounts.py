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

# Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
count_or_rpkm = "Tpms" # so that the results match for cross-gene comparisons
adata, phases_filt = read_counts_and_phases(dd, count_or_rpkm, False, "")
qc_filtering(adata, False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

allccd_transcript_regulated = set(ccd_gene_names(allccd_transcript_regulated))
dianaccd_transcript_regulated = set(ccd_gene_names(dianaccd_transcript_regulated))
dianaccd_nontranscript_regulated = set(ccd_gene_names(dianaccd_nontranscript_regulated))
genes_analyzed = set(ccd_gene_names(genes_analyzed))
ccd_regev_filtered = set(ccd_gene_names(ccd_regev_filtered))
ccd_filtered = set(ccd_gene_names(ccd_filtered))
nonccd_filtered = set(ccd_gene_names(nonccd_filtered))

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in the protein group results
search = "2019-06-19-15-54-09_CommonPtmsWithOccupancy\\Task1-SearchTask"
search = "2019-07-12_NoMetalOccupancyFix\\Task2-SearchTask"
filenameCommon = 'C:\\Users\\antho\\Box\ProjectData\\Variability\\U2OS_gptmd_search\\' + search + '\\AllProteinGroups.tsv'
filenameCommonAndLessCommon = 'C:\\Users\\antho\\Box\ProjectData\\Variability\\U2OS_gptmd_search\\U2OSGptmdSearchAllProteinGroups.tsv'
file = pd.read_csv(filenameCommon, sep="\t", index_col=False)
# targets=file[(file["Protein Decoy/Contaminant/Target"] == "T") & (file["Protein QValue"] <= 0.01)]
modifiedProteins = targets[targets["Sequence Coverage with Mods"].str.replace("[","") != targets["Sequence Coverage with Mods"]]
modifications = [re.findall('\[.*?\]',s) for s in modifiedProteins["Sequence Coverage with Mods"]]
unique_mods = set([item for sublist in modifications for item in sublist])

genes = list(targets["Gene"])
seqmods = list(targets["Sequence Coverage with Mods"])
seqs = list(targets["Sequence Coverage"])
modinfos = list(targets["Modification Info List"])
genemods = {}
blacklist = ["oxidation", "deamidation", "ammonia loss", "water loss", "carbamyl", "carbamidomethyl", # artifacts
    "fe[i", "zinc", "cu[i" , # metals
    "sodium", "potassium", "calcium", "magnesium",] # counterions
    
for idx in range(len(genes)):
    genesplit = str(genes[idx]).split("|")
    seqmod = seqmods[idx].split("|")
    seq = seqs[idx].split("|")
    modinfo = [x.strip(';') for x in str(modinfos[idx]).split("|")]
    for idxx, genegene in enumerate(genesplit):
        modstemp, occupancytemp = [], []
        modpepts, totpepts = [], []

        if idxx < len(modinfo):
            modis = [m for m in modinfo[idxx].split(";") if m != "nan"]
            aanums = [m[len("#aa"):m.index('[')] for m in modis]
            modstemp = [m[(m.index('[') + 1):m.index(",info")] for m in modis]
            occupancytemp = [m[(m.index("occupancy=") + len("occupancy=")):(m.index("occupancy=") + len("occupancy=") + len("#.##"))] for m in modis]
            modpepts = [m[(m.index("(", m.index("occupancy=")) + 1):(m.index("/", m.index("occupancy=")))] for m in modis]
            totpepts = [m[(m.index("/", m.index("occupancy=")) + 1):(m.index(")", m.index("occupancy=")))] for m in modis]
        else: continue

        # mods = [m for m in re.findall('\[.*?\]', seqmod[idxx]) if not any(mod in m.lower() for mod in blacklist)]
        covered_fraction = float(sum([1 for c in seq[idxx] if c.islower()])) / float(len(seq[idxx]))
        effective_length = float(len(seq[idxx])) * covered_fraction
        # mods_per_eff_base = math.log10(float(len(mods) + 1) / float(effective_length))
        mods_per_eff_base = float(len(mods)) / float(effective_length)
        mods, occupancies, modpeptsss, totpeptsss = [],[],[],[]
        for oidx, occcc in enumerate(occupancytemp):
            isblacklisted = any(mod in modstemp[oidx].lower() for mod in blacklist)
            isonehitwonder = int(modpepts[oidx]) <= 1
            if not "(" in occcc and not isblacklisted and not isonehitwonder: 
                mods.append(modstemp[oidx])
                occupancies.append(float(occcc))
                modpeptsss.append(int(modpepts[oidx]))
                totpeptsss.append(int(totpepts[oidx]))
        if genegene not in genemods: 
            genemods[genegene] = ([mods_per_eff_base], occupancies, mods, 
                modpeptsss, totpeptsss)
        else:
            genemods[genegene][0].append(mods_per_eff_base)
            genemods[genegene][1].extend(occupancies)
            genemods[genegene][2].extend(mods)
            genemods[genegene][3].extend(modpeptsss)
            genemods[genegene][4].extend(totpeptsss)

print(f"{str(len(targets))} proteins")
print(f"{str(len(modifiedProteins))} modified proteins ({str(round(float(len(modifiedProteins))/float(len(targets))*100,2))}%)")
print(f"{str(len([g for g in genemods.keys() if g in allccd_transcript_regulated]))}: number of all transcript regulated CCD genes of {len(allccd_transcript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in dianaccd_transcript_regulated]))}: number of transcript regulated CCD genes from Diana's study of {len(dianaccd_transcript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in dianaccd_nontranscript_regulated]))}: number of non-transcript regulated CCD genes from Diana's study of {len(dianaccd_nontranscript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in genes_analyzed]))}: number of proteins of {len(genes_analyzed)} detected.")

unambigenes, modcts, modsss = [], [], []
protmodgene, modmod, occocc = [], [], []
modpeps, totpeps = [], []
for gene in genemods.keys():
    unambigenes.append(gene)
    modcts.append(np.median(genemods[gene][0]))
    modsss.append(", ".join(genemods[gene][2]))
    for modidx, mod in enumerate(genemods[gene][2]):
        protmodgene.append(gene)
        modmod.append(mod)
        occocc.append(genemods[gene][1][modidx])
        modpeps.append(genemods[gene][3][modidx])
        totpeps.append(genemods[gene][4][modidx])
        
df = pd.DataFrame({
    "gene" : unambigenes, 
    "modcts" : modcts, 
    "modlist" : modsss})
occdf = pd.DataFrame({
    "gene" : protmodgene,
    "modification" : modmod,
    "occupancy" : occocc,
    "modpeps" : modpeps,
    "totpeps" : totpeps})

# all genes
all_modctsdf = df[df["gene"].isin(genes_analyzed)]
all_modctsdf.to_csv("output/ProteinModCounts.csv")
all_modcts = all_modctsdf["modcts"]
all_modcts.hist()
plt.title("All, modcts")
plt.show()
plt.close()
all_modctsoccdf = occdf[occdf["gene"].isin(genes_analyzed)]
all_modctsoccdf.to_csv("output/ProteinModOccupancies.csv")
all_modocc = all_modctsoccdf["occupancy"]
all_modocc.hist()
plt.title("All occupancies")
plt.show()
plt.close()
all_modctsoccdf[all_modctsoccdf.occupancy == 1]["modpeps"].hist(bins=60)
plt.title("All modified peptide counts (occupancy = 1)")
plt.xlabel("Modified peptide count")
plt.ylabel("Modified residues within bin")
plt.show(); plt.close()
all_modctsoccdf["totpeps"].hist()
plt.title("All total peptide counts")
plt.show(); plt.close()

# regev
ccd_rt_modctsdf = df[df["gene"].isin(ccd_regev_filtered)]
ccd_rt_modctsdf.to_csv("output/ProteinModCountsTransRegRegevCcd.csv")
ccd_rt_modcts = ccd_rt_modctsdf["modcts"]
ccd_rt_modcts.hist()
plt.title("Regev CCD Transcript Reg, modcts")
plt.show()
plt.close()
ccd_rt_modctsoccdf = occdf[occdf["gene"].isin(ccd_regev_filtered)]
ccd_rt_modctsoccdf.to_csv("output/ProteinModRegevCcdOccupancies.csv")
ccd_rt_modocc = ccd_rt_modctsoccdf["occupancy"]
ccd_rt_modocc.hist()
plt.title("Regev CCD Transcript Reg, occupancies")
plt.show()
plt.close()

ccd_at_modctsdf = df[df["gene"].isin(allccd_transcript_regulated)]
ccd_at_modctsdf.to_csv("output/ProteinModCountsTransRegAllCcd.csv")
ccd_at_modcts = ccd_at_modctsdf["modcts"]
ccd_at_modcts.hist()
plt.title("All CCD Transcript Reg, modcts")
plt.show()
plt.close()
ccd_at_modctsoccdf = occdf[occdf["gene"].isin(allccd_transcript_regulated)]
ccd_at_modctsoccdf.to_csv("output/ProteinModAllCcdOccupancies.csv")
ccd_at_modocc = ccd_at_modctsoccdf["occupancy"]
ccd_at_modocc.hist()
plt.title("All CCD Transcript Reg, occupancies")
plt.show()
plt.close()

ccd_t_modctsdf = df[df["gene"].isin(dianaccd_transcript_regulated)]
ccd_t_modctsdf.to_csv("output/ProteinModCountsTransRegDianaCcd.csv")
ccd_t_modcts = ccd_t_modctsdf["modcts"]
ccd_t_modcts.hist()
plt.title("Diana's CCD Transcript Reg, modcts")
plt.show()
plt.close()
ccd_t_modctsoccdf = occdf[occdf["gene"].isin(dianaccd_transcript_regulated)]
ccd_t_modctsoccdf.to_csv("output/ProteinModDianaCcdOccupancies.csv")
ccd_t_modocc = ccd_t_modctsoccdf["occupancy"]
ccd_t_modocc.hist()
plt.title("Diana's CCD Transcript Reg, occupancies")
plt.show()
plt.close()

ccd_n_modctsdf = df[df["gene"].isin(dianaccd_nontranscript_regulated)]
ccd_n_modctsdf.to_csv("output/ProteinModCountsNonTransRegDianaCcd.csv")
ccd_n_modcts = ccd_n_modctsdf["modcts"]
ccd_n_modcts.hist()
plt.title("Diana's CCD Non-Transcript Reg, modcts")
plt.show()
plt.close()
ccd_n_modctsoccdf = occdf[occdf["gene"].isin(dianaccd_nontranscript_regulated)]
ccd_n_modctsoccdf.to_csv("output/ProteinModDianaNonCcdOccupancies.csv")
ccd_n_modocc = ccd_n_modctsoccdf["occupancy"]
ccd_n_modocc.hist()
plt.title("Diana's CCD Non-Transcript Reg, occupancies")
plt.show()
plt.close()

# Just number of proteins with mods
print(f"{str(len(all_modcts[all_modcts == 0]))} modified proteins ({str(round(float(len(all_modcts[all_modcts == 0]))/float(len(all_modcts))*100,2))}%) in all genes detected")
print(f"{str(len(ccd_rt_modcts[ccd_rt_modcts == 0]))} modified proteins ({str(round(float(len(ccd_rt_modcts[ccd_rt_modcts == 0]))/float(len(ccd_rt_modcts))*100,2))}%) in Regev CCD")
print(f"{str(len(ccd_at_modcts[ccd_at_modcts == 0]))} modified proteins ({str(round(float(len(ccd_at_modcts[ccd_at_modcts == 0]))/float(len(ccd_at_modcts))*100,2))}%) in all transcript CCD genes")
print(f"{str(len(ccd_t_modcts[ccd_t_modcts == 0]))} modified proteins ({str(round(float(len(ccd_t_modcts[ccd_t_modcts == 0]))/float(len(ccd_t_modcts))*100,2))}%) in Diana's CCD proteins")
print(f"{str(len(ccd_n_modcts[ccd_n_modcts == 0]))} modified proteins ({str(round(float(len(ccd_n_modcts[ccd_n_modcts == 0]))/float(len(ccd_n_modcts))*100,2))}%) in Diana's non-CCD proteins")

# The t-test isn't good here because these counts don't form a normal distribution... it works with the log(x+1 / len) distribution
print("COUNT OF MODIFICATIONS PER COVERED RESIDUE PER PROTEIN")
print(f"mean number of mods / seq length for all proteins: {np.mean(all_modcts)}")
print(f"mean number of mods / seq length for all transcriptionally regulated CCD: {np.mean(ccd_at_modcts)}")
print(f"mean number of mods / seq length for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modcts)}")
print(f"mean number of mods / seq length for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modcts)}")
t, p = scipy.stats.ttest_ind(ccd_t_modcts, ccd_n_modcts)
print(f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}")
t, p = scipy.stats.ttest_ind(all_modcts, ccd_t_modcts)
print(f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}")
print()
print(f"median number of mods / seq length for all proteins: {np.median(all_modcts)}")
print(f"median number of mods / seq length for all transcriptionally regulated CCD: {np.median(ccd_at_modcts)}")
print(f"median number of mods / seq length for Diana's transcriptionally regulated CCD: {np.median(ccd_t_modcts)}")
print(f"median number of mods / seq length for Diana's non-transcriptionally regulated: {np.median(ccd_n_modcts)}")
t, p = scipy.stats.kruskal(ccd_at_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all transcriptionally reg. and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(ccd_t_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of Diana's transcriptionally reg. and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_t_modcts)
print(f"one-sided kruskal for same medians of all genes and Diana's transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_at_modcts)
print(f"one-sided kruskal for same medians of all genes and all transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_rt_modcts)
print(f"one-sided kruskal for same medians of all genes and regev transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_t_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_at_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}")
print()

# The t-test isn't good here because these counts don't form a normal distribution... it works with the log(x+1 / len) distribution
print("OCCUPANCY PER MODIFIED BASE")
print("Note: I excluded artifact modifications and required there to be")
print("at least 2 peptides with the modification to count it.")
print()
print(f"mean occupancy per modified residue for all proteins: {np.mean(all_modocc)}")
print(f"mean occupancy per modified residue for all transcriptionally regulated CCD: {np.mean(ccd_at_modocc)}")
print(f"mean occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modocc)}")
print(f"mean occupancy per modified residue for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modocc)}")
t, p = scipy.stats.ttest_ind(ccd_t_modocc, ccd_n_modocc)
print(f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}")
t, p = scipy.stats.ttest_ind(all_modocc, ccd_t_modocc)
print(f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}")
print()
print(f"median occupancy per modified residue for all proteins: {np.median(all_modocc)}")
print(f"median occupancy per modified residue for all transcriptionally regulated CCD: {np.median(ccd_at_modocc)}")
print(f"median occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.median(ccd_t_modocc)}")
print(f"median occupancy per modified residue for Diana's non-transcriptionally regulated: {np.median(ccd_n_modocc)}")
t, p = scipy.stats.kruskal(ccd_at_modocc, ccd_n_modocc)
print(f"one-sided kruskal for same medians of all transcriptionally reg. and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(ccd_t_modocc, ccd_n_modocc)
print(f"one-sided kruskal for same medians of Diana's transcriptionally reg. and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_t_modocc)
print(f"one-sided kruskal for same medians of all genes and Diana's transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_at_modocc)
print(f"one-sided kruskal for same medians of all genes and all transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_rt_modocc)
print(f"one-sided kruskal for same medians of all genes and regev transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_n_modocc)
print(f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_t_modocc, ccd_n_modocc)
print(f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(all_modocc, ccd_at_modocc, ccd_n_modocc)
print(f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}")
print()

# Generate a histogram of effective mod counts with bins normalized to 1'''
def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

# bar histogram (modcts)
bins=np.histogram(np.hstack((all_modcts, ccd_t_modcts, ccd_n_modcts)), bins=40)[1] #get the bin edges
plt.hist(all_modcts, bins=bins, weights=weights(all_modcts), histtype="bar", alpha=0.4, label="All Proteins")
plt.hist(ccd_t_modcts, bins=bins, weights=weights(ccd_t_modcts), histtype="bar", alpha=0.4, label="Transcript Reg. CCD")
plt.hist(ccd_n_modcts, bins=bins, weights=weights(ccd_n_modcts), histtype="bar", alpha=0.5, label="Non-Transcript Reg. CCD")
plt.legend(loc="upper left")
plt.xlabel("Modifications per Covered Base")
plt.savefig("figures/ModsPerCoveredBaseHistogram.png")
plt.show()
plt.close()

# bar histogram (modoccupancy)
bins=np.histogram(np.hstack((all_modocc, ccd_t_modocc, ccd_n_modocc)), bins=20)[1] #get the bin edges
plt.hist(all_modocc, bins=bins, weights=weights(all_modocc), histtype="bar", alpha=0.4, label="All Proteins")
plt.hist(ccd_t_modocc, bins=bins, weights=weights(ccd_t_modocc), histtype="bar", alpha=0.4, label="Transcript Reg. CCD")
plt.hist(ccd_n_modocc, bins=bins, weights=weights(ccd_n_modocc), histtype="bar", alpha=0.5, label="Non-Transcript Reg. CCD")
plt.legend(loc="upper left")
plt.xlabel("Occupancy per Modified Residue")
plt.savefig("figures/OccupancyPerModHistogram.png")
plt.show()
plt.close()

# line histogram (mod counts)
xx = [all_modcts, ccd_t_modcts, ccd_n_modcts]
ww = [weights(all_modcts), weights(ccd_t_modcts), weights(ccd_n_modcts)]
ll = ["All Proteins", "Trans Reg CCD", "Non-Trans Reg CCD"]
binEdges=np.histogram(np.hstack((all_modcts, ccd_t_modcts, ccd_n_modcts)), bins=40)[1] #get the bin edges
allHist = np.histogram(all_modcts, bins=binEdges, weights=weights(all_modcts))[0]
ccdtHist = np.histogram(ccd_t_modcts, bins=binEdges, weights=weights(ccd_t_modcts))[0]
ccdnHist = np.histogram(ccd_n_modcts, bins=binEdges, weights=weights(ccd_n_modcts))[0]
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.figure(figsize=(10,10))
plt.plot(bincenters, allHist, label="All Proteins")
plt.plot(bincenters, ccdtHist, label="Trans Reg CCD")
plt.plot(bincenters, ccdnHist, label="Non-Trans Reg CCD")
# plt.hist(ccd_t_modcts, bins=bins, weights=weights(ccd_t_modcts), histtype="bar", alpha=0.4, label="Transcript Reg. CCD")
# plt.hist(ccd_n_modcts, bins=bins, weights=weights(ccd_n_modcts), histtype="bar", alpha=0.5, label="Non-Transcript Reg. CCD")
plt.legend(loc="upper right")
plt.xlabel("Modifications per Covered Base")
plt.savefig("figures/ModsPerCoveredBaseLineHistogram.png")
plt.show()
plt.close()

# line histogram (mod occupancy)
xx = [all_modocc, ccd_t_modocc, ccd_n_modocc]
ww = [weights(all_modocc), weights(ccd_t_modocc), weights(ccd_n_modocc)]
ll = ["All Proteins", "Trans Reg CCD", "Non-Trans Reg CCD"]
binEdges=np.histogram(np.hstack((all_modocc, ccd_t_modocc, ccd_n_modocc)), bins=40)[1] #get the bin edges
allHist = np.histogram(all_modocc, bins=binEdges, weights=weights(all_modocc))[0]
ccdtHist = np.histogram(ccd_t_modocc, bins=binEdges, weights=weights(ccd_t_modocc))[0]
ccdnHist = np.histogram(ccd_n_modocc, bins=binEdges, weights=weights(ccd_n_modocc))[0]
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.figure(figsize=(10,10))
plt.plot(bincenters, allHist, label="All Proteins")
plt.plot(bincenters, ccdtHist, label="Trans Reg CCD")
plt.plot(bincenters, ccdnHist, label="Non-Trans Reg CCD")
# plt.hist(ccd_t_modocc, bins=bins, weights=weights(ccd_t_modocc), histtype="bar", alpha=0.4, label="Transcript Reg. CCD")
# plt.hist(ccd_n_modocc, bins=bins, weights=weights(ccd_n_modocc), histtype="bar", alpha=0.5, label="Non-Transcript Reg. CCD")
plt.legend(loc="upper right")
plt.xlabel("Modifications per Covered Base")
plt.savefig("figures/ModsPerCoveredBaseLineHistogram.png")
plt.show()
plt.close()

# Generate boxplots for effective mod counts
mmmm = np.concatenate((all_modcts, ccd_t_modcts, ccd_rt_modcts, ccd_at_modcts, ccd_n_modcts))
cccc = (["All"] * len(all_modcts))
cccc.extend(["Diana's Transcript CCD"] * len(ccd_t_modcts))
cccc.extend(["Regev Transcript CCD"] * len(ccd_rt_modcts))
cccc.extend(["All Transcript CCD"] * len(ccd_at_modcts))
cccc.extend(["Non-Transc CCD"] * len(ccd_n_modcts))
moddf = pd.DataFrame({"modcts": mmmm, "category" : cccc})
boxplot = moddf.boxplot("modcts", by="category", figsize=(12, 8), showfliers=False)
boxplot.set_xlabel("Protein Set", size=36,fontname='Arial')
boxplot.set_ylabel("Modifications per Covered Base", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/ModsPerCoveredBaseBoxplot.png")
plt.show()
plt.close()

# Generate boxplots for mod occupancies
mmmm = np.concatenate((all_modocc, ccd_t_modocc, ccd_rt_modocc, ccd_at_modocc, ccd_n_modocc))
cccc = (["All"] * len(all_modocc))
cccc.extend(["Diana's Transcript\nCCD"] * len(ccd_t_modocc))
cccc.extend(["Regev Transcript\nCCD"] * len(ccd_rt_modocc))
cccc.extend(["All Transcript\nCCD"] * len(ccd_at_modocc))
cccc.extend(["Non-Transc\nCCD"] * len(ccd_n_modocc))
moddf = pd.DataFrame({"modocc": mmmm, "category" : cccc})
boxplot = moddf.boxplot("modocc", by="category", figsize=(12, 8), showfliers=False)
boxplot.set_xlabel("Protein Set", size=36,fontname='Arial')
boxplot.set_ylabel("Occupancy per Modified Residue", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/ModsPerCoveredBaseBoxplot.png")
plt.show()
plt.close()

# Generate violin plots
# dddd = [sorted(x) for x in [all_modcts, ccd_t_modcts, ccd_n_modcts]]
# # cccc = (["All"] * len(all_modcts))
# # cccc.extend(["Transcript CCD"] * len(ccd_t_modcts))
# # cccc.extend(["Non-Transc CCD"] * len(ccd_n_modcts))
# ax = plt.violinplot(dddd, showextrema=False, showmedians=True)
# plt.ylim(top=0.06)
# # moddf = pd.DataFrame({"modcts": mmmm, "category" : cccc})
# # boxplot = moddf.boxplot("modcts", by="category", figsize=(12, 8), showfliers=False)
# ax.xlabel("Protein Set", size=36,fontname='Arial')
# ax.ylabel("Modifications per Covered Base", size=36,fontname='Arial')
# # boxplot.tick_params(axis="both", which="major", labelsize=16)
# # plt.title("")
# plt.savefig("figures/ModsPerCoveredBaseViolin.png")
# plt.show()
# plt.close()

# print(f"{str(len(modifiedPeptides))} interesting modified peptides ({str(round(float(len(modPeptides))/float(len(file))*100,2))}%)")

#%% [markdown]
# Bulk U2OS analysis results
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


#%%
