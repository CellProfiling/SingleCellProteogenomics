#%% Imports
from imports import *
from Bio import SeqIO
import requests
import sys
import re
import math

#%% Import the genes names we're analyzing
allccd_transcript_regulated = np.array(pd.read_csv("output/allccd_transcript_regulated.csv")["gene"])
dianaccd_transcript_regulated = np.array(pd.read_csv("output/ccd_transcript_regulated.csv")["gene"])
dianaccd_nontranscript_regulated = np.array(pd.read_csv("output/ccd_nontranscript_regulated.csv")["gene"])
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)


#%% Download uniprot PTM information
BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/uniprot/'
TOOL_ENDPOINT = '/uploadlists/'

def query_uniprot(gene_list, filename):
    query = f'proteome:UP000005640 AND (gene:{" OR gene:".join(gene_list)})'

    payload = {
        'query': query,
        'format': 'tab',
        # 'columns': 'id,entry_name,reviewed,protein_names,organism,comment(PTM),feature(GLYCOSYLATION),feature(SIGNAL)',
        'columns': 'id,entry_name,reviewed,genes(PREFERRED),comment(PTM),feature(GLYCOSYLATION),feature(SIGNAL)',
        }

    result = requests.get(BASE + KB_ENDPOINT, params=payload, stream=True)
    result.raise_for_status() # throw an error for bad status code
    with open(filename, "wb") as f:
        for block in result.iter_content(1024):
            f.write(block)

# query_uniprot(ccd_transcript_regulated, "output/ccd_transcript_regulated_uniprot.tsv")
# query_uniprot(ccd_nontranscript_regulated, "output/ccd_nontranscript_regulated_uniprot.tsv")
# query_uniprot(genes_analyzed, "output/genes_analyzed.tsv")

#%% Import UniProt PTM information for the proteins
# ccd_t_uniprot = pd.read_csv("output/ccd_transcript_regulated_uniprot.tsv", sep="\t")
# ccd_n_uniprot = pd.read_csv("output/ccd_nontranscript_regulated_uniprot.tsv", sep="\t")
# ccd_t_uniprot_ptreg = (pd.notnull(ccd_t_uniprot["Glycosylation"])) | (pd.notnull(ccd_t_uniprot["Post-translational modification"])) | (pd.notnull(ccd_t_uniprot["Signal peptide"]))
# print(len(np.array(ccd_t_uniprot)[np.array(ccd_t_uniprot_ptreg)]) / float(len(ccd_t_uniprot)))
# ccd_n_uniprot_ptreg = (pd.notnull(ccd_n_uniprot["Glycosylation"])) | (pd.notnull(ccd_n_uniprot["Post-translational modification"])) | (pd.notnull(ccd_n_uniprot["Signal peptide"]))
# print(len(np.array(ccd_n_uniprot)[np.array(ccd_n_uniprot_ptreg)]) / float(len(ccd_n_uniprot)))


#%% Yeah, that didn't work, but here's some C# code to get at the number of PTMs using mzLib

# string wd = @"C:\Users\antho\Documents\Projects\CellCycle\FucciRNA";
# Loaders.LoadElements();
# var psi = Loaders.LoadPsiMod(@"C:\Users\antho\Documents\Programs\SpritzSnake\PSI-MOD.obo.xml");
# var ptms = Loaders.LoadUniprot(@"C:\Users\antho\Documents\Programs\SpritzSnake\ptmlist.txt", Loaders.GetFormalChargesDictionary(psi));
# List<Protein> proteins = ProteinDbLoader.LoadProteinXML(@"C:\Users\antho\Documents\Programs\SpritzSnake\human.protein.xml.gz", true, DecoyType.None, ptms, false, null, out var un);
# Dictionary<string, int> ptm_counts = new Dictionary<string, int>();
# foreach (Protein p in proteins)
# {
#     string kk = p.GeneNames.Count() == 0 ? p.Accession : p.GeneNames.First().Item2;
#     if (ptm_counts.ContainsKey(kk)) { continue; }
#     else { ptm_counts.Add(kk, p.OneBasedPossibleLocalizedModifications.Sum(kv => kv.Value.Count)); }
# }
# File.WriteAllLines(wd + @"\output\all_ptmcounts.csv", (new List<string> { "gene,ptm_count" }).Concat(ptm_counts.Select(kv => $"{kv.Key},{kv.Value.ToString()}").ToArray()));


#%% Import UniProt PTM information for the proteins
# Idea: Non-transcriptionally regulated proteins with single cell variability
# may have more PTMs that regulate their variability
# Execution: build a table of ptm counts for each protein, 
# compare transcriptionally regulated to non-transcriptionally regulated
# Output: means and p-values for the two categories
all_ptmcounts = pd.read_csv("output/all_ptmcounts.csv")
ccd_t_uniprot = all_ptmcounts[np.isin(all_ptmcounts.gene, dianaccd_transcript_regulated)]
ccd_n_uniprot = all_ptmcounts[np.isin(all_ptmcounts.gene, dianaccd_nontranscript_regulated)]
print(f"mean number of ptms for transcriptionally regulated: {np.mean(ccd_t_uniprot['ptm_count'])}")
print(f"mean number of ptms for non-transcriptionally regulated: {np.mean(ccd_n_uniprot['ptm_count'])}")
t, p = scipy.stats.ttest_ind(ccd_t_uniprot["ptm_count"], ccd_n_uniprot["ptm_count"])
print(f"one-sided t-test for that non-transcriptionally regulated: {2*p}")


#%% Import and analyze phosphosite plus PTM information from table
# Idea: Determine if there are significantly higher counts of PTMs for the genes that 
# AREN'T regulated at the transcript level.
# Execution: check on the sequence annotations; analyze phosphosite dataframe
# Output: statistics summary comparing sites from those sets and all proteins

# 1) Import phosphosite plus PTM information from sequences
# human_proteins = []
# with open("C:\\Users\\antho\\Documents\\Projects\\PhosphositePlus\\Phosphosite_PTM_seq.fasta") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         if "human" in record.id: human_proteins.append(record)
# cp1 = sum([1 for p in human_proteins if any([c.islower() for c in p.seq])])
# print(f"number of {len(human_proteins)} total human proteins with a PTM in PhosphositePlus: {cp1}")
# sites = ''.join(set(''.join(set([str(r.seq) for r in human_proteins]))))
# print("sites with modifications: " + ''.join([c for c in sites if c.islower()]))

# 2) Import phosphosite dataframe
phosphosites = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\PhosphositePlus\\Phosphorylation_site_dataset", delimiter="\t")
phs_human = phosphosites[phosphosites.ORGANISM == "human"]
allphosphocts = phs_human.groupby("GENE").size()
ccd_t_psp = allphosphocts[np.isin(allphosphocts.index, dianaccd_transcript_regulated)]
ccd_n_psp = allphosphocts[np.isin(allphosphocts.index, dianaccd_nontranscript_regulated)]
print(f"mean number of phosphosites for all proteins: {np.mean(allphosphocts)}")
print(f"mean number of phosphosites for transcriptionally regulated: {np.mean(ccd_t_psp)}")
print(f"mean number of phosphosites for non-transcriptionally regulated: {np.mean(ccd_n_psp)}")
t, p = scipy.stats.ttest_ind(ccd_t_psp, ccd_n_psp)
print(f"one-sided t-test for same means of transcript and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.ttest_ind(allphosphocts, ccd_n_psp)
print(f"one-sided t-test for same means of all and non-transcriptionally regulated: {2*p}")

print()
print(f"median number of phosphosites for all proteins: {np.median(allphosphocts)}")
print(f"mean number of phosphosites for transcriptionally regulated: {np.median(ccd_t_psp)}")
print(f"mean number of phosphosites for non-transcriptionally regulated: {np.median(ccd_n_psp)}")
k, p = scipy.stats.kruskal(ccd_t_psp, ccd_n_psp)
print(f"one-sided Wilcoxon test for same median of transcript and non-transcriptionally regulated: {2*p}")
k, p = scipy.stats.kruskal(allphosphocts, ccd_n_psp)
print(f"one-sided Wilcoxon test for same median all and non-transcriptionally regulated: {2*p}")


#%% [markdown]
# # Conclusions for PTM count comparisons
# Among all proteins, the average number of PTMs per protein is 0.46;
# among transcriptionally regulated CCD proteins, there is an average of 3.3 PTMs per protein;
# among non-transcriptionally regulated CCD proteins, the average is 0.59 PTMs per protein.
# This difference between these two types of CCD proteins is significant, but not in the way we expected.
#
# This is recapitulated in both UniProt and phosphosite plus. In phosphosite plus, the median number
# of phosphosites is 7 for all proteins, 34.5 phosphosites for transcriptionally regulated CCD proteins,
# and 11 for non-transcriptionally regulated CCD proteins. These differences are significant.


#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in the protein group results
filenameCommon = 'C:\\Users\\antho\\Box\ProjectData\\Variability\\U2OS_gptmd_search\\2019-06-19-15-54-09_CommonPtmsWithOccupancy\\Task1-SearchTask\\AllProteinGroups.tsv'
filenameCommonAndLessCommon = 'C:\\Users\\antho\\Box\ProjectData\\Variability\\U2OS_gptmd_search\\U2OSGptmdSearchAllProteinGroups.tsv'
file = pd.read_csv(filenameCommon, sep="\t", index_col=False)
targets=file[(file["Protein Decoy/Contaminant/Target"] == "T") & (file["Protein QValue"] <= 0.01)]
modifiedProteins = targets[targets["Sequence Coverage with Mods"].str.replace("[","") != targets["Sequence Coverage with Mods"]]
modifications = [re.findall('\[.*?\]',s) for s in modifiedProteins["Sequence Coverage with Mods"]]
unique_mods = set([item for sublist in modifications for item in sublist])

genes = list(modifiedProteins["Gene"])
seqmods = list(modifiedProteins["Sequence Coverage with Mods"])
seqs = list(modifiedProteins["Sequence Coverage"])
coverage = list(modifiedProteins["Sequence Coverage %"])
genemods = {}
blacklist = ["oxidation", "deamidation", "ammonia loss", "water loss", "carbamyl", "carbamidomethyl", # artifacts
    "zinc", "fe[i", "magnesium", "cu[i", # metals
    "sodium", "potassium", "calcium"] # counterions
    
for idx in range(len(genes)):
    genesplit = str(genes[idx]).split("|")
    seqmod = seqmods[idx].split("|")
    seq = seqs[idx].split("|")
    cov = coverage[idx].split("|")
    for idxx in range(len(genesplit)):
        mods = [m for m in re.findall('\[.*?\]', seqmod[idxx]) if not any(mod in m.lower() for mod in blacklist)]
        effective_length = float(len(seq[idxx])) * float(cov[idxx][:-1]) / 100
        # mods_per_eff_base = math.log10(float(len(mods) + 1) / float(effective_length))
        mods_per_eff_base = float(len(mods)) / float(effective_length)
        genegene = genesplit[idxx]
        if genegene in genemods: 
            genemods[genegene][0].append(mods_per_eff_base)
            genemods[genegene][1].extend(mods)
        else: genemods[genegene] = ([mods_per_eff_base], mods)

print(f"{str(len(targets))} proteins")
print(f"{str(len(modifiedProteins))} modified proteins ({str(round(float(len(modifiedProteins))/float(len(targets))*100,2))}%)")
print(f"{str(len([g for g in genemods.keys() if g in allccd_transcript_regulated]))}: number of all transcript regulated CCD genes of {len(allccd_transcript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in dianaccd_transcript_regulated]))}: number of transcript regulated CCD genes from Diana's study of {len(dianaccd_transcript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in dianaccd_nontranscript_regulated]))}: number of non-transcript regulated CCD genes from Diana's study of {len(dianaccd_nontranscript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in genes_analyzed]))}: number of proteins of {len(genes_analyzed)} detected.")

unambigenes = []
modcts = []
modsss = []
for gene in genemods.keys():
    unambigenes.append(gene)
    modcts.append(np.median(genemods[gene][0]))
    modsss.append(", ".join(genemods[gene][1]))

df = pd.DataFrame({"gene" : unambigenes, "modcts" : modcts, "modlist" : modsss})

all_modctsdf = df[np.isin(df["gene"], genes_analyzed)]
all_modctsdf.to_csv("output/ProteinModCounts.csv")
all_modcts = all_modctsdf["modcts"]
all_modcts.hist()
plt.title("All, modcts")
plt.show()
plt.savefig("")
plt.close()

ccd_at_modctsdf = df[np.isin(df["gene"], allccd_transcript_regulated)]
ccd_at_modctsdf.to_csv("output/ProteinModCountsTransRegCcd.csv")
ccd_at_modcts = ccd_at_modctsdf["modcts"]
ccd_at_modcts.hist()
plt.title("All CCD Transcript Reg, modcts")
plt.show()
plt.close()

ccd_t_modctsdf = df[np.isin(df["gene"], dianaccd_transcript_regulated)]
ccd_t_modctsdf.to_csv("output/ProteinModCountsTransRegCcd.csv")
ccd_t_modcts = ccd_t_modctsdf["modcts"]
ccd_t_modcts.hist()
plt.title("Diana's CCD Transcript Reg, modcts")
plt.show()
plt.close()

ccd_n_modctsdf = df[np.isin(df["gene"], dianaccd_nontranscript_regulated)]
ccd_n_modctsdf.to_csv("output/ProteinModCountsNonTransRegCcd.csv")
ccd_n_modcts = ccd_n_modctsdf["modcts"]
ccd_n_modcts.hist()
plt.title("Diana's CCD Non-Transcript Reg, modcts")
plt.show()
plt.close()

# The t-test isn't good here because these counts don't form a normal distribution... it works with the log(x+1 / len) distribution
print(f"mean number of mods / seq length for all proteins: {np.mean(all_modcts)}")
print(f"mean number of mods / seq length for all transcriptionally regulated CCD: {np.mean(ccd_at_modcts)}")
print(f"mean number of mods / seq length for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modcts)}")
print(f"mean number of mods / seq length for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modcts)}")
# t, p = scipy.stats.ttest_ind(ccd_t_modcts, ccd_n_modcts)
# print(f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}")
# t, p = scipy.stats.ttest_ind(all_modcts, ccd_t_modcts)
# print(f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}")
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
t, p = scipy.stats.kruskal(all_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_t_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}")
t, p = scipy.stats.kruskal(all_modcts, ccd_at_modcts, ccd_n_modcts)
print(f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}")
print()

# Generate a histogram of effective mod counts with bins normalized to 1'''
def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

# bar histogram
bins=np.histogram(np.hstack((all_modcts, ccd_t_modcts, ccd_n_modcts)), bins=40)[1] #get the bin edges
plt.hist(all_modcts, bins=bins, weights=weights(all_modcts), histtype="bar", alpha=0.4, label="All Proteins")
plt.hist(ccd_t_modcts, bins=bins, weights=weights(ccd_t_modcts), histtype="bar", alpha=0.4, label="Transcript Reg. CCD")
plt.hist(ccd_n_modcts, bins=bins, weights=weights(ccd_n_modcts), histtype="bar", alpha=0.5, label="Non-Transcript Reg. CCD")
plt.legend(loc="upper right")
plt.xlabel("Modifications per Covered Base")
plt.savefig("figures/ModsPerCoveredBaseHistogram.png")
plt.show()
plt.close()

# line histogram
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

# Generate boxplots for effective mod counts
mmmm = np.concatenate((all_modcts, ccd_t_modcts, ccd_at_modcts, ccd_n_modcts))
cccc = (["All"] * len(all_modcts))
cccc.extend(["Diana's Transcript CCD"] * len(ccd_t_modcts))
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

#%% Okay, next I do want to incorporate the occupancy to see if it makes this more robust

#%% [markdown]
# Bulk U2OS analysis results
# It turns out the significance went away when I filtered out an artifact I missed (iron adducts),
# and then after checking all transcript-regulated CCD genes, rather than just the 
# set in Diana's study, it all went away.
# 
# On the bright side, there are some really interesting modifications even in the 
# bulk sequencing: plenty of phosphorylations and even a nitrosylation on BUB1B that's
# not in UniProt.


#%%
