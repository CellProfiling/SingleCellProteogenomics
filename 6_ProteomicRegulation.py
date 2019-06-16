#%% Imports
from imports import *
from Bio import SeqIO
# import requests
# import sys

#%% Import the genes names we're analyzing
ccd_transcript_regulated = np.array(pd.read_csv("output/ccd_transcript_regulated.csv")["gene"])
ccd_nontranscript_regulated = np.array(pd.read_csv("output/ccd_nontranscript_regulated.csv")["gene"])
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])

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
ccd_t_uniprot = all_ptmcounts[np.isin(all_ptmcounts.gene, ccd_transcript_regulated)]
ccd_n_uniprot = all_ptmcounts[np.isin(all_ptmcounts.gene, ccd_nontranscript_regulated)]
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
ccd_t_psp = allphosphocts[np.isin(allphosphocts.index, ccd_transcript_regulated)]
ccd_n_psp = allphosphocts[np.isin(allphosphocts.index, ccd_nontranscript_regulated)]
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

#%% Let's take a look at protein stability
# Idea: Is there a difference in protein stability or turnover for the proteins that are
# transcriptionally regulated CCD vs non-transcriptionally regulated?
# Execution: 
# Output:

# Using the categories
def plot_protein_stabilities(filename, title, splitname):
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    df1 = df["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    df2 = ccd_t_stab["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    df3 = ccd_n_stab["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    plt.plot(df1.index, df1 / np.array(df1).sum(axis=0, keepdims=True), label="All Proteins")
    plt.plot(df2.index, df2 / np.array(df2).sum(axis=0, keepdims=True), label="Transcript Reg. CCD")
    plt.plot(df3.index, df3 / np.array(df3).sum(axis=0, keepdims=True), label="Non-Transcript Reg. CCD")
    plt.title(title)
    plt.legend(loc="upper right")
    # plt.show()
    # plt.close()

plt.figure(figsize=(15,15))
plt.subplot(331)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
plt.subplot(332)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
plt.subplot(334)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
plt.subplot(335)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
plt.subplot(337)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
plt.subplot(338)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
plt.subplot(339)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
plt.savefig("figures/ProteinStabilities.png")
plt.show()
plt.close()

#%% [markdown]
# The stability classes don't show the greatest consistency between replicates
# or between cell lines. The melting points are a continuous category, so I might
# have more luck with it. It'll also omit the not-stable and non-melter categories,
# which seem like a bit of a crap shoot above

#%% Using melting temperatures (histograms)
# Idea: See if there is better consistency between replicates with the melting
#       points over the stability classes
# Execution: use melting points from data frames
# Output: histograms plots for the replicates and aggregate
#       statistical test for whether the transcript reg and non-transcript reg CCD
#       are different.

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))
    
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]

def add_temps(filename, title, splitname, merged):
    '''Adds melting temperature measurements from supp info files to lists'''
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    if len(merged) == 0: merged = df
    else: pd.merge(merged, df, on="ProteinName", suffixes=("", title))
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    all_temps.extend(df[pd.notna(df["Melting point [°C]"])]["Melting point [°C]"])
    transcript_reg.extend(ccd_t_stab[pd.notna(ccd_t_stab["Melting point [°C]"])]["Melting point [°C]"])
    nontranscr_reg.extend(ccd_n_stab[pd.notna(ccd_n_stab["Melting point [°C]"])]["Melting point [°C]"])

def stat_tests(title):
    print(f"{title} statistical tests")
    print("Testing whether transcript reg. CCD melting points are different than non-transcript reg. CCD")
    print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
    print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
    print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
    print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(transcript_reg, nontranscr_reg))
    print(scipy.stats.kruskal(transcript_reg, nontranscr_reg))
    print()
    print("Testing whether non-transcript reg. CCD is different than all proteins")
    print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
    print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
    print(f"{np.median(all_temps)}: All Proteins (median)")
    print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(all_temps, nontranscr_reg))
    print(scipy.stats.kruskal(all_temps, nontranscr_reg))
    print()
    print("Testing whether transcript reg. CCD is different than all proteins")
    print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
    print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
    print(f"{np.median(all_temps)}: All Proteins (median)")
    print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(all_temps, transcript_reg))
    print(scipy.stats.kruskal(all_temps, transcript_reg))
    print()

def temp_hist(title):
    '''Generates a histogram of melting points with bins normalized to 1'''
    bins=np.histogram(np.hstack((all_temps, transcript_reg, nontranscr_reg)), bins=40)[1] #get the bin edges
    plt.hist(all_temps, bins=bins, weights=weights(all_temps), alpha=0.5, label="All Proteins")
    plt.hist(transcript_reg, bins=bins, weights=weights(transcript_reg), alpha=0.5, label="Transcript Reg. CCD")
    plt.hist(nontranscr_reg, bins=bins, weights=weights(nontranscr_reg), alpha=0.6, label="Non-Transcript Reg. CCD")
    plt.legend(loc="upper right")
    plt.xlabel("Melting Point (°C)")
    plt.title(title)
    stat_tests(title)


# Individual histograms
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
plt.figure(figsize=(15,15))
plt.subplot(331)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True, merged)
temp_hist("A549_R1")
plt.subplot(332)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True, merged)
temp_hist("A549_R2")
plt.subplot(334)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True, merged)
temp_hist("HEK293_R1")
plt.subplot(335)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True, merged)
temp_hist("HEK293_R2")
plt.subplot(337)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False, merged)
temp_hist("HepG2_R1")
plt.subplot(338)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False, merged)
temp_hist("HepG2_R2")
plt.subplot(339)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False, merged)
temp_hist("HepG2_R3")
plt.savefig("figures/ProteinMeltingPointsIndivid.png")
plt.show()
plt.close()

# Aggregate histogram
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False, merged)
temp_hist("Aggregated Melting Points")
plt.savefig("figures/ProteinMeltingPoints.png")
plt.show()
plt.close()

# Statistical tests
print("Testing whether transcript reg. CCD melting points are different than non-transcript reg. CCD")
print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(transcript_reg, nontranscr_reg))
print(scipy.stats.kruskal(transcript_reg, nontranscr_reg))
print()
print("Testing whether non-transcript reg. CCD is different than all proteins")
print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
print(f"{np.median(all_temps)}: All Proteins (median)")
print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(all_temps, nontranscr_reg))
print(scipy.stats.kruskal(all_temps, nontranscr_reg))
print()
print("Testing whether transcript reg. CCD is different than all proteins")
print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
print(f"{np.median(all_temps)}: All Proteins (median)")
print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(all_temps, transcript_reg))
print(scipy.stats.kruskal(all_temps, transcript_reg))

#%% [markdown]
# The replicates do look better with the melting points, and the difference is
# significant between transcript reg CCD and non-transcript reg CCD. However, the
# non-transcript regulated CCD proteins are more similar to all proteins on average,
# and that difference is not significant.

#%% Using melting temperatures (histograms, combined protein measurements)
# Idea: The above conclusions are based on tests that require independent observations,
#       but some proteins measurements are represented more than once (max of 12 times)
# Execution: use melting points from data frames; combine protein measurements
# Output: histograms for the replicates and aggregate
#       statistical test for whether the transcript reg and non-transcript reg CCD
#       are different.

atp, trp, ntp = [], [], []
all_temps, transcript_reg, nontranscr_reg = [],[],[]

def add_temps_and_names(filename, title, splitname):
    '''Adds melting temperature measurements from supp info files to lists'''
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    notna = pd.notna(df["Melting point [°C]"])
    all_temps.extend(df[notna]["Melting point [°C]"])
    transcript_reg.extend(ccd_t_stab[notna]["Melting point [°C]"])
    nontranscr_reg.extend(ccd_n_stab[notna]["Melting point [°C]"])

    atp.extend(df[notna]["ProteinName"])
    trp.extend(ccd_t_stab[notna]["ProteinName"])
    ntp.extend(ccd_n_stab[notna]["ProteinName"])

def avg_prot_temps(temps, proteinNames):
    '''Finds median of the measurements from the same protein'''
    aaa = []
    for i in range(len(temps)):
        df = pd.DataFrame({"name" : proteinNames[i], "temps": temps[i]})
        aaa.append(list(df.groupby("name")["temps"].median()))
    return aaa

# Individual histograms
plt.figure(figsize=(15,15))
plt.subplot(331)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("A549_R1")
plt.subplot(332)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("A549_R2")
plt.subplot(334)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HEK293_R1")
plt.subplot(335)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
temp_hist("HEK293_R2")
plt.subplot(337)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R1")
plt.subplot(338)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R2")
plt.subplot(339)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R3")
plt.savefig("figures/ProteinMeltingPointsMedianedIndivid.png")
plt.show()
plt.close()

# Aggregate histogram
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("Aggregated Melting Points")
plt.savefig("figures/ProteinMeltingPointsMedianed.png")
plt.show()
plt.close()

#%% [markdown]
# The replicates do look better with the melting points, and the difference is
# significant between transcript reg CCD and non-transcript reg CCD. However, the
# non-transcript regulated CCD proteins are more similar to all proteins on average,
# and that difference is not significant.

#%% Making boxplots and violin plots for the melting point data
# Idea: trying out different visualizations for the above analysis
# (boxplots work; violin not so much yet)

# data = [transcript_reg, nontranscr_reg, all_temps]
# labels=["Transcript Reg.\nCCD", "Non-Transcript Reg.\nCCD", "All Proteins"]
# plt.boxplot(data, labels=labels)
# plt.savefig("figures/ProteinMeltingPointsBoxplot.png")
# plt.show()
# plt.close()

# def adjacent_values(vals, q1, q3):
#     upper_adjacent_value = q3 + (q3 - q1) * 1.5
#     upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

#     lower_adjacent_value = q1 - (q3 - q1) * 1.5
#     lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
#     return lower_adjacent_value, upper_adjacent_value

# def set_axis_style(ax, labels):
#     ax.get_xaxis().set_tick_params(direction='out')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.set_xticks(np.arange(1, len(labels) + 1))
#     ax.set_xticklabels(labels)
#     ax.set_xlim(0.25, len(labels) + 0.75)
#     ax.set_xlabel('Sample name')

# fig, ax1 = plt.subplots(1, 1, figsize=(9,4))
# parts = ax1.violinplot(data, showmeans=True, showmedians=False, showextrema=False)
# quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75])
# whiskers = np.array([adjacent_values(sorted_array, q1, q3) for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
# whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]
# set_axis_style(ax1, labels)
# plt.show()

#%%
