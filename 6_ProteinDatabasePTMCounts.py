#%% Imports
from imports import *
from Bio import SeqIO
import requests
import sys
import re
import math
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names


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


#%% Meh that might not help.
# Also, pulling in CORUM and looking for differences of PTMs within
# the modules of protein complexes.
# Idea: count number of modifications on genes and per base on the complexes
corum = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\ProteinClusters\\allComplexes.txt", sep="\t")
human_corum = corum[corum.organism == "Human"]
human_corum_genes = list(human_corum["subunits(Gene name)"])

# split the gene lists

# sum/average the mods counts for all of the genes

# plot the modcount plots and list the number of clusters with mods

#%%


#%%
