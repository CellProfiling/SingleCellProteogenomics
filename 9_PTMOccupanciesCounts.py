#%% Imports
from imports import *
from Bio import SeqIO
import sys
import re
import math
from methods_RNASeqData import read_counts_and_phases, qc_filtering
import seaborn as sbn

#%% Import the genes names we're analyzing
def ccd_gene_names(id_list_like):
    '''Convert gene ID list to gene name list'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    return gene_info[(gene_info["gene_id"].isin(id_list_like))]["name"]

def geneIdToHngc(id_list_like):
    '''Convert gene ID list to HNGC symbol if it exists'''
    gene_info = pd.read_csv("input/processed/python/ENSGToHGNC.csv", index_col=False, header=0)
    gene_info = gene_info[gene_info["hgnc_symbol"] != ""]
    return gene_info[(gene_info["gene_id"].isin(id_list_like))]["hgnc_symbol"]

def ccd_gene_lists(adata):
    '''Read in the published CCD genes / Diana's CCD / Non-CCD genes'''
    gene_info = pd.read_csv("input/processed/python/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    ccd_regev=pd.read_csv("input/processed/manual/ccd_regev.txt")    
    ccd=pd.read_csv("output/picklestxt/ccd_compartment_ensg.txt")#"input/processed/manual/ccd_genes.txt")
    nonccd=pd.read_csv("output/picklestxt/nonccd_compartment_ensg.txt")#("input/processed/manual/nonccd_genes.txt")
    ccd_regev_filtered = list(gene_info[(gene_info["name"].isin(ccd_regev["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    ccd_filtered = list(gene_info[(gene_info["name"].isin(ccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    nonccd_filtered = list(gene_info[(gene_info["name"].isin(nonccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    return ccd_regev_filtered, ccd_filtered, nonccd_filtered

ccdtranscript = np.load("output/pickles/ccdtranscript.npy", allow_pickle=True)
ccdprotein_transcript_regulated = np.load("output/pickles/ccdprotein_transcript_regulated.npy", allow_pickle=True)
ccdprotein_nontranscript_regulated = np.load("output/pickles/ccdprotein_nontranscript_regulated.npy", allow_pickle=True)
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])
bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures

wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)
nonccd_comp_ensg = wp_ensg[nonccd_comp]
ccd_comp_ensg = np.unique(wp_ensg[ccd_comp]) # just pseudotime
# ccd_comp_ensg = np.unique(np.concatenate((bioccd, wp_ensg[ccd_comp]))) # pseudotime and mitotic structures

# Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms"
biotype_to_use="protein_coding"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, False, biotype_to_use)
adata, phasesfilt = qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

isInterphaseCcd = np.isin(adata.var_names, ccd_comp_ensg)
nonccdtranscript = set(ccd_gene_names(adata.var_names[~ccdtranscript]))
ccdtranscript = set(ccd_gene_names(adata.var_names[ccdtranscript]))
ccdprotein_transcript_regulated = set(ccd_gene_names(adata.var_names[ccdprotein_transcript_regulated & isInterphaseCcd]))
ccdprotein_nontranscript_regulated = set(ccd_gene_names(adata.var_names[ccdprotein_nontranscript_regulated & isInterphaseCcd]))
diananonccd = set(ccd_gene_names(nonccd_filtered))
genes_analyzed = set(ccd_gene_names(genes_analyzed))
ccd_regev_filtered = set(ccd_gene_names(ccd_regev_filtered))
ccd_filtered = set(ccd_gene_names(ccd_filtered))
nonccdprotein = set(ccd_gene_names(nonccd_comp_ensg))
ccdprotein = set(ccd_gene_names(ccd_comp_ensg))

# Get the HGNC symbols for lists for GO analysis
bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
pd.DataFrame({"gene":list(set(geneIdToHngc(adata.var_names[np.load("output/pickles/ccdtranscript.npy", allow_pickle=True)])))}).to_csv("output/hgnc_ccdtranscript.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(geneIdToHngc(adata.var_names[np.load("output/pickles/ccdprotein_transcript_regulated.npy", allow_pickle=True)])))}).to_csv("output/hgnc_ccdprotein_transcript_regulated.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(geneIdToHngc(adata.var_names[np.load("output/pickles/ccdprotein_nontranscript_regulated.npy", allow_pickle=True)])))}).to_csv("output/hgnc_ccdprotein_nontranscript_regulated.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(geneIdToHngc(wp_ensg[~ccd_comp])))}).to_csv("output/hgnc_nonccdprotein.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(geneIdToHngc(np.concatenate((wp_ensg[ccd_comp], bioccd)))))}).to_csv("output/hgnc_ccdprotein.csv", index=False, header=False)

# Save the geneIds for each category
pd.DataFrame({"gene":list(set(adata.var_names[np.load("output/pickles/ccdtranscript.npy", allow_pickle=True)]))}).to_csv("output/ensg_ccdtranscript.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(adata.var_names[np.load("output/pickles/ccdprotein_transcript_regulated.npy", allow_pickle=True)]))}).to_csv("output/ensg_ccdprotein_transcript_regulated.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(adata.var_names[np.load("output/pickles/ccdprotein_nontranscript_regulated.npy", allow_pickle=True)]))}).to_csv("output/ensg_ccdprotein_nontranscript_regulated.csv", index=False, header=False)
pd.DataFrame({"gene":list(set(wp_ensg[nonccd_comp]))}).to_csv("output/ensg_nonccdprotein.csv", index=False, header=False)
pd.DataFrame({"gene": ccd_comp_ensg}).to_csv("output/ensg_ccdprotein.csv", index=False, header=False)


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
filenameCommon = 'C:\\Users\\antho\\Dropbox\\ProjectData\\Variability\\U2OS_gptmd_search\\' + search + '\\AllProteinGroups.tsv'
filenameCommonAndLessCommon = 'C:\\Users\\antho\\Dropbox\\ProjectData\\Variability\\U2OS_gptmd_search\\U2OSGptmdSearchAllProteinGroups.tsv'
file = pd.read_csv(filenameCommon, sep="\t", index_col=False)
targets=file[(file["Protein Decoy/Contaminant/Target"] == "T") & (file["Protein QValue"] <= 0.01)]
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
MIN_MOD_PEP = 1
MIN_TOT_PEP = 3

for idx in range(len(genes)):
    genesplit = str(genes[idx]).split("|")
    seqmod = seqmods[idx].split("|")
    seq = seqs[idx].split("|")
    modinfo = [x.strip(';') for x in str(modinfos[idx]).split("|")]
    for idxx, genegene in enumerate(genesplit):
        aatemp, modstemp, occupancytemp = [], [], []
        modpepts, totpepts = [], []

        if idxx < len(modinfo):
            modis = [m for m in modinfo[idxx].split(";") if m != "nan"]
            aatemp = [m[len("#aa"):m.index('[')] for m in modis]
            modstemp = [m[(m.index('[') + 1):m.index(",info")] for m in modis]
            modstemp = [re.sub("(Phospho[st]\w+)", "Phosphorylation", m) for m in modstemp] # name the phosphos consistently
            occupancytemp = [m[(m.index("occupancy=") + len("occupancy=")):(m.index("occupancy=") + len("occupancy=") + len("#.##"))] for m in modis]
            modpepts = [m[(m.index("(", m.index("occupancy=")) + 1):(m.index("/", m.index("occupancy=")))] for m in modis]
            totpepts = [m[(m.index("/", m.index("occupancy=")) + 1):(m.index(")", m.index("occupancy=")))] for m in modis]
        else: continue

        # mods = [m for m in re.findall('\[.*?\]', seqmod[idxx]) if not any(mod in m.lower() for mod in blacklist)]
        covered_fraction = float(sum([1 for c in seq[idxx] if c.islower()])) / float(len(seq[idxx]))
        effective_length = float(len(seq[idxx])) * covered_fraction
        # mods_per_eff_base = math.log10(float(len(modstemp) + 1) / float(effective_length))
        mods_per_eff_base = float(len(modstemp)) / float(effective_length)
        aanums, mods, occupancies, modpeptsss, totpeptsss = [],[],[],[],[]
        for oidx, occcc in enumerate(occupancytemp):
            isblacklisted = any(mod in modstemp[oidx].lower() for mod in blacklist)
            isonehitwonder = int(modpepts[oidx]) < MIN_MOD_PEP or int(totpepts[oidx]) < MIN_TOT_PEP
            if not "(" in occcc and not isblacklisted and not isonehitwonder: 
                aanums.append(aatemp[oidx])
                mods.append(modstemp[oidx])
                occupancies.append(float(occcc))
                modpeptsss.append(int(modpepts[oidx]))
                totpeptsss.append(int(totpepts[oidx]))
        if genegene not in genemods: 
            genemods[genegene] = ([mods_per_eff_base], occupancies, mods, 
                modpeptsss, totpeptsss, aanums)
        else:
            genemods[genegene][0].append(mods_per_eff_base)
            genemods[genegene][1].extend(occupancies)
            genemods[genegene][2].extend(mods)
            genemods[genegene][3].extend(modpeptsss)
            genemods[genegene][4].extend(totpeptsss)
            genemods[genegene][5].extend(aanums)

print(f"{str(len(targets))} proteins")
print(f"{str(len(modifiedProteins))} modified proteins ({str(round(float(len(modifiedProteins))/float(len(targets))*100,2))}%)")
print(f"{str(len([g for g in genemods.keys() if g in ccdtranscript]))}: number of all transcript regulated CCD genes of {len(ccdtranscript)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in ccdprotein_transcript_regulated]))}: number of transcript regulated CCD genes from Diana's study of {len(ccdprotein_transcript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in ccdprotein_nontranscript_regulated]))}: number of non-transcript regulated CCD genes from Diana's study of {len(ccdprotein_nontranscript_regulated)} detected.")
print(f"{str(len([g for g in genemods.keys() if g in genes_analyzed]))}: number of proteins of {len(genes_analyzed)} detected.")

unambigenes, modcts, avgocc, modsss = [], [], [], []
protmodgene, modmod, occocc, aanumnum = [], [], [], []
modpeps, totpeps = [], []
for gene in genemods.keys():
    unambigenes.append(gene)
    modcts.append(np.median(genemods[gene][0]))
    avgocc.append(np.median(genemods[gene][1]) if any(genemods[gene][1]) else 0)
    modsss.append(", ".join(genemods[gene][2]))
    for modidx, mod in enumerate(genemods[gene][2]):
        protmodgene.append(gene)
        modmod.append(mod)
        occocc.append(genemods[gene][1][modidx])
        modpeps.append(genemods[gene][3][modidx])
        totpeps.append(genemods[gene][4][modidx])
        aanumnum.append(genemods[gene][5][modidx])
        
df = pd.DataFrame({
    "gene" : unambigenes, 
    "modcts" : modcts, 
    "avgocc" : avgocc,
    "modlist" : modsss})
occdf = pd.DataFrame({
    "gene" : protmodgene,
    "modification" : modmod,
    "residue" : aanumnum,
    "occupancy" : occocc,
    "modpeps" : modpeps,
    "totpeps" : totpeps})


# all genes
all_modctsdf = df[df["gene"].isin(genes_analyzed)]
all_modcts = all_modctsdf["modcts"]
all_avgocc = all_modctsdf["avgocc"]#[all_modctsdf["avgocc"] != 0]
all_modctsoccdf = occdf[occdf["gene"].isin(genes_analyzed)]
all_modocc = all_modctsoccdf["occupancy"]

# regev
ccd_rt_modctsdf = df[df["gene"].isin(ccd_regev_filtered)]
ccd_rt_modcts = ccd_rt_modctsdf["modcts"]
ccd_rt_avgocc = ccd_rt_modctsdf["avgocc"]#[ccd_rt_modcts["avgocc"] != 0]
ccd_rt_modctsoccdf = occdf[occdf["gene"].isin(ccd_regev_filtered)]
ccd_rt_modocc = ccd_rt_modctsoccdf["occupancy"]

# transcript CCD (all)
ccd_at_modctsdf = df[df["gene"].isin(ccdtranscript)]
ccd_at_modcts = ccd_at_modctsdf["modcts"]
ccd_at_avgocc = ccd_at_modctsdf["avgocc"]#[ccd_at_modcts["avgocc"] != 0]
ccd_at_modctsoccdf = occdf[occdf["gene"].isin(ccdtranscript)]
ccd_at_modocc = ccd_at_modctsoccdf["occupancy"]

# transcript non-CCD (all)
ccd_nt_modctsdf = df[df["gene"].isin(nonccdtranscript)]
ccd_nt_modcts = ccd_nt_modctsdf["modcts"]
ccd_nt_avgocc = ccd_nt_modctsdf["avgocc"]#[ccd_nt_modcts["avgocc"] != 0]
ccd_nt_modctsoccdf = occdf[occdf["gene"].isin(nonccdtranscript)]
ccd_nt_modocc = ccd_nt_modctsoccdf["occupancy"]

# transcript reg. CCD proteins
ccd_t_modctsdf = df[df["gene"].isin(ccdprotein_transcript_regulated)]
ccd_t_modcts = ccd_t_modctsdf["modcts"]
ccd_t_avgocc = ccd_t_modctsdf["avgocc"]#[ccd_t_modcts["avgocc"] != 0]
ccd_t_modctsoccdf = occdf[occdf["gene"].isin(ccdprotein_transcript_regulated)]
ccd_t_modocc = ccd_t_modctsoccdf["occupancy"]

# non-transcript reg. CCD proteins
ccd_n_modctsdf = df[df["gene"].isin(ccdprotein_nontranscript_regulated)]
ccd_n_modcts = ccd_n_modctsdf["modcts"]
ccd_n_avgocc = ccd_n_modctsdf["avgocc"]#[ccd_n_modcts["avgocc"] != 0]
ccd_n_modctsoccdf = occdf[occdf["gene"].isin(ccdprotein_nontranscript_regulated)]
ccd_n_modocc = ccd_n_modctsoccdf["occupancy"]

# CCD proteins
ccd_modctsdf = df[df["gene"].isin(ccdprotein)]
ccd_modcts = ccd_modctsdf["modcts"]
ccd_avgocc = ccd_modctsdf["avgocc"]#[ccd_modcts["avgocc"] != 0]
ccd_modctsoccdf = occdf[occdf["gene"].isin(ccdprotein)]
ccd_modocc = ccd_modctsoccdf["occupancy"]

# non-CCD proteins
nonccd_modctsdf = df[df["gene"].isin(nonccdprotein)]
nonccd_modcts = nonccd_modctsdf["modcts"]
nonccd_avgocc = nonccd_modctsdf["avgocc"]#[nonccd_modcts["avgocc"] != 0]
nonccd_modctsoccdf = occdf[occdf["gene"].isin(nonccdprotein)]
nonccd_modocc = nonccd_modctsoccdf["occupancy"]


# Just number of proteins with mods
fom = "\n"
fom += "COUNTS OF MODIFIED PROTEINS\n"
fom += f"{str(len(all_modcts[all_modcts > 0]))} modified proteins ({str(round(float(len(all_modcts[all_modcts > 0]))/float(len(all_modcts))*100,2))}%) in all genes detected\n"
fom += f"{str(len(ccd_rt_modcts[ccd_rt_modcts > 0]))} modified proteins ({str(round(float(len(ccd_rt_modcts[ccd_rt_modcts > 0]))/float(len(ccd_rt_modcts))*100,2))}%) in Regev CCD\n"
fom += f"{str(len(ccd_at_modcts[ccd_at_modcts > 0]))} modified proteins ({str(round(float(len(ccd_at_modcts[ccd_at_modcts > 0]))/float(len(ccd_at_modcts))*100,2))}%) in all transcript CCD genes\n"
fom += f"{str(len(ccd_t_modcts[ccd_t_modcts > 0]))} modified proteins ({str(round(float(len(ccd_t_modcts[ccd_t_modcts > 0]))/float(len(ccd_t_modcts))*100,2))}%) in transcript regulated CCD proteins\n"
fom += f"{str(len(ccd_n_modcts[ccd_n_modcts > 0]))} modified proteins ({str(round(float(len(ccd_n_modcts[ccd_n_modcts > 0]))/float(len(ccd_n_modcts))*100,2))}%) in non-transcript regulated CCD proteins\n"
fom += f"{str(len(nonccd_modcts[nonccd_modcts > 0]))} modified proteins ({str(round(float(len(nonccd_modcts[nonccd_modcts > 0]))/float(len(nonccd_modcts))*100,2))}%) in non-CCD proteins\n"
fom += f"{str(len(ccd_modcts[ccd_modcts > 0]))} modified proteins ({str(round(float(len(ccd_modcts[ccd_modcts > 0]))/float(len(ccd_modcts))*100,2))}%) in CCD proteins\n"
fom += "\n"
fom += "COUNTS OF MODIFIED PROTEINS WITH MOD OCCUPANCY CUTOFFS\n"
fom += f"{str(len(set(all_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(all_modctsoccdf['gene'])))/float(len(all_modcts))*100,2))}%) in all genes detected\n"
fom += f"{str(len(set(ccd_rt_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_rt_modctsoccdf['gene'])))/float(len(ccd_rt_modcts))*100,2))}%) in Regev CCD\n"
fom += f"{str(len(set(ccd_at_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_at_modctsoccdf['gene'])))/float(len(ccd_at_modcts))*100,2))}%) in all transcript CCD genes\n"
fom += f"{str(len(set(ccd_t_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_t_modctsoccdf['gene'])))/float(len(ccd_t_modcts))*100,2))}%) in transcript regulated CCD proteins\n"
fom += f"{str(len(set(ccd_nt_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_nt_modctsoccdf['gene'])))/float(len(ccd_nt_modcts))*100,2))}%) non-CCD transcript proteins\n"
fom += f"{str(len(set(ccd_n_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_n_modctsoccdf['gene'])))/float(len(ccd_n_modcts))*100,2))}%) in non-transcript regulated CCD proteins\n"
fom += f"{str(len(set(ccd_n_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_n_modctsoccdf['gene'])))/float(len(nonccd_modcts))*100,2))}%) in non-CCD proteins\n"
fom += f"{str(len(set(ccd_modctsoccdf['gene'])))} modified proteins ({str(round(float(len(set(ccd_modctsoccdf['gene'])))/float(len(ccd_modcts))*100,2))}%) in CCD proteins\n"
fom += "\n"


# The t-test isn't good here because these counts don't form a normal distribution... it works with the log(x+1 / len) distribution
fom += "COUNT OF MODIFICATIONS PER COVERED RESIDUE PER PROTEIN\n"
fom += f"mean number of mods / seq length for all proteins: {np.mean(all_modcts)}\n"
fom += f"mean number of mods / seq length for all transcriptionally regulated CCD: {np.mean(ccd_at_modcts)}\n"
fom += f"mean number of mods / seq length for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modcts)}\n"
fom += f"mean number of mods / seq length for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modcts)}\n"
t, p = scipy.stats.ttest_ind(ccd_t_modcts, ccd_n_modcts)
fom += f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}\n"
t, p = scipy.stats.ttest_ind(all_modcts, ccd_t_modcts)
fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
fom += "\n"

fom += f"median number of mods / seq length for all proteins: {np.median(all_modcts)}\n"
fom += f"median number of mods / seq length for all transcriptionally regulated CCD: {np.median(ccd_at_modcts)}\n"
fom += f"median number of mods / seq length for Diana's transcriptionally regulated CCD: {np.median(ccd_t_modcts)}\n"
fom += f"median number of mods / seq length for Diana's non-transcriptionally regulated: {np.median(ccd_n_modcts)}\n"
t, p = scipy.stats.kruskal(ccd_at_modcts, ccd_n_modcts)
fom += f"one-sided kruskal for same medians of all transcriptionally reg. and non-transcriptionally regulated: {2*p}\n"
t, p = scipy.stats.kruskal(ccd_t_modcts, ccd_n_modcts)
fom += f"one-sided kruskal for same medians of Diana's transcriptionally reg. and non-transcriptionally regulated: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_t_modcts)
fom += f"one-sided kruskal for same medians of all genes and Diana's transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_at_modcts)
fom += f"one-sided kruskal for same medians of all genes and all transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_rt_modcts)
fom += f"one-sided kruskal for same medians of all genes and regev transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_n_modcts)
fom += f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_t_modcts, ccd_n_modcts)
fom += f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}\n"
t, p = scipy.stats.kruskal(all_modcts, ccd_at_modcts, ccd_n_modcts)
fom += f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}\n"
fom += "\n"

# The t-test isn't good here because these counts don't form a normal distribution... it works with the log(x+1 / len) distribution
fom += "OCCUPANCY PER MODIFIED BASE\n"
fom += "Note: I excluded artifact modifications and required there to be\n"
fom += f"at least {MIN_MOD_PEP} peptides with the modification to count it and\n"
fom += f"at least {MIN_TOT_PEP} peptides covering the base total.\n"
fom += "\n"
fom += f"mean occupancy per modified residue for all proteins: {np.mean(all_modocc)}\n"
fom += f"mean occupancy per modified residue for all transcriptionally regulated CCD: {np.mean(ccd_at_modocc)}\n"
fom += f"mean occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modocc)}\n"
fom += f"mean occupancy per modified residue for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modocc)}\n"
t, p = scipy.stats.ttest_ind(ccd_t_modocc, ccd_n_modocc)
fom += f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}\n"
t, p = scipy.stats.ttest_ind(all_modocc, ccd_t_modocc)
fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
t, p = scipy.stats.ttest_ind(all_modocc, ccd_t_modocc)
fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
t, p = scipy.stats.ttest_ind(all_modocc, ccd_t_modocc)
fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
fom += "\n"
fom += f"median occupancy per modified residue for all proteins: {np.median(all_modocc)}\n"
fom += f"median occupancy per modified residue for all transcriptionally regulated CCD: {np.median(ccd_at_modocc)}\n"
fom += f"median occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.median(ccd_t_modocc)}\n"
fom += f"median occupancy per modified residue for Diana's non-transcriptionally regulated: {np.median(ccd_n_modocc)}\n"
t, p = scipy.stats.kruskal(ccd_at_modocc, ccd_n_modocc)
fom += f"one-sided kruskal for same medians of all transcriptionally reg. and non-transcript regulated: {2*p}\n"
t, p = scipy.stats.kruskal(ccd_t_modocc, ccd_n_modocc)
fom += f"one-sided kruskal for same medians of Diana's transcriptionally reg. CCD and non-transcript regulated: {2*p}\n"
t, p = scipy.stats.kruskal(np.concatenate((ccd_t_modocc, ccd_n_modocc)), nonccd_modocc)
fom += f"one-sided kruskal for same medians of Diana's CCD and non-CCD regulated: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_t_modocc)
fom += f"one-sided kruskal for same medians of all genes and Diana's transcriptionally reg CCD: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_at_modocc)
fom += f"one-sided kruskal for same medians of all genes and all transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_rt_modocc)
fom += f"one-sided kruskal for same medians of all genes and regev transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_n_modocc)
fom += f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, nonccd_modocc)
fom += f"one-sided kruskal for same medians of all genes and non-CCD: {2*p}\n"
t, p = scipy.stats.kruskal(ccd_modocc, nonccd_modocc)
fom += f"one-sided kruskal for same medians of CCD genes and non-CCD: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_t_modocc, ccd_n_modocc)
fom += f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}\n"
t, p = scipy.stats.kruskal(all_modocc, ccd_at_modocc, ccd_n_modocc)
fom += f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}\n"
fom += "\n"

# save and print results
print(fom)
with open("output/modificationCountsOccResults.txt", 'w') as file:
    file.write(fom)
    
# Generate a histogram of effective mod counts with bins normalized to 1'''
def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

# Generate boxplots for effective mod counts
plt.figure(figsize=(10,10))
mmmm = np.concatenate((all_modcts, ccd_t_modcts, ccd_rt_modcts, ccd_at_modcts, ccd_n_modcts, nonccd_modcts))
cccc = (["All"] * len(all_modcts))
cccc.extend(["Diana's Transcript CCD"] * len(ccd_t_modcts))
cccc.extend(["Regev Transcript CCD"] * len(ccd_rt_modcts))
cccc.extend(["All Transcript CCD"] * len(ccd_at_modcts))
cccc.extend(["Non-Transc CCD"] * len(ccd_n_modcts))
cccc.extend(["Non-CCD"] * len(nonccd_modcts))
boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=False)
boxplot.set_xlabel("Protein Set", size=36,fontname='Arial')
boxplot.set_ylabel("Modifications per Covered Base", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major")#, labelsize=16)
plt.title("")
plt.savefig("figures/ModsPerCoveredBaseBoxplot.png")
plt.show()
plt.close()

# Generate boxplots for mod occupancies
plt.figure(figsize=(10,10))
mmmm = np.concatenate((all_modocc, ccd_t_modocc, ccd_rt_modocc, ccd_at_modocc, ccd_n_modocc, nonccd_modocc))
cccc = (["All"] * len(all_modocc))
cccc.extend(["Diana's\nTranscript\nCCD"] * len(ccd_t_modocc))
cccc.extend(["Regev\nTranscript\nCCD"] * len(ccd_rt_modocc))
cccc.extend(["All\nTranscript\nCCD"] * len(ccd_at_modocc))
cccc.extend(["Diana's\nNon-Transc\nCCD"] * len(ccd_n_modocc))
cccc.extend(["Non-\nCCD"] * len(nonccd_modocc))
boxplot = sbn.boxplot(x=cccc, y=mmmm)#, showfliers=False)
boxplot.set_xlabel("Protein Set", size=36,fontname='Arial')
boxplot.set_ylabel("Occupancy per Modified Residue", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/ModsOccupancyBoxplot.pdf")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
mmmm = np.concatenate((all_modocc, ccd_at_modocc, ccd_nt_modocc, ccd_t_modocc, ccd_n_modocc, ccd_modocc, nonccd_modocc))
cccc = (["All"] * len(all_modocc))
cccc.extend(["All\nTranscript\nCCD"] * len(ccd_at_modocc))
cccc.extend(["All\nTranscript\nNon-CCD"] * len(ccd_nt_modocc))
cccc.extend(["Transcript\nReg. CCD"] * len(ccd_t_modocc))
cccc.extend(["Non-\nTranscript\nReg. CCD"] * len(ccd_n_modocc))
cccc.extend(["CCD"] * len(ccd_modocc))
cccc.extend(["Non-\nCCD"] * len(nonccd_modocc))
boxplot = sbn.boxplot(x=cccc, y=mmmm)#, showfliers=False)
boxplot = sbn.stripplot(x=cccc, y=mmmm, alpha=0.3, color=".3", size=5, jitter=0.25)#, showfliers=False)
boxplot.set_ylabel("Occupancy per Modified Residue", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
boxplot.set(ylim=(-0.01,1.01))
plt.savefig("figures/ModsOccupancyBoxplotSelected.pdf")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
mmmm = np.concatenate((ccd_at_modocc, ccd_nt_modocc, all_modocc, ccd_modocc, nonccd_modocc))
cccc = (["All"] * len(all_modocc))
cccc.extend(["All\nTranscript\nCCD"] * len(ccd_at_modocc))
cccc.extend(["All\nTranscript\nNon-CCD"] * len(ccd_nt_modocc))
# cccc.extend(["Transcript\nReg. CCD"] * len(ccd_t_modocc))
# cccc.extend(["Non-\nTranscript\nReg. CCD"] * len(ccd_n_modocc))
cccc.extend(["CCD"] * len(ccd_modocc))
cccc.extend(["Non-\nCCD"] * len(nonccd_modocc))
boxplot = sbn.boxplot(x=cccc, y=mmmm)#, showfliers=False)
boxplot = sbn.stripplot(x=cccc, y=mmmm, alpha=0.2, color=".3", size=7, jitter=0.25)#, showfliers=False)
boxplot.set_ylabel("Occupancy per Modified Residue", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
boxplot.set(ylim=(-0.01,1.01))
plt.savefig("figures/ModsOccupancyBoxplotSelected2.pdf")
plt.show()
plt.close()

# plt.figure(figsize=(10,10))
# mmmm = np.concatenate((all_avgocc, ccd_at_avgocc, ccd_nt_avgocc, ccd_t_avgocc, ccd_n_avgocc, ccd_avgocc, nonccd_avgocc))
# cccc = (["All"] * len(all_avgocc))
# cccc.extend(["All\nTranscript\nCCD"] * len(ccd_at_avgocc))
# cccc.extend(["All\nTranscript\nNon-CCD"] * len(ccd_nt_avgocc))
# cccc.extend(["Transcript\nReg. CCD"] * len(ccd_t_avgocc))
# cccc.extend(["Non-\nTranscript\nReg. CCD"] * len(ccd_n_avgocc))
# cccc.extend(["CCD"] * len(ccd_avgocc))
# cccc.extend(["Non-\nCCD"] * len(nonccd_avgocc))
# boxplot = sbn.violinplot(x=cccc, y=mmmm)#, showfliers=False)
# boxplot = sbn.stripplot(x=cccc, y=mmmm)#, showfliers=False)
# boxplot.set_ylabel("Avg Occupancy per Modified Residue per Protein", size=36,fontname='Arial')
# boxplot.tick_params(axis="both", which="major", labelsize=16)
# boxplot.set(ylim=(0,1))
# plt.savefig("figures/ModsAvgOccupancyBoxplotSelected.pdf")
# plt.show()
# plt.close()


#%% Pickle the results
pd.DataFrame({
    "all_modocc":all_modocc,
    "all_modocc_genes":all_modctsoccdf["gene"],
    "ccd_t_modocc":ccd_t_modocc,
    "ccd_t_modocc_genes":ccd_t_modctsoccdf["gene"],
    "ccd_rt_modocc":ccd_rt_modocc,
    "ccd_rt_modocc_genes":ccd_rt_modctsoccdf["gene"],
    "ccd_at_modocc":ccd_at_modocc,
    "ccd_at_modocc_genes":ccd_at_modctsoccdf["gene"],
    "ccd_n_modocc":ccd_n_modocc,
    "ccd_n_modocc_genes":ccd_n_modctsoccdf["gene"],
    "nonccd_modocc":nonccd_modocc,
    "nonccd_modocc_genes":nonccd_modctsoccdf["gene"],
}).to_pickle("output/modocc.pkl")


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



#%%


#%%
