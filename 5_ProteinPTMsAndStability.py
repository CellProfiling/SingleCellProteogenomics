# -*- coding: utf-8 -*-
"""
Investigation of the properties of proteins with different cell cycle regulation using PTMs and stability measurements
-  PTMs are observed using previously published bulk and phospho-enriched mass spectrometry (MS) proteomic data
-  Differences in PTM regulation is inferred using PTM occupancy for each PTM site
-  Protein stability was measured by MS thermal profiling in an external study
-  Differences in thermal shifts indicate different stabilities and propensity for unfolding

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

#%% Imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, Loaders, FucciCellCycle, RNADataPreparation, PTMAnalysis, ProteinStabilityAnalysis
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

fucci = FucciCellCycle.FucciCellCycle()

#%% Import the genes names we're analyzing
# Read in RNA-Seq data again and the CCD gene lists
valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(valuetype, use_spikeins, biotype_to_use)
adata, phasesfilt = RNADataPreparation.qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)

import_dict = Loaders.load_ptm_and_stability(adata)
wp_ensg, ccd_comp, ccdtranscript, wp_max_pol = import_dict["wp_ensg"], import_dict["ccd_comp"], import_dict["ccdtranscript"], import_dict["wp_max_pol"]
name_results = utils.save_gene_names_by_category(adata, wp_ensg, ccd_comp, ccdtranscript)
names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_genes_analyzed, names_ccd_regev_filtered, names_genes_analyzed, names_ccd_filtered = name_results

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in MetaMorpheus results
genemodsBulk = PTMAnalysis.analyze_ptms("input/raw/U2OSBulkAllProteinGroups.tsv", names_ccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_genes_analyzed)
genemodsPhospho = PTMAnalysis.analyze_ptms("input/raw/U2OSPhosphoAllProteinGroups.tsv", names_ccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_genes_analyzed)

# Analyze the modification occupancy for each PTM site
dfBulk, occdfBulk = PTMAnalysis.process_genemods(genemodsBulk)
dfPhospho, occdfPhospho = PTMAnalysis.process_genemods(genemodsPhospho)

# Investigate whether PTM occupancies have correlations with cell division pseudotime (they don't with this dataset)
PTMAnalysis.temporal_ptm_regulation_not_observed(dfBulk, occdfBulk, "Bulk", wp_max_pol, wp_ensg, names_ccdprotein)
PTMAnalysis.temporal_ptm_regulation_not_observed(dfPhospho, occdfPhospho, "Phospho", wp_max_pol, wp_ensg, names_ccdprotein)

# Compare the PTM regulation using PTM site occupancies between transcript and non-transcript regulated groups
PTMAnalysis.compare_ptm_regulation(dfBulk, occdfBulk, "Bulk", names_genes_analyzed, names_ccd_regev_filtered, names_ccdtranscript, names_nonccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_ccdprotein, names_nonccdprotein)
PTMAnalysis.compare_ptm_regulation(dfPhospho, occdfPhospho, "Phospho", names_genes_analyzed, names_ccd_regev_filtered, names_ccdtranscript, names_nonccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_ccdprotein, names_nonccdprotein)

#%% Perform melting point analysis
# Idea: A lower melting point for a protein indicates a higher propensity to unfold, 
#    and so the melting points measured by MS proteomics serve as a useful way to study the protein stability of CCD proteins
# Output: Boxplots illustrating differences in stability for different classes of CCD proteins
ProteinStabilityAnalysis.melting_point_analysis(names_ccdtranscript, names_nonccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_nonccdprotein)

#%% Analyze properties of the different groups
import gzip
import iupred2a.iupred2a_lib
from itertools import groupby
from scipy import stats
import re
proteome = {}
header = ""
with gzip.open("input/raw/uniprot-proteome_UP000005640.fasta.gz", mode="rt") as file_handler:
    for line in file_handler:
        if line.startswith(">"):
            header = line
            proteome[header] = ""
        elif header:
            proteome[header] += line.strip()

gene_name_rgx = re.compile("GN=([^ ]+)")
protein_name_rgx = re.compile("\|([A-Z0-9_]+)")
protein_names = []
gene_names = []
for header, seq in proteome.items():
    gn = gene_name_rgx.findall(header)
    if len(gn) == 0: 
        gene_names.append("")
    else: 
        gene_names.append(gn[0])
    pn = protein_name_rgx.findall(header)
    if len(pn) < 2: 
        protein_names.append("")
    else: 
        protein_names.append(pn[1].split("_")[0])
print(f"{sum(np.array(gene_names) == '')}: number of empty gene names")
print(f"{sum(np.array(protein_names) == '')}: number of empty protein names")
print(f"{sum([gene_names[ii] == protein_names[ii] for ii in np.arange(len(protein_names))])}: number of gene/protein name matches")

hydrophilic = ['D','E','K','R','H','S','C','Y','T','Q','N','U']
hydrophobic = ['G','A','V','L','I','F','W','M','P']
protein_disorder = {} # prot_name, (isDisordered, maxDisorderedSegmentLength, numDisorderedRes, cys, phobic, philic, length)
# protein_names = np.array([protein_name_rgx.findall(header)[1].split('_') for header in proteome.keys()])
# nCysteines = np.array([sum(np.isin(list(seq), ['C','U'])) for seq in proteome.values()])
# nHydrophilic = np.array([sum(np.isin(list(seq), hydrophilic)) for seq in proteome.values()])
# nHydrophobic = np.array([sum(np.isin(list(seq), hydrophobic)) for seq in proteome.values()])
# lengths = np.array([len(seq) for seq in proteome.values()])
totalDisorderedResidues = 0
disordered_proteins = 0
totalResidues = 0

aebersoldNumbers = pd.read_csv("C:/Users/antho/Dropbox/Projects/Nucleoli/AebersoldNumbers.csv", index_col=False)
all_intensities = aebersoldNumbers["Copies Per Cell"]
aebGenes = list(aebersoldNumbers["GeneName"])
nucleoli_intensities = aebersoldNumbers["Copies Per Cell"][np.isin(aebersoldNumbers["GeneName"], np.concatenate(nucleoli))]

for ii, xx in enumerate(proteome.items()):
    header, seq = xx
    if ii % 1000 == 0: print(f"processing protein {ii}")
    if len(seq) == 0: continue
    pn = protein_name_rgx.findall(header)
    if len(pn) < 2: continue
    pn = pn[1].split('_')[0]
    iupred_scores = iupred2a_lib.iupred(seq, "long")
    maxDisorderedLen = max([sum([1 for _ in y]) if x == 1 else 0 for x, y in groupby([1 if x>=0.5 else 0 for x in iupred_scores])])
    isDisordered = maxDisorderedLen > 30    
    if isDisordered: disordered_proteins += int(isDisordered)
    currDisorderedResidues = sum([1 for x in iupred_scores if x>=0.5])
    totalDisorderedResidues += currDisorderedResidues
    totalResidues += len(seq)
    seq_list = list(seq)
    currCysteines = sum(np.isin(seq_list, ['C','U']))
    currHydrophilic = sum(np.isin(seq_list, hydrophilic))
    currHydrophobic = sum(np.isin(seq_list, hydrophobic))
    if pn in protein_disorder: print(f"{pn} already in dictionary")
    protein_disorder[pn] = (isDisordered, maxDisorderedLen, currDisorderedResidues, currCysteines, currHydrophilic, currHydrophobic, len(seq))
print("Fraction of disordered residues: {:.2f}".format(totalDisorderedResidues/totalResidues))
print("Fraction of disordered proteins: {:.2f}".format(disordered_proteins/len(proteome)))

names_bioccd = utils.ccd_gene_names(bioccd)
disorderedMitotic = sum([protein_disorder[nn][0] for nn in names_bioccd if nn in protein_disorder])
totalMitotic = sum([nn in protein_disorder for nn in names_bioccd])
print(f"{disorderedMitotic / totalMitotic}: fraction of disordered mitotic proteins")

disorderedCcdProtein = sum([protein_disorder[nn][0] for nn in names_ccdprotein if nn in protein_disorder and not nn in names_bioccd])
totalCcdProtein = sum([nn in protein_disorder for nn in names_ccdprotein if not nn in names_bioccd])
print(f"{disorderedCcdProtein / totalCcdProtein}: fraction of disordered CCD pseuduotime-only proteins")

def analyze(a, b, a_lab, b_lab, lab_lab, showfliers, filename):
    utils.general_boxplot((a, b), (a_lab, b_lab), "", lab_lab, "", showfliers, filename)
    print(stats.kruskal(a, b))

# disorder
fractDisorderTransreg = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_ccdprotein_transcript_regulated if nn in protein_disorder]
fractDisorderNontransreg = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_ccdprotein_nontranscript_regulated if nn in protein_disorder]
analyze(fractDisorderTransreg, fractDisorderNontransreg, "transreg", "nontransreg", "Fraction Disordered", True, "figures/fractDisordered.png")

fractDisorderBioccd = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_bioccd if nn in protein_disorder]
fractDisorderCcdPseudotime = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder and nn not in names_bioccd]
fractDisorderNonCcd = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_nonccdprotein if nn in protein_disorder]
fractDisorderCcd = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder]
fractDisorderAll = [protein_disorder[nn][2] / protein_disorder[nn][6] for nn in protein_disorder.keys()]
analyze(fractDisorderBioccd, fractDisorderCcdPseudotime, "mitotic", "interphase", "Fraction Disordered", True, "figures/fractDisorderedMitotic.png")
analyze(fractDisorderCcd, fractDisorderNonCcd, "ccd", "nonccd", "Fraction Disordered", True, "figures/fractDisorderedCcd.png")
analyze(fractDisorderCcd, fractDisorderAll, "ccd", "all", "Fraction Disordered", True, "figures/fractDisorderedCcdVsAll.png")

# cysteines
fractCysteineTransreg = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_ccdprotein_transcript_regulated if nn in protein_disorder]
fractCysteineNontransreg = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_ccdprotein_nontranscript_regulated if nn in protein_disorder]
analyze(fractCysteineTransreg, fractCysteineNontransreg, "transreg", "nontransreg", "Fraction Cysteine", True, "figures/fractCysteine.png")

fractCysteineBioccd = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_bioccd if nn in protein_disorder]
fractCysteineInterphase = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder and nn not in names_bioccd]
fractCysteineNonCcd = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_nonccdprotein if nn in protein_disorder]
fractCysteineCcd = [protein_disorder[nn][3] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder]
analyze(fractCysteineBioccd, fractCysteineInterphase, "mitotic", "interphase", "Fraction Cysteine", True, "figures/fractDisorderedMitotic.png")
analyze(fractCysteineCcd, fractCysteineNonCcd, "ccd", "nonccd", "Fraction Cysteine", True, "figures/fractCysteineCcd.png")

# hydrophilic
fractHydrophobicTransreg = [protein_disorder[nn][4] / protein_disorder[nn][6] for nn in names_ccdprotein_transcript_regulated if nn in protein_disorder]
fractHydrophobicNontransreg = [protein_disorder[nn][4] / protein_disorder[nn][6] for nn in names_ccdprotein_nontranscript_regulated if nn in protein_disorder]
analyze(fractHydrophobicTransreg, fractHydrophobicNontransreg, "transreg", "nontransreg", "Fraction Hydrophilic", True, "figures/fractCysteine.png")

fractHydrophobicBioccd = [protein_disorder[nn][4] / protein_disorder[nn][6] for nn in names_bioccd if nn in protein_disorder]
fractHydrophobicInterphase = [protein_disorder[nn][4] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder and nn not in names_bioccd]
analyze(fractHydrophobicBioccd, fractHydrophobicInterphase, "mitotic", "interphase", "Fraction Hydrophilic", True, "figures/fractDisorderedMitotic.png")

# hydrophobic
fractHydrophilicTransreg = [protein_disorder[nn][5] / protein_disorder[nn][6] for nn in names_ccdprotein_transcript_regulated if nn in protein_disorder]
fractHydrophilicNontransreg = [protein_disorder[nn][5] / protein_disorder[nn][6] for nn in names_ccdprotein_nontranscript_regulated if nn in protein_disorder]
analyze(fractHydrophilicTransreg, fractHydrophilicNontransreg, "transreg", "nontransreg", "Fraction Hydrophobic", True, "figures/fractCysteine.png")

fractHydrophilicBioccd = [protein_disorder[nn][5] / protein_disorder[nn][6] for nn in names_bioccd if nn in protein_disorder]
fractHydrophilicInterphase = [protein_disorder[nn][5] / protein_disorder[nn][6] for nn in names_ccdprotein if nn in protein_disorder and nn not in names_bioccd]
analyze(fractHydrophilicBioccd, fractHydrophilicInterphase, "mitotic", "interphase", "Fraction Hydrophobic", True, "figures/fractDisorderedMitotic.png")

# abundance
aebGenesSet = set(aebGenes)
abundanceTransreg = [all_intensities[aebGenes.index(nn)] for nn in names_ccdprotein_transcript_regulated if nn in aebGenesSet]
abundanceNontransreg = [all_intensities[aebGenes.index(nn)] for nn in names_ccdprotein_nontranscript_regulated if nn in aebGenesSet]
abundanceBioccd = [all_intensities[aebGenes.index(nn)] for nn in names_bioccd if nn in aebGenesSet]
abundanceCcdPseudotime = [all_intensities[aebGenes.index(nn)] for nn in names_ccdprotein if nn in aebGenesSet]
abundanceNonccd = [all_intensities[aebGenes.index(nn)] for nn in names_nonccdprotein if nn in aebGenesSet]
abundanceCcd = [all_intensities[aebGenes.index(nn)] for nn in names_ccdprotein if nn in aebGenesSet]
abundanceAll = all_intensities
analyze(abundanceTransreg, abundanceNontransreg, "transreg", "nontransreg", "Abundance", False, "figures/fractDisordered.png")
analyze(abundanceBioccd, abundanceCcdPseudotime, "mitotic", "interphase", "Abundance", False, "figures/fractDisorderedMitotic.png")
analyze(abundanceCcd, abundanceNonccd, "ccd", "nonccd", "Abundance", False, "figures/fractDisorderedCcd.png")
analyze(abundanceCcd, abundanceAll, "ccd", "all", "Abundance", False, "figures/fractDisorderedCcdVsAll.png")

# variants
import re, gzip
anncomma = re.compile(b',(?!\()')
transcript_variant = {}
with gzip.open(filepath, 'rb') as file:
    for line in file:
        if line.startswith(b"#"): continue
        infos = line.split(b'\t')[7].split(b';')
        annotation = [info for info in infos if info.startswith(b"ANN=")]
        if len(annotation) != 1: continue
        annotations = anncomma.split(annotation[0].split(b'ANN=')[1])
        for ann in annotations:
            allele, efftypes, putative_impact, geneName, geneId, featureType, featureId, biotype, exonIntronRank, hgvsDna, hgvsProtein, cdnaPosition, cdsPosition, protPos, distToFeature, warnings = ann.split(b'|')
            if putative_impact in [b"MODERATE", b"HIGH"] and featureId.startswith(b'transcript:ENST') and hgvsProtein: # skip splice acceptor/donator variations by requiring hgvsProtein
                transcriptId = featureId.strip(b'transcript:')
                if transcriptId in transcript_variant: transcript_variant[transcriptId].append(ann)
                else: transcript_variant[transcriptId] = [ann]
                
gene_info = pd.read_csv(f"input/processed/transcriptId_geneName.txt", sep="\t", index_col=False)
transcriptId_proteinName = dict((info[1]["Transcript stable ID"], info[1]["Gene name"]) for info in gene_info.iterrows())
proteinName_variant = {}
for item in transcript_variant.items():
    baseProteinName = transcriptId_proteinName[item[0].decode('ascii')]#.split('-')[0]
    if baseProteinName not in proteinName_variant: proteinName_variant[baseProteinName] = item[1]
    else: proteinName_variant[baseProteinName].extend(item[1])
    
fractWithVariant_nonTransReg = sum([nn in proteinName_variant for nn in names_ccdprotein_nontranscript_regulated]) / len(names_ccdprotein_nontranscript_regulated)
fractWithVariant_transReg = sum([nn in proteinName_variant for nn in names_ccdprotein_transcript_regulated]) / len(names_ccdprotein_transcript_regulated)
print(f"{fractWithVariant_nonTransReg}: fraction of nontranscript regulated genes with variant")
print(f"{fractWithVariant_transReg}: fraction of transcript regulated genes with variant")