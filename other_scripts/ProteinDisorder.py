# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 20:58:29 2020

@author: antho
"""

import gzip
from iupred2a import iupred2a_lib
from itertools import groupby
from scipy import stats
import re
proteome = {}
header = ""
with gzip.open("iupred2a/data/uniprot-proteome_UP000005640.fasta.gz", mode="rt") as file_handler:
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
hydrophobic = ['A','V','L','I','F','W','M']  # 'G' and 'P' weren't included in this set in the meltome atlas paper 
polar = ['S','T','Y','N','Q']
protein_disorder = {} # prot_name, (isDisordered, maxDisorderedLen, currDisorderedResidues, currCysteines, currHydrophilic, currHydrophobic, len(seq))
totalDisorderedResidues = 0
disordered_proteins = 0
totalResidues = 0

aebersoldNumbers = pd.read_csv("C:/Users/antho/Dropbox/Projects/Nucleoli/AebersoldNumbers.csv", index_col=False)
all_intensities = aebersoldNumbers["Copies Per Cell"]
aebGenes = list(aebersoldNumbers["GeneName"])
# nucleoli_intensities = aebersoldNumbers["Copies Per Cell"][np.isin(aebersoldNumbers["GeneName"], np.concatenate(nucleoli))]

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
    currPolar = sum(np.isin(seq_list, polar))
    if pn in protein_disorder: print(f"{pn} already in dictionary")
    protein_disorder[pn] = (isDisordered, maxDisorderedLen, currDisorderedResidues, currCysteines, currHydrophilic, currHydrophobic, currPolar, len(seq))

print("Fraction of disordered residues: {:.2f}".format(totalDisorderedResidues/totalResidues))
print("Fraction of disordered proteins: {:.2f}".format(disordered_proteins/len(proteome)))

protein_names = list(protein_disorder.keys())
protein_properties = list(protein_disorder.values())
pd.DataFrame({
    "ProteinName" : protein_names,
    "IsDisordered" : [x[0] for x in protein_properties],
    "MaxDisorderedLen" : [x[1] for x in protein_properties],
    "DisorderedResidueCount" : [x[2] for x in protein_properties],
    "CysteineCount" : [x[3] for x in protein_properties],
    "HydrophilicResidueCount" : [x[4] for x in protein_properties],
    "HydrophobicResidueCount" : [x[5] for x in protein_properties],
    "PolarResidueCount" : [x[6] for x in protein_properties],
    "Length" : [x[7] for x in protein_properties]
    }).to_csv("input/processed/ProteinDisorderProperties.csv.gz")