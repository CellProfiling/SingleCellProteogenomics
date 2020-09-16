# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 23:50:33 2020

@author: antho
"""

from Bio import SeqIO
import sys, re, math, itertools, os, json
import seaborn as sbn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from SingleCellProteogenomics import utils, FucciCellCycle

MM_BLACKLIST = ["oxidation", "deamidation", "ammonia loss", "water loss", "carbamyl", "carbamidomethyl", # artifacts
    "fe[i", "zinc", "cu[i" , # metals
    "sodium", "potassium", "calcium", "magnesium", # counterions
    "tmt6-plex"] # isotopic labels
MQ_BLACKLIST = ["ox","de", "gl","hy","ac"]
fucci = FucciCellCycle.FucciCellCycle() # Object representing FUCCI cell cycle phase durations

class PeptideModifications:
    def __init__(self, modifications, clearedMods, peptideString, modifiedPeptideSeq, protAcc, protGeneName):
        self.modifications = modifications
        self.clearedMods = set(clearedMods)
        self.peptideString = peptideString
        self.modifiedPeptideSeq = modifiedPeptideSeq
        self.proteinAccession = protAcc
        self.proteinGeneName = protGeneName
        self.modCountSet = False
    def isSamePep(self, modPep):
        return self.peptideString == modPep.peptideString and self.proteinAccession == modPep.proteinAccession and self.proteinGeneName == modPep.proteinGeneName
    def setModCountsMQ(self):
        if self.modCountSet: return
        self.modCountSet = True
        self.phCount = sum([x[1] == 'ph' for x in self.modifications])
        self.oxCount = sum([x[1] == 'ox' for x in self.modifications])
        self.deCount = sum([x[1] == 'de' for x in self.modifications])
        self.glCount = sum([x[1] == 'gl' for x in self.modifications])
        self.hyCount = sum([x[1] == 'hy' for x in self.modifications])
        self.acCount = sum([x[1] == 'ac' for x in self.modifications])
    def isSameWithoutPhosphoMQ(self, modPep):
        self.setModCountsMQ()
        modPep.setModCountsMQ()
        sameOtherMods = self.oxCount == modPep.oxCount \
            and self.deCount == modPep.deCount \
            and self.glCount == modPep.glCount \
            and self.hyCount == modPep.hyCount \
            and self.acCount == modPep.acCount
        return self.isSamePep(modPep) and sameOtherMods

def get_seq_without_mods(seq):
    returnseq = []
    isMod = False
    openBracketCt = 0
    for c in seq:
        if c == "[" and not isMod: isMod = True
        elif c == "[" and isMod: openBracketCt += 1
        elif c == "]" and isMod and openBracketCt > 0: openBracketCt -= 0
        if not isMod: returnseq.append(c)
        if c == "]" and isMod and openBracketCt == 0: isMod = False
    return "".join(returnseq)
 
def get_mod_starts(seq, start, protAcc, protGeneName, blacklist, modBeginChar, modEndChar):
    if not start: return []
    modifications = []
    positions = []
    clearedMods, clearedPos = [], []
    modStr = []
    peptideStr = []
    isMod = False
    openBracketCt = 0
    peptideIdx = int(start)
    for idx, c in enumerate(seq):
        if c == modBeginChar and not isMod: 
            isMod = True
        elif c == modBeginChar and isMod: 
            openBracketCt += 1
        elif c == modEndChar and isMod and openBracketCt > 0: 
            openBracketCt -= 0
        elif not isMod:
            if len(peptideStr) > 0: 
                peptideIdx += 1
            peptideStr.append(c)
        elif c == modEndChar and isMod and openBracketCt == 0: 
            isMod = False
            modifications.append((peptideIdx, "".join(modStr)))
            isBlacklisted = any(mod in "".join(modStr).lower() for mod in blacklist)
            if not isBlacklisted:
                clearedMods.append((peptideIdx, "".join(modStr)))
            modStr = []
        elif c != modBeginChar and isMod: 
            modStr.append(c)
    return PeptideModifications(modifications, clearedMods, "".join(peptideStr), seq, protAcc, protGeneName)

def mmToMqSequence(seq, modBeginChar="[", modEndChar="]"):
    modifications = []
    clearedMods, clearedPos = [], []
    modStr = []
    peptideStr = []
    modPeptideStr = []
    isMod = False
    openBracketCt = 0
    for idx, c in enumerate(seq):
        if c == modBeginChar and not isMod: 
            isMod = True
        elif c == modBeginChar and isMod: 
            openBracketCt += 1
        elif c == modEndChar and isMod and openBracketCt > 0: 
            openBracketCt -= 0
        elif not isMod:
            peptideStr.append(c)
            modPeptideStr.append(c)
        elif c == modEndChar and isMod and openBracketCt == 0: 
            isMod = False
            if "oxidation" in "".join(modStr).lower(): modPeptideStr.extend(list("(ox)"))
            elif "deamidation" in "".join(modStr).lower(): modPeptideStr.extend(list("(de)"))
            elif "acetyl" in "".join(modStr).lower(): modPeptideStr.extend(list("(ac)"))
            elif "hydroxylation" in "".join(modStr).lower(): modPeptideStr.extend(list("(hy)"))
            elif "phospho"  in "".join(modStr).lower(): modPeptideStr.extend(list("(ph)"))
            # else: modPeptideStr.extend(modStr) # makes sure peptides not in MQ search space don't match
            modStr = []
        elif c != modBeginChar and isMod: 
            modStr.append(c)
    return "".join(modPeptideStr)

MOD_PEPTIDE_METHOD = "justmod"
MOD_RATIO_METHOD = "ratioratio"
def get_mod_ratios(modificationsPerPeptide, pepSumIntensities, intensitiesNorm, proteinIds, proteinIntensitiesNorm):
    '''
    Gets the occupancy given intensities for each peptide
    
    Methods:
    anymod method: compare the intensities with mod to those without (default)
    justmod method: compare the intensities with mod to those without; only consider modified-unmodified comparisons and nothing more complicated (similar to Olsen 2010)
    
    Methods:
    sum: use raw sum intensities to calculate modification ratios
    ratioratio: use the differences of ratios from protein ratios
    '''
    modRatios = {}
    modRatioPeptideToKey = {}
    peptidesWithModCounts = []
    peptidesWithoutModCounts = []
    for iii, modPep in enumerate(modificationsPerPeptide):
        if iii % 1000 == 0: print(f"{iii} of {len(modificationsPerPeptide)}")
        if not modPep or MOD_PEPTIDE_METHOD == "justmod" and len(modPep.modifications) > 1: continue
        # currProteinIntensities = proteinIntensities[[modPep.proteinAccession in pp.split(';') for pp in proteinIds],:][0]
        for jjj, clearedMod in enumerate(modPep.clearedMods):
            key=(modPep.proteinAccession, modPep.proteinGeneName, clearedMod)
            if key in modRatios: continue
            intensity = pepSumIntensities[iii]
            peptidesWithMod = [False if not modmod else modPep.isSamePep(modmod) and clearedMod in modmod.clearedMods for modmod in modificationsPerPeptide]
            peptidesWithoutMod = [False if not modmod else modPep.isSamePep(modmod) and not clearedMod in modmod.clearedMods for modmod in modificationsPerPeptide]
            # if MOD_RATIO_METHOD == "sum":
            #     intensityWithMod = np.sum(pepSumIntensities[peptidesWithMod])
            #     intensityWithoutMod = np.sum(pepSumIntensities[peptidesWithoutMod])
            #     if sum(peptidesWithoutMod) > 0:
            #         if intensityWithMod + intensityWithoutMod > 0:
            #             modRatios[key] = intensityWithMod / (intensityWithMod + intensityWithoutMod)
            #             modRatioPeptideToKey[modPep.modifiedPeptideSeq] = key
            # elif all(currProteinIntensities == 0):
            #     continue
            # else:
            #     ratiosWithMod = intensitiesNorm[peptidesWithMod][:,:3] / intensitiesNorm[peptidesWithMod][:,3:]
            #     ratiosWithoutMod = intensitiesNorm[peptidesWithoutMod][:,:3] / intensitiesNorm[peptidesWithoutMod][:,3:]
            #     protRatio = currProteinIntensities[:3] / currProteinIntensities[3:]
            #     lProportion = (protRatio - ratiosWithoutMod) / (ratiosWithMod - protRatio)
            #     modRatios[key] = lProportion / (1 + lProportion)
            #     modRatioPeptideToKey[modPep.modifiedPeptideSeq] = key
            peptidesWithModCounts.append(sum(peptidesWithMod))
            peptidesWithoutModCounts.append(sum(peptidesWithoutMod))
            
    # How many peptides had the mod and didn't have the mod
    plt.hist([peptidesWithModCounts, peptidesWithoutModCounts], label=["Peptides With Mod", "Peptides Without Mod"])
    plt.legend()
    plt.xlabel("Peptide Count Per Modification Site")
    plt.ylabel("Count")
    plt.savefig(f"figures/PeptideCountPerModSite{experiment}.png")
    plt.show()
    plt.close()
    
    return modRatios, modRatioPeptideToKey

def getTmtIntensities(ionString):
    ions = ionString.replace('[','').replace(']','').replace(';',', ').split(', ')
    diagnosticIons = np.array([ion.strip('D').split(':') for ion in ions if ion.startswith("D")])
    intensities = []
    for mass in np.arange(126, 132):
        isDiagnosticIon = [d[0].startswith(str(mass)) for d in diagnosticIons]
        if any(isDiagnosticIon): intensities.append(diagnosticIons[isDiagnosticIon][0][1])
        else: intensities.append(0)
    return np.array(intensities)

# class IntensityPTMAnalysis:
#     def __init__(self, ensg_ccdtranscript, ensg_nonccdtranscript, ensg_ccdprotein, ensg_nonccdprotein, ensg_ccdprotein_transcript_regulated, 
#                ensg_ccdprotein_nontranscript_regulated, genes_analyzed, ccd_regev_filtered, ccd_filtered,
#                names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, names_ccdprotein_transcript_regulated,
#                names_ccdprotein_nontranscript_regulated, names_genes_analyzed, names_ccd_regev_filtered, names_ccd_filtered,
#                bioccd):
#         # Gene IDs
#         self.ensg_ccdtranscript = ensg_ccdtranscript
#         self.ensg_nonccdtranscript = ensg_nonccdtranscript
#         self.ensg_ccdprotein = ensg_ccdprotein
#         self.ensg_nonccdprotein = ensg_nonccdprotein
#         self.ensg_ccdprotein_transcript_regulated = ensg_ccdprotein_transcript_regulated
#         self.ensg_ccdprotein_nontranscript_regulated = ensg_ccdprotein_nontranscript_regulated
#         self.genes_analyzed = genes_analyzed
#         self.ccd_regev_filtered = ccd_regev_filtered
#         self.ccd_filtered = ccd_filtered
#         self.bioccd = bioccd
        
#         # Gene names
#         self.names_ccdtranscript = names_ccdtranscript
#         self.names_nonccdtranscript = names_nonccdtranscript
#         self.names_ccdprotein = names_ccdprotein
#         self.names_nonccdprotein = names_nonccdprotein
#         self.names_ccdprotein_transcript_regulated = names_ccdprotein_transcript_regulated
#         self.names_ccdprotein_nontranscript_regulated = names_ccdprotein_nontranscript_regulated
#         self.names_genes_analyzed = names_genes_analyzed
#         self.names_ccd_regev_filtered = names_ccd_regev_filtered
#         self.names_ccd_filtered = names_ccd_filtered
    
#     def get_olsen_stoichiometry(self):
#         '''Looking at the results of the Olsen 2010 paper in context of these CCD results'''
#         olsen=pd.read_csv("input/raw/Olsen2010OccupancyResults.csv")
        
#         # Look at the distribution of scores & PTM scores
#         plt.scatter(olsen["Best Mascot Score"], olsen["PTM Score"])
#         plt.xlabel("Best Mascot Score")
#         plt.ylabel("PTM Score")
#         plt.show()
#         plt.close()
        
#         # Are the proteins part of any particular group?
#         olsen_ccdprotein = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein and ensg not in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_nonccdprotein = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_nonccdprotein for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_transreg = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein_transcript_regulated and ensg not in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_nontransreg = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein_nontranscript_regulated and ensg not in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_ccdprotein_mitotic = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein and ensg in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_transreg_mitotic = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein_transcript_regulated and ensg in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
#         olsen_nontransreg_mitotic = np.array([False if pd.isna(ensg_list) else any([ensg in self.ensg_ccdprotein_nontranscript_regulated and ensg in self.bioccd for ensg in ensg_list.split(';')]) for ensg_list in olsen["ENSEMBL ID"]])
        
#         # is there a difference between asynchronous phosphosite occupancies?
#         hasStoich = ~pd.isna(olsen["Stoichiometry Async Average [%]"])
#         stoich = olsen["Stoichiometry Async Average [%]"]
#         print("Comparing asynchronous cell phosphosite occupancies")
#         utils.general_boxplot([stoich[hasStoich & olsen_transreg], stoich[hasStoich & olsen_nontransreg], stoich[hasStoich]], (
#             "olsen_transreg", "olsen_nontransreg", "olsen"), "", "Phosphosite Stoichiometry", 
#             "", True, "figures/ccdtransreg_stoicholsen.png")
#         print(f"CCD vs all: {stats.kruskal(stoich[hasStoich & olsen_ccdprotein], stoich[hasStoich])}")
#         print(f"non-CCD vs all: {stats.kruskal(stoich[hasStoich & olsen_nonccdprotein], stoich[hasStoich])}")
#         print(f"CCD vs non-CCD: {stats.kruskal(stoich[hasStoich & olsen_ccdprotein], stoich[hasStoich & olsen_nonccdprotein])}")
#         utils.general_boxplot([stoich[hasStoich & olsen_ccdprotein], stoich[hasStoich & olsen_nonccdprotein], stoich[hasStoich]], (
#             "ccd", "nonccd", "all"), "", "Phosphosite Stoichiometry", 
#             "", True, "figures/ccdstoicholsen.png")
#         print(f"CCD-transreg vs all: {stats.kruskal(stoich[hasStoich & olsen_transreg], stoich[hasStoich])}")
#         print(f"CCD-nontransreg vs all: {stats.kruskal(stoich[hasStoich & olsen_nontransreg], stoich[hasStoich])}")
#         print(f"CCD-transreg vs CCD-nontransreg: {stats.kruskal(stoich[hasStoich & olsen_transreg], stoich[hasStoich & olsen_nontransreg])}")
#         print()
        
#         # are there differences in stoich between phases?
#         for twophases in itertools.combinations(["Stoichiometry - G1  [%]", "Stoichiometry G1/S  [%]", "Stoichiometry - G2  [%]", 
#                                      "Stoichiometry - Early S  [%]", "Stoichiometry - Late S  [%]"], 2):
#             both= ~pd.isna(olsen[twophases[0]]) & ~pd.isna(olsen[twophases[1]])
#             diff = olsen[twophases[0]][both] - olsen[twophases[1]][both]
#             print(f"Comparing {' '.join(twophases)}")
#             print(f"normtest: {stats.normaltest(diff)}")
#             utils.general_boxplot([diff[olsen_ccdprotein[both]], diff[olsen_nonccdprotein[both]], diff], ["ccd_diff", "nonccd_diff", "all"],
#                                   "\n".join(twophases), "Phosphosite Stoichiometry", "", True, "figures/StoichG1G2Diff.png")
#             print(f"CCD vs all: {stats.kruskal(diff[olsen_ccdprotein[both]], diff)}")
#             print(f"non-CCD vs all: {stats.kruskal(diff[olsen_nonccdprotein[both]], diff)}")
#             print(f"CCD vs non-CCD: {stats.kruskal(diff[olsen_ccdprotein[both]], diff[olsen_nonccdprotein[both]])}")
#             utils.general_boxplot([diff[olsen_transreg[both]], diff[olsen_nontransreg[both]], diff], ["ccd_transreg_diff", "ccd_nontransreg_diff", "all"],
#                   "\n".join(twophases), "Phosphosite Stoichiometry G1-G2", "", True, "figures/StoichG1G2Diff.png")
#             print(f"CCD-transreg vs CCD-nontransreg: {stats.kruskal(diff[olsen_transreg[both]], diff[olsen_nontransreg[both]])}")
#             print(f"all vs CCD-nontransreg: {stats.kruskal(diff, diff[olsen_nontransreg[both]])}")
#             print(f"CCD-transreg vs all: {stats.kruskal(diff[olsen_transreg[both]], diff)}")
#             print()
        
#         # are the trajectories different
#         for label in ["Stoichiometry Async Average [%]", "Stoichiometry - G1  [%]", "Stoichiometry G1/S  [%]", "Stoichiometry - G2  [%]", 
#                                      "Stoichiometry - Early S  [%]", "Stoichiometry - Late S  [%]"]:
#             hasStoich = ~pd.isna(olsen[label])
#             stoich = olsen[label]
#             plt.scatter(np.arange(sum(hasStoich)) / sum(hasStoich), sorted(stoich[hasStoich]))
#             plt.scatter(np.arange(sum(hasStoich & olsen_transreg)) / sum(hasStoich & olsen_transreg), sorted(stoich[hasStoich & olsen_transreg]))
#             plt.scatter(np.arange(sum(hasStoich & olsen_nontransreg)) / sum(hasStoich & olsen_nontransreg), sorted(stoich[hasStoich & olsen_nontransreg]))
#             plt.title(f"All, transreg, nontransreg: {label}"); plt.xlabel("Phosphosite fraction"); plt.ylabel("Phosphosite stoichiometry"); plt.show(); plt.close()
#             plt.scatter(np.arange(sum(hasStoich)) / sum(hasStoich), sorted(stoich[hasStoich]))
#             plt.scatter(np.arange(sum(hasStoich & olsen_ccdprotein)) / sum(hasStoich & olsen_ccdprotein), sorted(stoich[hasStoich & olsen_ccdprotein]))
#             plt.scatter(np.arange(sum(hasStoich & olsen_nonccdprotein)) / sum(hasStoich & olsen_nonccdprotein), sorted(stoich[hasStoich & olsen_nonccdprotein]))
#             plt.title(f"All, ccd, nonccd: {label}"); plt.xlabel("Phosphosite fraction"); plt.ylabel("Phosphosite stoichiometry"); plt.show(); plt.close()

def analyze_mod_ratios(modRatios, experiment):
    '''Analyze modification ratios for an experiment'''
    # What is the distribution of mod ratios
    plt.hist(modRatios.values())
    plt.xlabel("Modification Occupancy")
    plt.ylabel("Modification Site Count")
    plt.savefig(f"figures/ModRatioHist{experiment}.png"); plt.show(); plt.close()
            
    # CCD vs Non-CCD
    ccdVsNonccdKruskal = stats.kruskal([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_nonccdprotein])
    ccdVsAllKruskal = stats.kruskal([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed])
    nonccdVsAllKruskal = stats.kruskal([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed])
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_nonccdprotein]), 
        ("ccd", "nonccd"), "", f"Occupancy by {experiment} Intensity", 
        f"p={ccdVsNonccdKruskal[1]}", True, f"figures/boxplotMs1{experiment}.png")
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed]), # True]),
        ("ccd", "all"), "", f"Occupancy by {experiment} Intensity", 
        f"p={ccdVsAllKruskal[1]}",  True, f"figures/boxplotMs1{experiment}_all.png")
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_nonccdprotein],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed]), 
        ("nonccd", "all"), "", f"Occupancy by {experiment} Intensity", 
        f"p={nonccdVsAllKruskal[1]}", True, f"figures/boxplotMs1{experiment}.png")
    
    # Transreg CCD vs nonTransreg CCD
    transregVsNontransregKruskal = stats.kruskal(
        [x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_transcript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_nontranscript_regulated])
    transregVsAllKruskal = stats.kruskal(
        [x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_transcript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed])
    nontransregVsAllKruskal = stats.kruskal(
        [x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_nontranscript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed])
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_transcript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_nontranscript_regulated]), 
        ("transreg", "nontransreg"), "", f"Occupancy by {experiment} Intensity", 
        f"p={transregVsNontransregKruskal[1]}", True, f"figures/boxplotMs1{experiment}.png")
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_transcript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed]), 
        ("transreg", "all"), "", f"Occupancy by {experiment} Intensity", 
        f"p={transregVsAllKruskal[1]}", True, f"figures/boxplotMs1{experiment}.png")
    utils.boxplot_with_stripplot(
        ([x[1] for x in modRatios.items() if x[0][1] in names_ccdprotein_nontranscript_regulated],
        [x[1] for x in modRatios.items() if x[0][1] in names_genes_analyzed]), 
        ("nontransreg", "all"), "", f"Occupancy by {experiment} Intensity", 
        f"p={nontransregVsAllKruskal[1]}", True, f"figures/boxplotMs1{experiment}.png")

def analyze_mm_ms1_ratios(filename):
    '''Uses MetaMorpheus MS1 intensities to look at modification occupancy ratios'''
    # get the start values for each peptide and protein combination
    filePep = pd.read_csv(filename, sep="\t", index_col=False)
    targetsPep = filePep[(filePep["Decoy/Contaminant/Target"] == "T") & (filePep["QValue"] <= 0.01)]
    unambigTargetPeps = np.asarray([not "|" in x for x in targetsPep["Protein Accession"]])
    unambigStarts = np.asarray([not "|" in x for x in targetsPep["Start and End Residues In Protein"]])
    pepSeqProt_start = {}
    for idx, pep in targetsPep[unambigTargetPeps & unambigStarts].iterrows():
        key = (pep[list(targetsPep.columns).index("Base Sequence")], pep[list(targetsPep.columns).index("Protein Accession")])
        value = pep[list(targetsPep.columns).index("Start and End Residues In Protein")]
        if key in pepSeqProt_start: pepSeqProt_start[key].add(value)
        else: pepSeqProt_start[key] = set([value])
    
    # get the modifications and their positions corresponding to MS1 quantified peptides
    quantPeptides = pd.read_csv(os.path.join(os.path.dirname(filename), "AllQuantifiedPeptides.tsv.gz"), sep="\t")
    unambigPeptides = np.asarray([not "|" in x for x in quantPeptides["Protein Groups"]])
    quantUnambigPeptides = np.array(quantPeptides[unambigPeptides])
    startEndInProtein = [list(pepSeqProt_start[(pep[1], pep[2])])[0] if (pep[1], pep[2]) in pepSeqProt_start and len(pepSeqProt_start[(pep[1], pep[2])]) == 1 else "" for pep in quantUnambigPeptides]
    startInProtein = [x.split(" to ")[0].strip("[") if x else "" and len(x[0])  for x in startEndInProtein]
    modificationsPerPeptide = np.array([get_mod_starts(pep[0], startInProtein[iii], pep[2], pep[3], MM_BLACKLIST, "[", "]") for iii, pep in enumerate(quantUnambigPeptides)], dtype=object)
    
    # normalize the intensities
    intensities = quantUnambigPeptides[:, ["Intensity" in col for col in quantPeptides.columns]]
    columnSums = np.sum(intensities, axis=0)
    intensitiesNorm = intensities / columnSums
    pepSumIntensities = np.sum(intensitiesNorm, 1)

    # compare sum of intensity for peptides with the mod to those without the mod
    modRatios_mm_ms1, modRatioPeptides_mm_ms1 = get_mod_ratios(modificationsPerPeptide, pepSumIntensities)
    return modRatios_mm_ms1, modRatioPeptides_mm_ms1, targetsPep, unambigTargetPeps, unambigStarts

def analyze_mm_ms1_ratios_bulk():
    '''Uses MetaMorpheus MS1 intensities to look at modification ratios for bulk PTMs'''
    experiment = "MM_MS1_Bulk"
    result = analyze_mm_ms1_ratios("input/raw/v305ProteinResults/Bulk/search/AllPeptides.psmtsv.gz")
    modRatios_mm_ms1, modificationsPerPeptide, targetsPep, unambigTargetPeps, unambigStarts = result
    analyze_mod_ratios(modRatios_mm_ms1, experiment)
    return result

def analyze_mm_ms2_ratios(targetsPep, unambigTargetPeps, unambigStarts):
    '''Use MS2 ion intensities to get modification occupancy ratios based on TMT intensities from MetaMorpheus'''
    # Get MS2 intensity ratios for phospho
    pepCols = list(targetsPep.columns)
    seqCol, accCol, gncol, startEndCol, miiCol = pepCols.index("Essential Sequence"), pepCols.index("Protein Accession"), pepCols.index("Gene Name"), pepCols.index("Start and End Residues In Protein"), pepCols.index("Matched Ion Intensities")
    targetUnambigPeps = np.array(targetsPep[unambigTargetPeps & unambigStarts])
    startInProtein = [x.strip("[").strip("]").split(" to ")[0] if x else "" and len(x[0]) for x in targetUnambigPeps[:,startEndCol]]
    modificationsPerPeptide = np.array([get_mod_starts(pep[seqCol], startInProtein[iii], pep[accCol], str(pep[gncol]).split(', ')[0].strip('primary:'), MM_BLACKLIST, "[", "]") for iii, pep in enumerate(targetUnambigPeps)], dtype=object)
    
    # normalize the intensities
    matchedIonIntensities = np.array([getTmtIntensities(ions) for ions in targetUnambigPeps[:,miiCol]], dtype=int)
    columnSums = np.sum(matchedIonIntensities, axis=0)
    intensitiesNorm = matchedIonIntensities / columnSums
    pepSumIntensities = np.sum(intensitiesNorm, 1)

    # compare sum of intensity for peptides with the mod to those without the mod
    modRatios_mm_ms2, modRatioPeptides_mm_ms2 = get_mod_ratios(modificationsPerPeptide, pepSumIntensities)
    
    # prepare for comparison with MaxQuant
    mqSeqStrings = np.array([mmToMqSequence(pep[seqCol]) for pep in targetUnambigPeps])
    seqToIntensities_mqMs2 = dict((mqSeqStrings[iii], intensitiesNorm[iii]) for iii in np.arange(len(mqSeqStrings)))
    return modRatios_mm_ms2, modRatioPeptides_mm_ms2, seqToIntensities_mqMs2

def get_mm_ratios_phospho():
    '''gets dictionary of (protein, mod, aa#) : MM MS1 intensity ratio'''
    experiment = "MM_MS1_Phospho"
    result_ms1 = analyze_mm_ms1_ratios("input/raw/v305ProteinResults/phospho/search/AllPeptides.psmtsv.gz")
    modRatios_mm_ms1_phospho, modPeptides_mm_ms1_phospho, targetsPep_mm_ms1_phospho, unambigTargetPeps_mm_ms1_phospho, unambigStarts_mm_ms1_phospho = result_ms1
    analyze_mod_ratios(modRatios_mm_ms1_phospho,  experiment)
    
    experiment = "MM_MS2_Phospho"
    result_ms2 = analyze_mm_ms2_ratios(targetsPep_mm_ms1_phospho, unambigTargetPeps_mm_ms1_phospho, unambigStarts_mm_ms1_phospho)
    modRatios_mm_ms2_phospho, modPeptides_mm_ms2_phospho, seqToIntensities_mqMs2 = result_ms2
    analyze_mod_ratios(modRatios_mm_ms2_phospho, experiment)
    
    return (modRatios_mm_ms1_phospho, modPeptides_mm_ms1_phospho, 
            modRatios_mm_ms2_phospho, modPeptides_mm_ms2_phospho, seqToIntensities_mqMs2)

def get_mq_ratios_phospho():
    # gets dictionary of (protein, aa#) : MQ MS1 intensity ratio
    MIN_PIF, MIN_PEP = 0.75, 0.01
    evidenceAllPepRaw = pd.read_csv("input/raw/v305ProteinResults/Phospho/Maxquant_Search_CDKL5_TMT_phosphoproteomic/evidence.txt.gz", sep="\t")
    peptides = pd.read_csv("input/raw/v305ProteinResults/Phospho/Maxquant_Search_CDKL5_TMT_phosphoproteomic/peptides.txt.gz", sep="\t")
    proteins = pd.read_csv("input/raw/v305ProteinResults/Phospho/Maxquant_Search_CDKL5_TMT_phosphoproteomic/proteinGroups.txt.gz", sep="\t")
    isReversed = np.array([xx.startswith("REV_") for xx in evidenceAllPepRaw["Leading Proteins"]])
    evidenceAllPep = np.array(evidenceAllPepRaw[~isReversed & ~pd.isna(evidenceAllPepRaw["Score"]) & ~pd.isna(evidenceAllPepRaw["PIF"]) & (evidenceAllPepRaw["PIF"] > MIN_PIF) & (evidenceAllPepRaw["PEP"] < MIN_PEP)])
    pepCols = list(peptides.columns)
    accGeneName = pd.read_csv("input/raw/v305ProteinResults/uniprot-proteome_UP000005640_accGeneName.txt.gz", sep="\t")
    accToGn = dict([(xx[0], str(xx[2]).split(' ')[0]) for idx, xx in accGeneName.iterrows()])
    peptidesProtStartDict = dict([(xx[pepCols.index("Sequence")], (xx[pepCols.index("Start position")], xx[pepCols.index("Leading razor protein")])) for idx, xx in peptides.iterrows()])
    evidenceProtStart = [peptidesProtStartDict[xx[0]] for xx in evidenceAllPep] # use sequence to get prot and start
    modificationsPerPeptide = np.array([get_mod_starts(xx[3], peptidesProtStartDict[xx[0]][0], peptidesProtStartDict[xx[0]][1], 
               accToGn[peptidesProtStartDict[xx[0]][1]] if peptidesProtStartDict[xx[0]][1] in accToGn else "", 
               MQ_BLACKLIST, "(", ")") for xx in evidenceAllPep])
    
    modsGenenameSite = {}
    for m in modificationsPerPeptide:
        if m.proteinGeneName not in modsGenenameSite: modsGenenameSite[m.proteinGeneName] = set()
        for mod in m.modifications:
            if mod[1] == "ph":
                modsGenenameSite[m.proteinGeneName].add(mod)
    phoscount = [len(x) for x in modsGenenameSite.values()]
    
    phoscount_ccdprotein = [phoscount[idx] for idx in np.arange(len(phoscount)) if list(modsGenenameSite.keys())[idx] in names_ccdprotein and list(modsGenenameSite.keys())[idx] not in names_bioccd]
    phoscount_nonccdprotein = [phoscount[idx] for idx in np.arange(len(phoscount)) if list(modsGenenameSite.keys())[idx] in names_nonccdprotein and list(modsGenenameSite.keys())[idx] not in names_bioccd]
    phoscount_ccdprotein_transcript_regulated = [phoscount[idx] for idx in np.arange(len(phoscount)) if list(modsGenenameSite.keys())[idx] in names_ccdprotein_transcript_regulated and list(modsGenenameSite.keys())[idx] not in names_bioccd]
    phoscount_ccdprotein_nontranscript_regulated = [phoscount[idx] for idx in np.arange(len(phoscount)) if list(modsGenenameSite.keys())[idx] in names_ccdprotein_nontranscript_regulated and list(modsGenenameSite.keys())[idx] not in names_bioccd]
    phoscount_allgenes = [phoscount[idx] for idx in np.arange(len(phoscount)) if list(modsGenenameSite.keys())[idx] in names_genes_analyzed and list(modsGenenameSite.keys())[idx] not in names_bioccd]
    utils.general_boxplot([phoscount_ccdprotein, phoscount_nonccdprotein, phoscount_ccdprotein_transcript_regulated, phoscount_ccdprotein_nontranscript_regulated, phoscount_allgenes],
                          ["ccd", "nonccd", "ccdtransreg", "ccdnontransreg", "all"],
                          "", "Number of phosphosites\nper protein in U2OS", "", False, "figures/phosphositeCounts.png")
    utils.boxplot_with_stripplot([phoscount_ccdprotein_transcript_regulated, phoscount_ccdprotein_nontranscript_regulated],
                      ["ccdtransreg", "ccdnontransreg"],
                      "", "Number of phosphosites in U2OS", "", False, "figures/phosphositeCounts.png")
    
    # with ms1 intensities
    experiment = "MQ_MS1_Phospho"
    intensities = evidenceAllPep[:, ["Intensity" in col and not "not corrected" in col and not "count" in col for col in evidenceAllPepRaw.columns]]
    columnSums = np.sum(intensities, axis=0)
    intensitiesNorm = intensities / columnSums
    pepSumIntensities = np.sum(intensitiesNorm, 1)
    result_mq_ms1 = get_mod_ratios(modificationsPerPeptide, pepSumIntensities)
    modRatios_mq_ms1, modPeptides_mq_ms1 = result_mq_ms1
    analyze_mod_ratios(modRatios_mq_ms1, experiment)
    
    # normalize the ms2 intensities
    experiment = "MQ_MS2_Phospho"
    intensities = evidenceAllPep[:, ["Reporter intensity" in col and not "not corrected" in col and not "count" in col for col in evidenceAllPepRaw.columns]]
    proteinIds = np.array(proteins["Protein IDs"])
    proteinIntensities = np.array(proteins)[:, ["Reporter intensity" in col and not "not corrected" in col and not "count" in col and not "ppTMT_CDK5L" in col for col in proteins.columns]]
    columnSums = np.sum(intensities, axis=0)
    intensitiesNorm = intensities / columnSums
    pepSumIntensities = np.sum(intensitiesNorm, 1)        
    result_mq_ms2 = get_mod_ratios(modificationsPerPeptide, pepSumIntensities, intensitiesNorm, proteinIds, proteinIntensitiesNorm)
    
    # get the ratios based on the MS2 intensities
    columnSumsProt = np.sum(proteinIntensities)
    proteinIntensitiesNorm = proteinIntensities / columnSumsProt
    proteinCtrlIntNorm = np.sum(proteinIntensitiesNorm[:,:3], axis=1)
    proteinConditionIntNorm = np.sum(proteinIntensitiesNorm[:,3:], axis=1)
    proteinNonzeroIntensity = (proteinCtrlIntNorm != 0) & (proteinConditionIntNorm != 0)
    proteinRatios = proteinCtrlIntNorm[proteinNonzeroIntensity] / proteinConditionIntNorm[proteinNonzeroIntensity]
    proteinRatioIds = proteinIds[proteinNonzeroIntensity]
    proteinRatioIdSet = set(np.unique(np.concatenate([p.split(';') for p in proteinRatioIds])))
    peptideCtrlIntNorm = np.sum(intensitiesNorm[:,:3], axis=1)
    peptideConditionIntNorm = np.sum(intensitiesNorm[:,3:], axis=1)
    peptideNonzeroIntensity = (peptideCtrlIntNorm != 0) & (peptideConditionIntNorm != 0)
    peptideRatios = peptideCtrlIntNorm[peptideNonzeroIntensity] / peptideConditionIntNorm[peptideNonzeroIntensity]

    modRatios = {}
    modRatioPeptideToKey = {}
    peptidesWithModCounts = []
    peptidesWithoutModCounts = []
    for modPep in modificationsPerPeptide:
        modPep.setModCountsMQ()
    isMod = np.array([not mod == False for mod in modificationsPerPeptide])
    phCounts = np.array([modPep.phCount for modPep in modificationsPerPeptide])
    otherCounts = np.array([(modPep.oxCount, modPep.oxCount, modPep.deCount, modPep.glCount, modPep.hyCount, modPep.acCount) for modPep in modificationsPerPeptide])
    sameOtherCounts = np.array([np.sum(otherCounts[iii] - otherCounts, axis=1) == 0 for iii in np.arange(len(otherCounts))])
    pepInfo = np.array([modPep.peptideString for modPep in modificationsPerPeptide])
    samePhDict = dict((phct, phCounts == phct) for phct in np.unique(phCounts))
    withoutPh = dict((phct, phCounts == 0) for phct in np.unique(phCounts))
    indices = np.arange(len(isMod))
    indicesIsMod = set(indices[isMod])
    samePepInfo = np.array([set(indices[pepInfo == pepInfo[iii]]) for iii in indices])
    
    # checking that I'm getting the right number of sites
    evidenceCols = list(evidenceAllPepRaw.columns)
    evidenceAllPepRaw["prottag"] = [f"{e[1][evidenceCols.index('Sequence')]}_O{e[1][evidenceCols.index('Oxidation (M)')]}_" + \
                                    f"A{e[1][evidenceCols.index('Acetyl (Protein N-term)')]}_" + 
                                    f"D{e[1][evidenceCols.index('Deamidation (NQ)')]}_G{e[1][evidenceCols.index('Gln->pyro-Glu')]}_" + 
                                    f"H{e[1][evidenceCols.index('Hydroxyproline')]}" \
                                        for e in evidenceAllPepRaw.iterrows()]
    evidenceAllPepPass = evidenceAllPepRaw[~isReversed & ~pd.isna(evidenceAllPepRaw["Score"]) & ~pd.isna(evidenceAllPepRaw["PIF"]) & (evidenceAllPepRaw["PIF"] > MIN_PIF) & (evidenceAllPepRaw["PEP"] < MIN_PEP)]
    numberOfSitesWithAndWithoutPhospho = len([p[1] for p in evidenceAllPepPass.groupby("prottag") if np.array([len(p[1]) > 1]) and 0 in set(p[1]["Phospho (STY)"]) and any(p[1]["Phospho (STY)"] > 0)])
    print(f"{numberOfSitesWithAndWithoutPhospho}: num sites with and without phospho")

    for iii, modPep in enumerate(modificationsPerPeptide):
        if iii % 1000 == 0: print(f"{iii} of {len(modificationsPerPeptide)}")
        if not modPep or modPep.phCount != 1 or modPep.proteinAccession not in proteinRatioIdSet: continue
        sameOthers = indicesIsMod.intersection(samePepInfo[iii]).intersection(indices[sameOtherCounts[iii]])
        peptidesWithmod = sameOthers.intersection(indices[samePhDict[modPep.phCount]])
        peptidesWithoutmod = sameOthers.intersection(indices[withoutPh[modPep.phCount]])
        peptidesWithModCounts.append(len(peptidesWithmod))
        peptidesWithoutModCounts.append(len(peptidesWithoutmod))
        if (len(peptidesWithmod) > 0 and len(peptidesWithoutmod) > 0):
            key=(modPep.proteinAccession, modPep.proteinGeneName, tuple(np.array(modPep.modifications)[[m[1]=='ph' for m in modPep.modifications]][0]))
            phosPepRatio = np.median(peptideRatios[list(peptidesWithmod.intersection(indices[peptideNonzeroIntensity]))])
            unmodPepRatio = np.median(peptideRatios[list(peptidesWithoutmod.intersection(indices[peptideNonzeroIntensity]))]) #sum(peptideCtrlIntNorm[np.array(list(peptidesWithoutmod), dtype=int)]) / sum(peptideConditionIntNorm[np.array(list(peptidesWithoutmod), dtype=int)])
            prots = [modPep.proteinAccession in acc.split(';') for acc in proteinRatioIds]
            protRatio = proteinRatios[prots][0]
            a = (1-unmodPepRatio/protRatio)/(phosPepRatio/protRatio-1)
            modRatios[key] = a / (1+a)
        
    # How many peptides had the mod and didn't have the mod
    plt.hist([peptidesWithModCounts, peptidesWithoutModCounts], label=["Peptides With Mod", "Peptides Without Mod"])
    plt.legend()
    plt.xlabel("Peptide Count Per Modification Site")
    plt.ylabel("Count")
    plt.savefig(f"figures/PeptideCountPerModSite{experiment}.png")
    plt.show()
    plt.close()
    
    modRatios_mq_ms2, modPeptides_mq_ms2 = result_mq_ms2
    analyze_mod_ratios(modRatios_mq_ms2, experiment)
    
    return modRatios_mq_ms1, modPeptides_mq_ms1, modRatios_mq_ms2, modPeptides_mq_ms2, evidenceAllPep

def compare_ms1_to_ms2_quant():
    '''Use the MaxQuant results to compare MS1 to MS2 quantification for each site'''
    # 1. Get the intensities for each TMT channel from the 6-plex experiment and get the ratio between the peptide with and without the site
    # 2. Get the intensities for each MS1 peak and get the ratio for each site with and without the phospho
    # 3. compare the ratios observed for each.
    plt.scatter(modRatios_mq_ms2.values(), modRatios_mq_ms1.values())
    plt.xlabel("Modification site occupancy, MQ MS2 Phospho")
    plt.ylabel("Modification site occupancy, MQ MS1 Phospho")
    plt.savefig("figures/Ms2Ms1Comparison_MQ.png"); plt.show(); plt.close()
    
    mmkeys = set(modRatios_mm_ms1_phospho.keys()).intersection(modRatios_mm_ms2_phospho.keys())
    modRatio_mm_ms1_ph_vals = [modRatios_mm_ms1_phospho[xx] for xx in mmkeys]
    modRatio_mm_ms2_ph_vals = [modRatios_mm_ms2_phospho[xx] for xx in mmkeys]
    plt.scatter(modRatio_mm_ms2_ph_vals, modRatio_mm_ms1_ph_vals)
    plt.xlabel("Modification site occupancy, MM MS2 Phospho")
    plt.ylabel("Modification site occupancy, MM MS1 Phospho")
    plt.savefig("figures/Ms2Ms1Comparison_MM.png"); plt.show(); plt.close()

def compare_psm_to_ms1_quant():
    '''Compare the PSM level quantificaiton to MS1 quantification in MetaMorheus for each site'''
    # 1. read in the MM quantified peptides
    # 2. take the peak intensity ratio for the peptides with the site and those without
    # 3. compare to the ratios observed for the same site on PSM level
    

def compare_mqMs1_to_mmMs1_quant(modPeptides_mm_ms1_phospho, modPeptides_mq_ms1, bfrr):
    '''Compare the MS1 quantification from MaxQuant to the MS1 quantification from MetaMorpheus for each site'''
    # 1. Compare the ratios quantified for MS1 from MM and MQ for the same peptides. Should have high correlation.
    mqSeqStrings = np.array([mmToMqSequence(pep) for pep in modPeptides_mm_ms1_phospho.keys()])
    mqSeqs = np.array([x.strip('_') for x in modPeptides_mq_ms1.keys()])
    phosphoPepKeys = set(mqSeqStrings).intersection(mqSeqs)
    mmMqOrdering = [list(mqSeqStrings).index(x) for x in phosphoPepKeys if x in phosphoPepKeys]

    plt.scatter([modRatios_mm_ms1_phospho[list(modPeptides_mm_ms1_phospho.items())[idx][1]] for idx in mmMqOrdering], 
                [modRatios_mq_ms1[modPeptides_mq_ms1[f"_{key}_"]] for key in phosphoPepKeys])
    plt.xlabel("Modification site occupancy, MM MS1")
    plt.ylabel("Modification site occupancy, MQ MS1")
    plt.savefig("figures/MmMqMs1Comparison.png"); plt.show(); plt.close()

def compare_mqMs2_to_mmMs2_quant(modRatios_mm_ms2, modRatios_mq_ms2):
    '''Compare the MS1 quantification from MaxQuant to the MS1 quantification from MetaMorpheus for each site'''
    # Compare the MS2 intensities for MQ and MM
    mmSumControl = np.array([sum(seqToIntensities_mqMs2[x.strip('_')][:3]) for x in evidenceAllPep[:,3][hasIntensity]])
    mqSumControl = np.sum(intensitiesNorm[hasIntensity][:,:3], axis=1)
    aboveCutoff = (mmSumControl > 0.001) & (mqSumControl > 0.001)
    plt.scatter(mmSumControl[aboveCutoff], mqSumControl[aboveCutoff])
    
    # Compare the MS2 modification ratios
    mqSeqStrings = np.array([mmToMqSequence(pep) for pep in modPeptides_mm_ms2_phospho.keys()])
    mqSeqs = np.array([x.strip('_') for x in modPeptides_mq_ms2.keys()])
    phosphoPepKeys = set(mqSeqStrings).intersection(mqSeqs)
    mmMqOrdering = [list(mqSeqStrings).index(x) for x in phosphoPepKeys if x in phosphoPepKeys]

    plt.scatter([modRatios_mm_ms2_phospho[list(modPeptides_mm_ms2_phospho.items())[idx][1]] for idx in mmMqOrdering], 
                [modRatios_mq_ms2[modPeptides_mq_ms2[f"_{key}_"]] for key in phosphoPepKeys])
    plt.xlabel("Modification site occupancy, MM MS2")
    plt.ylabel("Modification site occupancy, MQ MS2")
    plt.savefig("figures/MmMqMs2Comparison.png"); plt.show(); plt.close()