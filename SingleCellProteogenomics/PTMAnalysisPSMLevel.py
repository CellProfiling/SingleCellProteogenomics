# -*- coding: utf-8 -*-
"""
Analysis of PTM regulation using bulk and phospho-enriched mass spectrometry (MS) proteomic data:
    - PTM site occupancy was used to infer differences in PTM regulation
    - MS proteomic data was analyzed by MetaMorpheus
    - Protein level results with PTM occupancies are read and processed here

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from Bio import SeqIO
import sys, re, math
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
MIN_MOD_PEP = 1
MIN_TOT_PEP = 5
fucci = FucciCellCycle.FucciCellCycle() # Object representing FUCCI cell cycle phase durations

class PTMAnalysis:
    def __init__(self, ccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, genes_analyzed,
                 ccdprotein, ccd_regev_filtered, nonccdtranscript, nonccdprotein):
        self.ccdtranscript = ccdtranscript
        self.ccdprotein_transcript_regulated = ccdprotein_transcript_regulated
        self.ccdprotein_nontranscript_regulated = ccdprotein_nontranscript_regulated
        self.genes_analyzed = genes_analyzed
        self.ccdprotein = ccdprotein
        self.ccd_regev_filtered = ccd_regev_filtered
        self.nonccdtranscript = nonccdtranscript
        self.nonccdprotein = nonccdprotein

    def analyze_ptms(self, filename):
        '''Load MetaMorpheus protein results and store information regarding protein PTMs'''
        print(f"Loading {filename} ...")
        file = pd.read_csv(filename, sep="\t", index_col=False)
        targets = file[(file["Protein Decoy/Contaminant/Target"] == "T") & (file["Protein QValue"] <= 0.01)]
        modifiedProteins = targets[targets["Sequence Coverage with Mods"].str.replace("[","") != targets["Sequence Coverage with Mods"]]
        modifications = [re.findall('\[.*?\]',s) for s in modifiedProteins["Sequence Coverage with Mods"]]
        unique_mods = set([item for sublist in modifications for item in sublist])
        
        accessions = list(targets["Protein Accession"])
        genes = list(targets["Gene"])
        seqmods = list(targets["Sequence Coverage with Mods"])
        seqs = list(targets["Sequence Coverage"])
        modinfos = list(targets["Modification Info List"])
        genemods = {}
        for idx in range(len(genes)):
            accesplit = str(accessions[idx]).split("|")
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
                covered_fraction = float(sum([1 for c in seq[idxx] if c.isupper()])) / float(len(seq[idxx]))
                effective_length = float(len(seq[idxx])) * covered_fraction
                # mods_per_eff_base = math.log10(float(len(modstemp) + 1) / float(effective_length))
                modstempfiltered = []
                aanums, mods, occupancies, modpeptsss, totpeptsss = [],[],[],[],[]
                for oidx, occcc in enumerate(occupancytemp):
                    isblacklisted = any(mod in modstemp[oidx].lower() for mod in blacklist)
                    isonehitwonder = int(modpepts[oidx]) < MIN_MOD_PEP or int(totpepts[oidx]) < MIN_TOT_PEP
                    if not isblacklisted:
                        modstempfiltered.append(modstemp[oidx])
                    if not "(" in occcc and not isblacklisted and not isonehitwonder: 
                        aanums.append(aatemp[oidx])
                        mods.append(modstemp[oidx])
                        occupancies.append(float(occcc))
                        modpeptsss.append(int(modpepts[oidx]))
                        totpeptsss.append(int(totpepts[oidx]))
                mods_per_eff_base = float(len(modstempfiltered)) / float(effective_length)
                if genegene not in genemods: 
                    genemods[genegene] = ([mods_per_eff_base], [(covered_fraction, effective_length)], [len(modstemp)], modstempfiltered, occupancies, mods, 
                        modpeptsss, totpeptsss, aanums, [accesplit[idxx]])
                else:
                    genemods[genegene][0].append(mods_per_eff_base)
                    genemods[genegene][1].append((covered_fraction, effective_length))
                    genemods[genegene][2].append(len(modstempfiltered))
                    genemods[genegene][3].extend(modstempfiltered)
                    genemods[genegene][4].extend(occupancies)
                    genemods[genegene][5].extend(mods)
                    genemods[genegene][6].extend(modpeptsss)
                    genemods[genegene][7].extend(totpeptsss)
                    genemods[genegene][8].extend(aanums)
                    genemods[genegene][9].append(accesplit[idxx])
        print(f"{str(len(targets))} proteins")
        print(f"{str(len(modifiedProteins))} modified proteins ({str(round(float(len(modifiedProteins))/float(len(targets))*100,2))}%)")
        print(f"{str(len([g for g in genemods.keys() if g in self.ccdtranscript]))}: number of all transcript regulated CCD genes of {len(self.ccdtranscript)} detected.")
        print(f"{str(len([g for g in genemods.keys() if g in self.ccdprotein_transcript_regulated]))}: number of transcript regulated CCD genes from Diana's study of {len(self.ccdprotein_transcript_regulated)} detected.")
        print(f"{str(len([g for g in genemods.keys() if g in self.ccdprotein_nontranscript_regulated]))}: number of non-transcript regulated CCD genes from Diana's study of {len(self.ccdprotein_nontranscript_regulated)} detected.")
        print(f"{str(len([g for g in genemods.keys() if g in self.genes_analyzed]))}: number of proteins of {len(self.genes_analyzed)} detected.")
        print(f"{np.median([genemods[g][1][0][0] for g in genemods.keys() if len(genemods[g][5]) >= 0])}: median coverage for all proteins")
        print(f"{np.median([genemods[g][1][0][0] for g in genemods.keys() if len(genemods[g][5]) > 0])}: median coverage for modified proteins")
        return genemods

    def process_genemods(self, genemods):
        '''Convert information packaged from loading modification information into counts and occupancies'''
        unambigenes, modcts, efflength, coveredfract, modrawcts, modrawlist, avgocc, modsss = [],[],[],[],[],[],[],[]
        protmodgene, protmodacc, modmod, occocc, aanumnum, modpeps, totpeps = [],[],[],[],[], [], []
        for gene in genemods.keys():
            unambigenes.append(gene)
            modcts.append(np.median(genemods[gene][0]))
            coveredfract.append(np.median([g[0] for g in genemods[gene][1]]))
            efflength.append(np.median([g[1] for g in genemods[gene][1]]))
            modrawcts.append(np.median(genemods[gene][2]))
            modrawlist.append(", ".join(genemods[gene][3]))
            avgocc.append(np.median(genemods[gene][4]) if any(genemods[gene][4]) else 0)
            modsss.append(", ".join(genemods[gene][5]))
            for modidx, mod in enumerate(genemods[gene][5]):
                protmodgene.append(gene)
                protmodacc.append(genemods[gene][9])
                modmod.append(mod)
                occocc.append(genemods[gene][4][modidx])
                modpeps.append(genemods[gene][6][modidx])
                totpeps.append(genemods[gene][7][modidx])
                aanumnum.append(genemods[gene][8][modidx])
        df = pd.DataFrame({
            "gene" : unambigenes, 
            "modcts" : modcts, 
            "coveredfract" : coveredfract, 
            "efflength" : efflength, 
            "modrawcts" : modrawcts, 
            "modrawlist" : modrawlist, 
            "avgocc" : avgocc,
            "modlist" : modsss})
        occdf = pd.DataFrame({
            "accession" : protmodacc,
            "gene" : protmodgene,
            "modification" : modmod,
            "residue" : aanumnum,
            "occupancy" : occocc,
            "modpeps" : modpeps,
            "totpeps" : totpeps})
        return df, occdf

    def count_mods(self, analyzethese, df, occdf):
        '''Gather modification observations from a certain subset of proteins'''
        modctsdf = df[df["gene"].isin(analyzethese)]
        modcts = modctsdf["modcts"]
        avgocc = modctsdf["avgocc"]
        modctsoccdf = occdf[occdf["gene"].isin(analyzethese)]
        modocc = modctsoccdf["occupancy"]
        return modctsdf, modcts, avgocc, modctsoccdf, modocc

    def temporal_ptm_regulation_not_observed(self, df, occdf, analysisTitle, wp_max_pol, wp_ensg):
        '''Investigate whether there is a correlation between cell division time and modification occupancies'''
        # Conclusion: there's not much going on in terms of correlation of modification occupancy to the peak expression of the protein in this dataset
        # Time of peak expression
        hngcWithGaps = list(utils.geneIdToHngc_withgaps(wp_ensg))
        ccd_modctsdf, ccd_modcts, ccd_avgocc, ccd_modctsoccdf, ccd_modocc = self.count_mods(self.ccdprotein, df, occdf)
        maxpol_per_gene = np.array([wp_max_pol[hngcWithGaps.index(gene)] if gene in hngcWithGaps else -1 for gene in ccd_modctsoccdf["gene"]])
        # Phase of peak expression  
        phase = np.array(["g1" if pol * fucci.TOT_LEN < fucci.G1_LEN else "g1s" if pol * fucci.TOT_LEN >= fucci.G1_LEN and pol * fucci.TOT_LEN < fucci.G1_LEN + fucci.G1_S_TRANS else "g2" for pol in maxpol_per_gene if pol >= 0])
        g1 = ccd_modocc[maxpol_per_gene >= 0][phase == "g1"]
        g1s = ccd_modocc[maxpol_per_gene >= 0][phase == "g1s"]
        g2 = ccd_modocc[maxpol_per_gene >= 0][phase == "g2"]
        
        # Plot time of peak expression (scatterplot) and phase of peak expression (boxplot)
        utils.general_scatter(maxpol_per_gene[maxpol_per_gene >= 0] * fucci.TOT_LEN, ccd_modocc[maxpol_per_gene >= 0], 
             "Cell Division Time, hrs", "PTM Site Occupancy", f"figures/ModOccVsTimeOfPeakExpression{analysisTitle}.png")
        utils.boxplot_with_stripplot((g1, g1s, g2), ("G1", "S-trans", "G2"),
             "", "Occupancy per Modified Residue", "", False, f"figures/ModsOccupancyBoxplotSelected{analysisTitle}.pdf",alpha=0.2, size=7, jitter=0.25, ylim=(-0.01,1.01))

    def compare_ptm_regulation(self, df, occdf, analysisTitle):
        '''Investigate modification occupancies on differently regulated CCD proteins using MetaMorpheus results'''
        # all genes; regev; mRNA CCD (all); mRNA non-CCD (all); mRNA reg. CCD proteins; non-mRNA reg. CCD proteins; CCD proteins; non-CCD proteins
        all_modctsdf, all_modcts, all_avgocc, all_modctsoccdf, all_modocc = self.count_mods(self.genes_analyzed, df, occdf) 
        ccd_rt_modctsdf, ccd_rt_modcts, ccd_rt_avgocc, ccd_rt_modctsoccdf, ccd_rt_modocc = self.count_mods(self.ccd_regev_filtered, df, occdf)
        ccd_at_modctsdf, ccd_at_modcts, ccd_at_avgocc, ccd_at_modctsoccdf, ccd_at_modocc = self.count_mods(self.ccdtranscript, df, occdf)
        ccd_nt_modctsdf, ccd_nt_modcts, ccd_nt_avgocc, ccd_nt_modctsoccdf, ccd_nt_modocc = self.count_mods(self.nonccdtranscript, df, occdf)
        ccd_t_modctsdf, ccd_t_modcts, ccd_t_avgocc, ccd_t_modctsoccdf, ccd_t_modocc = self.count_mods(self.ccdprotein_transcript_regulated, df, occdf)
        ccd_n_modctsdf, ccd_n_modcts, ccd_n_avgocc, ccd_n_modctsoccdf, ccd_n_modocc = self.count_mods(self.ccdprotein_nontranscript_regulated, df, occdf)
        ccd_modctsdf, ccd_modcts, ccd_avgocc, ccd_modctsoccdf, ccd_modocc = self.count_mods(self.ccdprotein, df, occdf)
        nonccd_modctsdf, nonccd_modcts, nonccd_avgocc, nonccd_modctsoccdf, nonccd_modocc = self.count_mods(self.nonccdprotein, df, occdf)
        
        all_modctsoccdf.to_csv(f"output/all_modctsoccdf_{analysisTitle}.csv", index=False)
        ccd_modctsoccdf.to_csv(f"output/ccd_modctsoccdf_{analysisTitle}.csv", index=False)
        nonccd_modctsoccdf.to_csv(f"output/nonccd_modctsoccdf_{analysisTitle}.csv", index=False)
        ccd_at_modctsoccdf.to_csv(f"output/ccd_at_modctsoccdf_{analysisTitle}.csv", index=False)
        ccd_nt_modctsoccdf.to_csv(f"output/ccd_nt_modctsoccdf_{analysisTitle}.csv", index=False)
        ccd_t_modctsoccdf.to_csv(f"output/ccd_t_modctsoccdf_{analysisTitle}.csv", index=False)
        ccd_n_modctsoccdf.to_csv(f"output/ccd_n_modctsoccdf_{analysisTitle}.csv", index=False)
        
        categories = {
            "all genes detected" : self.count_mods(self.genes_analyzed, df, occdf), 
            # "Regev CCD" : self.count_mods(ccd_regev_filtered, df, occdf), 
            "all transcript CCD genes" : self.count_mods(self.ccdtranscript, df, occdf), 
            "all transcript non-CCD genes" : self.count_mods(self.nonccdtranscript, df, occdf), 
            "transcript regulated CCD proteins" : self.count_mods(self.ccdprotein_transcript_regulated, df, occdf),
            "non-transcript regulated CCD proteins" : self.count_mods(self.ccdprotein_nontranscript_regulated, df, occdf),
            "non-CCD proteins" : self.count_mods(self.nonccdprotein, df, occdf),
            "CCD proteins" : self.count_mods(self.ccdprotein, df, occdf)}
        
        counts_of_modified_proteins = pd.DataFrame({
            "category" : list(categories.keys()),
            # Modification counts
            "modifiedProteins" : [len(v[1][v[1] > 0]) for v in categories.values()],
            "modifiedProtein%" : [float(len(v[1][v[1] > 0])) / float(len(v[1])) * 100 for v in categories.values()],
            "meanModSites" : [np.mean(v[0]["modrawcts"][v[1] > 0]) for v in categories.values()],    
            "medianModSites" : [np.median(v[0]["modrawcts"][v[1] > 0]) for v in categories.values()],    
            "totalModSites" : [sum(v[0]["modrawcts"][v[1] > 0]) for v in categories.values()],
            "meanModPerEffBase, mean mods / seq length" : [np.mean(v[1][v[1] > 0]) for v in categories.values()],
            "medianModPerEffBase, mean mods / seq length" : [np.median(v[1][v[1] > 0]) for v in categories.values()],
            "twoSidedKruskal_vs_allProtsModPerEffBase" : [stats.kruskal(v[1][v[1] > 0], categories[list(categories.keys())[0]][1][categories[list(categories.keys())[0]][1] > 0])[1] for v in categories.values()],
            # Modification occupancies
            "modifiedOccProteins" : [len(set(v[3]['gene'])) for v in categories.values()],
            "modifiedOccProtein%" : [float(len(set(v[3]['gene']))) / float(len(v[1])) * 100 for v in categories.values()],
            "meanOccModSites" : [np.mean([sum(v[3]["gene"] == gene) for gene in set(v[3]['gene'])]) for v in categories.values()],
            "medianOccModSites" : [np.median([sum(v[3]["gene"] == gene) for gene in set(v[3]['gene'])]) for v in categories.values()],
            "totalOccModSites" : [sum([sum(v[3]["gene"] == gene) for gene in set(v[3]['gene'])]) for v in categories.values()],
            "twoSidedKruskal_vs_allProtsOccupancy" : [stats.kruskal(v[4], categories[list(categories.keys())[0]][4])[1] for v in categories.values()],
            # Some figures of merit for the experiment
            "medianCoverage" : [np.median(v[0]["coveredfract"]) for v in categories.values()],
            "medianCoverageOfModProteins" : [np.median(v[0]["coveredfract"][v[1] > 0]) for v in categories.values()]})   
        
        counts_of_modified_proteins.to_csv(f"output/counts_of_modified_proteins_{analysisTitle}.csv", index=False)
        
        fom = "OCCUPANCY PER MODIFIED BASE\n"
        fom += "Note: I excluded artifact modifications and required there to be\n"
        # fom += f"at least {MIN_MOD_PEP} peptides with the modification to count it and\n"
        # fom += f"at least {MIN_TOT_PEP} peptides covering the base total.\n"
        fom += "\n"
        fom += f"mean occupancy per modified residue for all proteins: {np.mean(all_modocc)}\n"
        fom += f"mean occupancy per modified residue for all transcriptionally regulated CCD: {np.mean(ccd_at_modocc)}\n"
        fom += f"mean occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.mean(ccd_t_modocc)}\n"
        fom += f"mean occupancy per modified residue for Diana's non-transcriptionally regulated CCD: {np.mean(ccd_n_modocc)}\n"
        t, p = stats.ttest_ind(ccd_t_modocc, ccd_n_modocc)
        fom += f"two-sided t-test for same means of transcript and non-transcriptionally regulated: {p}\n"
        t, p = stats.ttest_ind(all_modocc, ccd_t_modocc)
        fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
        t, p = stats.ttest_ind(all_modocc, ccd_t_modocc)
        fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
        t, p = stats.ttest_ind(all_modocc, ccd_t_modocc)
        fom += f"two-sided t-test for same means of all and non-transcriptionally regulated: {p}\n"
        fom += "\n"
        fom += f"median occupancy per modified residue for all proteins: {np.median(all_modocc)}\n"
        fom += f"median occupancy per modified residue for all transcriptionally regulated CCD: {np.median(ccd_at_modocc)}\n"
        fom += f"median occupancy per modified residue for Diana's transcriptionally regulated CCD: {np.median(ccd_t_modocc)}\n"
        fom += f"median occupancy per modified residue for Diana's non-transcriptionally regulated: {np.median(ccd_n_modocc)}\n"
        t, p = stats.kruskal(ccd_at_modocc, ccd_n_modocc)
        fom += f"one-sided kruskal for same medians of all transcriptionally reg. and non-transcript regulated: {2*p}\n"
        t, p = stats.kruskal(ccd_t_modocc, ccd_n_modocc)
        fom += f"one-sided kruskal for same medians of Diana's transcriptionally reg. CCD and non-transcript regulated: {2*p}\n"
        t, p = stats.kruskal(np.concatenate((ccd_t_modocc, ccd_n_modocc)), nonccd_modocc)
        fom += f"one-sided kruskal for same medians of Diana's CCD and non-CCD regulated: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_t_modocc)
        fom += f"one-sided kruskal for same medians of all genes and Diana's transcriptionally reg CCD: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_at_modocc)
        fom += f"one-sided kruskal for same medians of all genes and all transcriptionally reg: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_rt_modocc)
        fom += f"one-sided kruskal for same medians of all genes and regev transcriptionally reg: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_n_modocc)
        fom += f"one-sided kruskal for same medians of all genes and non-transcriptionally reg: {2*p}\n"
        t, p = stats.kruskal(all_modocc, nonccd_modocc)
        fom += f"one-sided kruskal for same medians of all genes and non-CCD: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_modocc)
        fom += f"one-sided kruskal for same medians of all genes and CCD: {2*p}\n"
        t, p = stats.kruskal(ccd_modocc, nonccd_modocc)
        fom += f"one-sided kruskal for same medians of CCD genes and non-CCD: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_t_modocc, ccd_n_modocc)
        fom += f"one-sided kruskal for same medians of all, Diana's transcript CCD, and non-transcriptionally regulated: {2*p}\n"
        t, p = stats.kruskal(all_modocc, ccd_at_modocc, ccd_n_modocc)
        fom += f"one-sided kruskal for same medians of all, all transcript CCD, and non-transcriptionally regulated: {2*p}\n"
        fom += "\n"
        
        # Generate boxplots for mod occupancies
        utils.general_boxplot((all_modocc, ccd_t_modocc, ccd_rt_modocc, ccd_at_modocc, ccd_n_modocc, nonccd_modocc), 
            ("All", "Transcript\nReg. CCD", "Regev\nTranscript\nCCD", "All\nTranscript\nCCD", "Non-\nTranscript\nReg. CCD", "Non-\nCCD"),
            "Protein Set", "Occupancy per Modified Residue", "", True, f"figures/ModsOccupancyBoxplot{analysisTitle}.pdf", ylim=(-0.01,1.01))
        # Generate a the same boxplot with a stripplot
        utils.boxplot_with_stripplot((all_modocc, ccd_t_modocc, ccd_rt_modocc, ccd_at_modocc, ccd_n_modocc, nonccd_modocc), 
            ("All", "Transcript\nReg. CCD", "Regev\nTranscript\nCCD", "All\nTranscript\nCCD", "Non-\nTranscript\nReg. CCD", "Non-\nCCD"),
            "Protein Set", "Occupancy per Modified Residue", "", True, f"figures/ModsOccupancyBoxplotSelected{analysisTitle}.pdf",alpha=0.3, size=5, jitter=0.25, ylim=(-0.01,1.01))
        # Generate boxplot for the significant difference higlighted in text
        utils.general_boxplot((all_modocc, ccd_modocc, nonccd_modocc), ("All", "CCD", "Non-\nCCD"),
            "", "Occupancy per Modified Residue", "", False, f"figures/ModsOccupancyBoxplotSelected3{analysisTitle}.pdf", ylim=(-0.01,1.01))
        utils.general_boxplot((ccd_modocc, all_modocc, ccd_at_modocc), ("CCD\nProtein PTMs", "All\nProtein PTMs", "CCD\nTranscript\nProtein PTMs"),
            "", "Occupancy per Modified Residue", "", False, f"figures/ModsOccupancyBoxplotSelected4{analysisTitle}.pdf", ylim=(-0.01,1.01))
        
        # save and print results
        print(fom)
        with open(f"output/modificationCountsOccResults{analysisTitle}.txt", 'w') as file:
            file.write(fom)
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
        }).to_pickle(f"output/modocc_{analysisTitle}.pkl")
   
        