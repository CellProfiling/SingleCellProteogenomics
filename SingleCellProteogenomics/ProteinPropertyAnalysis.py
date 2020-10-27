# -*- coding: utf-8 -*-
"""
Analysis of mass spectrometry (MS) proteomic profiling of protein melting points:
    - Allows evaluation of protein stability, since lower melting points indicate higher propensity for unfolding
    - The results from ten different human cell line and primary samples were averaged for each protein

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

import re, gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from SingleCellProteogenomics import utils

class ProteinProperties:
    '''Analyzes protein properties for evaluating differences in melting temperatures'''
    def __init__(self, wp_ensg, ensg_ccdprotein, 
            ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated, 
            bioccd, ensg_nonccdprotein, ensg_ccdtranscript,
            names_bioccd, names_ccdprotein, 
            names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, 
            names_nonccdprotein, names_ccdtranscript):
        '''Reads in information about protein properties for evaluating differences in melting temps'''
        self.proteinDisorder = dict([(r[1][1], r[1][2:]) for r in pd.read_csv("input/ProteinProperties/ProteinDisorderProperties.csv.gz").iterrows()])
        self.wp_ensg = wp_ensg
        self.ensg_ccdprotein_transcript_regulated = ensg_ccdprotein_transcript_regulated
        self.ensg_ccdprotein_nontranscript_regulated = ensg_ccdprotein_nontranscript_regulated
        self.ensg_ccdprotein = ensg_ccdprotein
        self.ensg_nonccdprotein = ensg_nonccdprotein
        self.ensg_bioccd = bioccd
        self.ensg_ccdtranscript = ensg_ccdtranscript
        self.names_bioccd = names_bioccd
        self.names_ccdprotein = names_ccdprotein
        self.names_ccdprotein_transcript_regulated = names_ccdprotein_transcript_regulated
        self.names_ccdprotein_nontranscript_regulated = names_ccdprotein_nontranscript_regulated
        self.names_nonccdprotein = names_nonccdprotein
        self.names_ccdtranscript = names_ccdtranscript
        self.names_hpamapped = pd.read_csv("input/ProteinProperties/proteinatlas.tsv.gz", sep="\t")["Gene"]
        self.test_results = {}

    def analyze2(self, aa, bb, aa_lab, bb_lab, property_label, testNonzeroDistribution=False):
        if testNonzeroDistribution:
            self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}_NonzeroDistribution")] = ("Median", np.median(aa[aa!=0]), np.median(bb[bb!=0]), 
                    "Kruskal", stats.kruskal(aa[aa!=0], bb[bb!=0])[1] * 2) # one-sided
            self.test_results[(property_label, f"{aa_lab}_eqvar_{bb_lab}_NonzeroDistribution")] = ("Std", np.std(aa[aa!=0]), np.std(bb[bb!=0]), 
                   "EqualVarianceLevene", stats.levene(aa[aa!=0], bb[bb!=0])[1])
        else:
            self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}")] = ("Median", np.median(aa), np.median(bb), 
                    "Kruskal", stats.kruskal(aa, bb)[1] * 2) # one-sided
            self.test_results[(property_label, f"{aa_lab}_eqvar_{bb_lab}")] = ("Std", np.std(aa), np.std(bb), 
                   "EqualVarianceLevene", stats.levene(aa, bb)[1])

    def analyze3(self, aa, bb, cc, aa_lab, bb_lab, cc_lab, property_label, testNonzeroDistribution=False):#, lab_lab, showfliers, fileprefix):
        '''General boxplot and kruskal-wallis test'''
        self.analyze2(aa, bb, aa_lab, bb_lab, property_label, testNonzeroDistribution)
        self.analyze2(bb, cc, bb_lab, cc_lab, property_label, testNonzeroDistribution)
        self.analyze2(aa, cc, aa_lab, cc_lab, property_label, testNonzeroDistribution)
        
    def analyzeAll(self, property_label, *groups, testNonzeroDistribution=False):
        self.test_results[(property_label), "allgroups"] = ("N/A", "N/A", "N/A", 
                        "Kruskal", stats.kruskal(*[group[group != 0] if testNonzeroDistribution else group for group in groups])[1])
    
    def analyze_melting_points(self):
        '''Gather measurements of melting points for each protein across 10 cell lines; get the median per protein and evaluate for differences between CCD groups'''
        # Load melting points for each protein across 10 samples
        meltingDf = pd.read_csv("input/ProteinProperties/human.proteinstability.csv.gz")
        meltingDf = meltingDf[pd.notna(meltingDf["quan_norm_meltPoint"])]
        print(f"{','.join(np.unique(meltingDf['cell_line_or_type']))}: unique human cell samples\n")
        self.all_temps, allnonccdtranscript, allccdtranscript, transcript_reg, nontranscr_reg, nonccd_temps = [],[],[],[],[],[]
        self.all_protnames, ant, att, trp, ntp, nnp = [],[],[],[],[],[]
        med = meltingDf.groupby("gene_name")["quan_norm_meltPoint"].median()
        self.all_temps = np.asarray(med)
        self.all_protnames = np.asarray(med.index)
        self.all_protnames_list = list(self.all_protnames)
        self.all_protnames_set = set(self.all_protnames)
        self.ccdprotein_temps = np.asarray(med[np.isin(self.all_protnames, self.names_ccdprotein)])
        self.nonccdprotein_temps = np.asarray(med[np.isin(self.all_protnames, self.names_nonccdprotein)])
        self.ccdprotein_transreg_temps = np.asarray(med[np.isin(self.all_protnames, self.names_ccdprotein_transcript_regulated)])
        self.ccdprotein_nontransreg_temps = np.asarray(med[np.isin(self.all_protnames, self.names_ccdprotein_nontranscript_regulated)])
        self.ccdprotein_mitotic_temps = np.asarray(med[np.isin(self.all_protnames, self.names_bioccd)])
        self.ccdprotein_interphase_temps = np.asarray(med[np.isin(self.all_protnames, self.names_ccdprotein[~np.isin(self.names_ccdprotein, self.names_bioccd)])])
        self.ccdtranscript_temps = np.asarray(med[np.isin(self.all_protnames, self.names_ccdtranscript)])
        self.allmapped_temps = np.asarray(med[np.isin(self.all_protnames, self.names_hpamapped)])
        self.variable_temps = np.concatenate((self.ccdprotein_temps, self.nonccdprotein_temps))
        
        # Aggregate measurements per protein
        property_label = "Tm"
        self.analyze3(self.ccdprotein_temps, self.nonccdprotein_temps, self.allmapped_temps, "ccd", "nonccd", "allmapped", property_label)
        self.analyze3(self.ccdprotein_transreg_temps, self.ccdprotein_nontransreg_temps, self.allmapped_temps, "transreg", "nontransreg", "allmapped", property_label)
        self.analyze2(self.variable_temps, self.allmapped_temps, "variable", "allmapped", property_label)

    def getSequenceFractions(self, idx, isLength, names):
        '''Gets the feature as fraction of length (or the log2-length if returning length)'''
        protLen_idx = 7
        sequenceFractions = np.array([np.log2(self.proteinDisorder[nn][idx]) if isLength else self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in names if nn in self.proteinDisorder])
        return sequenceFractions

    def analyzeSequenceFraction(self, idx, property_label, isLength=False, testNonzeroDistribution=False):
        '''Performs statistical tests on fraction of sequence with feature (or just length of sequence)'''
        fractPropertyAll = self.getSequenceFractions(idx, isLength, self.proteinDisorder.keys())
        fractPropertyAllMapped = self.getSequenceFractions(idx, isLength, self.names_hpamapped)
        fractPropertyTransreg = self.getSequenceFractions(idx, isLength, self.names_ccdprotein_transcript_regulated)
        fractPropertyNontransreg = self.getSequenceFractions(idx, isLength, self.names_ccdprotein_nontranscript_regulated)
        fractPropertyBioccd = self.getSequenceFractions(idx, isLength,  self.names_bioccd)
        fractPropertyCcdPseudotime = self.getSequenceFractions(idx, isLength, self.names_ccdprotein[~np.isin(self.names_ccdprotein, self.names_bioccd)])
        fractPropertyNonCcd = self.getSequenceFractions(idx, isLength,  self.names_nonccdprotein)
        fractPropertyCcd = self.getSequenceFractions(idx, isLength,  self.names_ccdprotein)
        fractPropertyTranscriptCCD = self.getSequenceFractions(idx, isLength,  self.names_ccdtranscript)
        fractPropertyVariable = np.concatenate((fractPropertyCcd, fractPropertyNonCcd))
        self.analyze3(fractPropertyTransreg, fractPropertyNontransreg, fractPropertyAllMapped, 
                      "transreg", "nontransreg", "allmapped", property_label, 
                      testNonzeroDistribution=testNonzeroDistribution)
        self.analyze3(fractPropertyCcd, fractPropertyNonCcd, fractPropertyAllMapped, 
                      "ccd", "nonccd", "allmapped", property_label, 
                      testNonzeroDistribution=testNonzeroDistribution)
        self.analyze2(fractPropertyVariable, fractPropertyAllMapped, "variable", "allmapped", property_label, testNonzeroDistribution=testNonzeroDistribution)
        if property_label == "Disorder":
            self.analyze3(fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyAllMapped, 
                          "mitoticCCD", "interphaseCCD", "allmapped", property_label, 
                          testNonzeroDistribution=testNonzeroDistribution)
            self.analyze2(fractPropertyCcdPseudotime, fractPropertyNonCcd, 
                          "interphaseCCD", "nonccd", property_label, testNonzeroDistribution=testNonzeroDistribution)
            self.analyze2(fractPropertyBioccd, fractPropertyNonCcd, 
                          "mitoticCCD", "nonccd", property_label, testNonzeroDistribution=testNonzeroDistribution)
        return fractPropertyAll, fractPropertyAllMapped, fractPropertyVariable, fractPropertyTransreg, fractPropertyNontransreg, fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyNonCcd, fractPropertyCcd, fractPropertyTranscriptCCD

    def analyze_disorder(self):
        disorderedMitotic = sum([self.proteinDisorder[nn][0] for nn in self.names_bioccd if nn in self.proteinDisorder])
        totalMitotic = sum([nn in self.proteinDisorder for nn in self.names_bioccd])
        print(f"{disorderedMitotic / totalMitotic}: fraction of disordered mitotic proteins")
        disorderedCcdProtein = sum([self.proteinDisorder[nn][0] for nn in self.names_ccdprotein if nn in self.proteinDisorder and not nn in self.names_bioccd])
        totalCcdProtein = sum([nn in self.proteinDisorder for nn in self.names_ccdprotein if not nn in self.names_bioccd])
        print(f"{disorderedCcdProtein / totalCcdProtein}: fraction of disordered CCD pseuduotime-only proteins")
        results = self.analyzeSequenceFraction(2, "Disorder", testNonzeroDistribution=True)
        self.fractDisorderAll, self.fractDisorderAllMapped, self.fractDisorderVariable, self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderBioccd, self.fractDisorderCcdPseudotime, self.fractDisorderNonCcd, self.fractDisorderCcd, self.fractDisorderTranscriptCCD = results
        
    def analyze_cysteines(self):
        results = self.analyzeSequenceFraction(3, "Cysteines", testNonzeroDistribution=True)
        self.fractCysteinesAll, self.fractCysteinesAllMapped, self.fractCysteinesVariable, self.fractCysteinesTransreg, self.fractCysteinesNontransreg, self.fractCysteinesBioccd, self.fractCysteinesCcdPseudotime, self.fractCysteinesNonCcd, self.fractCysteinesCcd, self.fractCysteinesTranscriptCCD = results

    def analyze_hydrophilic(self):
        results = self.analyzeSequenceFraction(4, "Hydrophilic")
        self.fractHydrophilicAll, self.fractHydrophilicAllMapped, self.fractHydrophilicVariable, self.fractHydrophilicTransreg, self.fractHydrophilicNontransreg, self.fractHydrophilicBioccd, self.fractHydrophilicCcdPseudotime, self.fractHydrophilicNonCcd, self.fractHydrophilicCcd, self.fractHydrophilicTranscriptCCD = results

    def analyze_hydrophobic(self):
        results = self.analyzeSequenceFraction(5, "Hydrophobic")
        self.fractHydrophobicAll, self.fractHydrophobicAllMapped, self.fractHydrophobicVariable, self.fractHydrophobicTransreg, self.fractHydrophobicNontransreg, self.fractHydrophobicBioccd, self.fractHydrophobicCcdPseudotime, self.fractHydrophobicNonCcd, self.fractHydrophobicCcd, self.fractHydrophobicTranscriptCCD = results

    def analyze_polar(self):
        results = self.analyzeSequenceFraction(6, "Polar")
        self.fractPolarAll, self.fractPolarAllMapped, self.fractPolarVariable, self.fractPolarTransreg, self.fractPolarNontransreg, self.fractPolarBioccd, self.fractPolarCcdPseudotime, self.fractPolarNonCcd, self.fractPolarCcd, self.fractPolarTranscriptCCD = results

    def analyze_length(self):
        results = self.analyzeSequenceFraction(7, "Log2 Length", True)
        self.lengthAll, self.lengthAllMapped, self.lengthVariable, self.lengthTransreg, self.lengthNontransreg, self.lengthBioccd, self.lengthCcdPseudotime, self.lengthNonCcd, self.lengthCcd, self.lengthTranscriptCCD = results

    def get_property_and_temps(self, propertyValues, names_filter):
        '''Filters temperatures given a list of names'''
        self.all_protnames_set = set(self.all_protnames)
        self.all_protnames_list = list(self.all_protnames)
        filterNames = [nn in self.all_protnames_set for nn in names_filter if nn in self.proteinDisorder] 
        filteredPropertyValues = np.asarray(propertyValues)[filterNames]
        filteredTempValues = [self.all_temps[self.all_protnames_list.index(nn)] for nn in names_filter if nn in self.all_protnames_set and nn in self.proteinDisorder]
        return filteredPropertyValues, filteredTempValues
        
    def temp_property_scatter(self, transregProperty, nontransregProperty, nonccdProperty, propertyLabel):
        plt.scatter(*self.get_property_and_temps(transregProperty, self.names_ccdprotein_transcript_regulated), label="Trans Reg")
        plt.scatter(*self.get_property_and_temps(nontransregProperty, self.names_ccdprotein_nontranscript_regulated), label="Non-Trans Reg")
        plt.scatter(*self.get_property_and_temps(nonccdProperty, self.names_nonccdprotein), label="Non-CCD")
        plt.xlabel(f"Fraction {propertyLabel} Residues")
        plt.ylabel("Melting Point (°C)")
        plt.legend()
        plt.savefig(f"figures/{propertyLabel}VsTm.png")
        plt.close()
    
    def apply_linregress(self, allProperty, propertyLabel):
        names = self.proteinDisorder.keys()
        propertyValues, propertyTemps = self.get_property_and_temps(allProperty, names)
        linModel = stats.linregress(propertyValues, propertyTemps)
        plt.scatter(propertyValues, propertyTemps, alpha=0.1)
        xfit = np.min(propertyValues) + np.arange(100) / 100 * (np.max(propertyValues) - np.min(propertyValues))
        yfit = linModel.intercept + xfit * linModel.slope
        plt.scatter(xfit, yfit)
        plt.xlabel(f"Fraction {propertyLabel} Residues")
        plt.ylabel("Melting Point (°C)")
        plt.savefig(f"figures/All{propertyLabel}VsTm.png")
        plt.close()
        print(f"{linModel.slope}: slope for all melting temperatures vs fract {propertyLabel} residues")
        print(f"{linModel.rvalue**2}: r-squared for all melting temperatures vs fract {propertyLabel} residues")
        print(f"{linModel.pvalue}: p-value for nonzero slope for all melting temperatures vs fract {propertyLabel} residues")
    
    def tm_scatters(self):
        self.temp_property_scatter(self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderNonCcd, "Disorder")
        self.apply_linregress(self.fractDisorderAll, "Disorder")
        self.apply_linregress(self.lengthAll, "Length")
        
    def generate_properties_table(self):
        '''Make a table with all the values analyzed'''
        protDisorderGeneList = list(self.proteinDisorder.keys())
        genes = np.sort(np.unique(np.concatenate((protDisorderGeneList, self.all_protnames_list))))
        lengths, disordered, hydrophob = [],[],[]
        for nn in genes:
            isGene = nn in self.proteinDisorder
            idx = -1 if not isGene else protDisorderGeneList.index(nn)                
            lengths.append("" if not isGene else self.lengthAll[idx])
            disordered.append("" if not isGene else self.fractDisorderAll[idx])
        pd.DataFrame({
            "Protein" : genes,
            "IsCCD" : np.isin(genes, self.names_ccdprotein),
            "IsNonCCD" : np.isin(genes, self.names_nonccdprotein),
            "IsInterphaseCCD" : np.isin(genes, self.names_ccdprotein) & ~np.isin(genes, self.names_bioccd),
            "IsMitoticCCD" : np.isin(genes, self.names_bioccd),
            "IsTranscriptRegCCD" : np.isin(genes, self.names_ccdprotein_transcript_regulated),
            "IsNonTranscriptRegCCD" : np.isin(genes, self.names_ccdprotein_nontranscript_regulated),
            "Melting Temperature (deg C)" : [self.all_temps[self.all_protnames_list.index(nn)] if nn in self.all_protnames_set else "" for nn in genes],
            "Log2 Length" : lengths,
            "Disordered Residue Fraction" : disordered,
            }).to_csv("output/PropertiesTable.csv", index=False)

    def statistical_properties_table(self):
        '''Multiple testing correction; summarize the statistical results'''
        #self.test_results # (comparisonName, propertyName),  (valuename, valueA, valueB, testName, testPvalue)
        self.pvalues = np.array([x[1][4] for x in self.test_results.items()], dtype=float)
        self.pvalues_bonf, self.pass_bonf = np.zeros(len(self.pvalues)), np.zeros(len(self.pvalues))
        
        # Correct equal variance tests
        testnames = np.array([x[1][3] for x in self.test_results.items()])
        eqvartests = testnames == "EqualVarianceLevene"
        eqvar_bonf, eqvar_bonfpass = utils.bonf(0.05, self.pvalues[eqvartests])
        self.pvalues_bonf[eqvartests] = eqvar_bonf
        self.pass_bonf[eqvartests] = eqvar_bonfpass
        
        # Correct property-wise tests
        properties = np.array([x[0][0] for x in self.test_results.items()])
        for prop in np.unique(properties):
            isprop = properties == prop
            prop_bonf, prop_bonfpass = utils.bonf(0.05, self.pvalues[isprop])
            self.pvalues_bonf[isprop] = prop_bonf
            self.pass_bonf[isprop] = prop_bonfpass
            
        pd.DataFrame({
            "Property" : properties,
            "Comparison" : [x[0][1] for x in self.test_results.items()],
            "ValueName" : [x[1][0] for x in self.test_results.items()],
            "ValueA" : [x[1][1] for x in self.test_results.items()],
            "ValueB" : [x[1][2] for x in self.test_results.items()],
            "TestName" : testnames,
            "PValue" : self.pvalues,
            "PValue_Bonf": self.pvalues_bonf,
            "PValue_BonfAlpha0.05Pass" : self.pass_bonf
            }).to_csv("output/PropertyStatistics.csv", index=False)
    
    def boxplot_summary(self, allmapped, variable, nonccd, ccd, transreg, nontransreg, ylabel, showfliers):
        tempsforboxes = (allmapped, variable, nonccd, ccd, transreg, nontransreg)
        utils.general_boxplot_setup(tempsforboxes, [str(len(x)) for x in tempsforboxes], "", ylabel, "", showfliers)
        
    def boxplot_saveclose(self, filename):
        plt.savefig(filename)
        plt.close()
    
    def generate_statistical_boxplots(self):
        # Temps (ccd vs nonccd; transreg vs nontransreg)
        values = (self.allmapped_temps, self.variable_temps, 
                       self.nonccdprotein_temps, self.ccdprotein_temps, 
                       self.ccdprotein_transreg_temps, self.ccdprotein_nontransreg_temps)
        self.boxplot_summary(*values, "Melting Point (°C)", True)
        utils.barplot_annotate_brackets(2, 3, f"p={'%.1e' % self.test_results['Tm','ccd_vs_nonccd'][4]}", np.arange(len(values)), [max(vv) for vv in values]) # highlight CCD vs nonccd
        utils.barplot_annotate_brackets(4, 5, f"p={'%.1e' % self.test_results['Tm','transreg_vs_nontransreg'][4]}", np.arange(len(values)), [5+max(vv) for vv in values]) # highlight CCD vs nonccd        
        self.boxplot_saveclose("figures/ProteinMeltingPointBox.pdf") # highlight transreg vs nontransreg
        
        # Disorder ()
        values = (self.fractDisorderAllMapped, self.fractDisorderVariable, 
                       self.fractDisorderNonCcd, self.fractDisorderCcd, 
                       self.fractDisorderTransreg, self.fractDisorderNontransreg)
        self.boxplot_summary(*values, "Disorder", False)
        utils.barplot_annotate_brackets(0, 1, f"p={'%.1e' % self.test_results['Disorder','variable_vs_allmapped_NonzeroDistribution'][4]}", np.arange(len(values)), [max(vv) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/DisorderBox.pdf")
        
        # Length ()
        values = (self.lengthAllMapped, self.lengthVariable, 
                       self.lengthNonCcd, self.lengthCcd, 
                       self.lengthTransreg, self.lengthNontransreg)
        self.boxplot_summary(*values, "Log2 Length", False)
        utils.barplot_annotate_brackets(0, 1, f"p={'%.1e' % self.test_results['Log2 Length','variable_vs_allmapped'][4]}", np.arange(len(values)), [np.percentile(vv, 98) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/LengthBox.pdf")
                
    
    def proportion_test(self, phosHuman, names_subset, names_alt=[]):
        '''Test the proportions of kinase families with Fisher exact test'''
        subset = np.isin(phosHuman["SUBSTRATE"], names_subset)
        mapped_minus = np.isin(phosHuman["SUBSTRATE"], self.names_hpamapped if len(names_alt) == 0 else names_alt) & ~subset # mutually exclusive mapped proteome
        counts_subset = phosHuman["KINASE_FAMILY"][subset].value_counts()
        counts_mappedminus = phosHuman["KINASE_FAMILY"][mapped_minus].value_counts()
        labels_subset = counts_subset.index
        labels_mappedminus = counts_mappedminus.index
        proportion_subset = np.array([phosHuman["KINASE_FAMILY"][subset].value_counts()[idx] / sum(phosHuman["KINASE_FAMILY"][subset].value_counts()) for idx in np.arange(len(labels_subset))])
        proportion_mappedminus = np.array([counts_mappedminus[idx] / sum(phosHuman["KINASE_FAMILY"][mapped_minus].value_counts()) for idx in np.arange(len(labels_mappedminus))])
        
        labels_comb = np.unique(np.concatenate((labels_subset, labels_mappedminus)))
        counts_comb_subset = np.array([counts_subset[np.array(labels_subset) == x][0] if x in labels_subset else 0 for x in labels_comb])
        counts_comb_mappedminus = np.array([counts_mappedminus[np.array(labels_mappedminus) == x][0] if x in labels_mappedminus else 0 for x in labels_comb])
        values_comb_subset = np.array([proportion_subset[np.array(labels_subset) == x][0] if x in labels_subset else 0 for x in labels_comb])
        values_comb_mappedminus = np.array([proportion_mappedminus[np.array(labels_mappedminus) == x][0] if x in labels_mappedminus else 0 for x in labels_comb])
        fisher_comb_subset = np.array([stats.fisher_exact([
            [counts_comb_subset[idx], sum(counts_comb_subset) - counts_comb_subset[idx]], 
            [counts_comb_mappedminus[idx], sum(counts_comb_mappedminus) - counts_comb_mappedminus[idx]]], "greater") for idx in np.arange(len(labels_comb))])
        return labels_comb, counts_comb_subset, values_comb_subset, counts_comb_mappedminus, values_comb_mappedminus, fisher_comb_subset
        
    def kinase_families(self):
        '''Investigate whether there are differences in upstream kinases from mapped proteins'''
        print("Running kinase family analysis")
        kinaseFams = pd.read_csv("input/ProteinProperties/KinHubKinaseFamilies.csv")
        kinFamDict = dict([(row[7], row[4]) for idx, row in kinaseFams.iterrows()])
        phosphositeplus = pd.read_csv("input/ProteinProperties/KinaseSubstrateDatasetWithoutHeader.txt.gz", sep="\t")
        phosphositeplus["KINASE_FAMILY"] = [kinFamDict[x] if x in kinFamDict else "Other" for x in phosphositeplus["KIN_ACC_ID"]]
        phosHuman = phosphositeplus[phosphositeplus["SUB_ORGANISM"] == "human"]
        
        # Are there any overrepresented kinase families upstream phosphosites on CCD or non-CCD proteins?
        labels_comb_ccd, counts_comb_ccd, values_comb_ccd, counts_comb_mappedminus_ccd, values_comb_mappedminus_ccd, fisher_comb_ccd = self.proportion_test(phosHuman, self.names_ccdprotein)
        labels_comb_nonccd, counts_comb_nonccd, values_comb_nonccd, counts_comb_mappedminus_nonccd, values_comb_mappedminus_nonccd, fisher_comb_nonccd = self.proportion_test(phosHuman, self.names_nonccdprotein)
        labels_comb_transreg, counts_comb_transreg, values_comb_transreg, counts_comb_mappedminus_transreg, values_comb_mappedminus_transreg, fisher_comb_transreg = self.proportion_test(phosHuman, self.names_ccdprotein_transcript_regulated)
        labels_comb_nontransreg, counts_comb_nontransreg, values_comb_nontransreg, counts_comb_mappedminus_nontransreg, values_comb_mappedminus_nontransreg, fisher_comb_nontransreg = self.proportion_test(phosHuman, self.names_ccdprotein_nontranscript_regulated)
        
        allfisher = np.array(np.concatenate((fisher_comb_ccd[:,1], fisher_comb_nonccd[:,1], 
                                            fisher_comb_transreg[:,1], fisher_comb_nontransreg[:,1])), dtype=float)
        pvals_corrected_bonf, reject_bonf = utils.bonf(0.05, allfisher)
        pd.DataFrame({
            "Group" : np.concatenate((["ccd"] * len(labels_comb_ccd), 
                                      ["nonccd"] * len(labels_comb_nonccd),
                                      ["transregCCD"] * len(labels_comb_nonccd),
                                      ["nontransregCCD"] * len(labels_comb_nontransreg))),
            "KinaseFamily" : np.concatenate((labels_comb_ccd, labels_comb_nonccd, 
                                              labels_comb_transreg, labels_comb_nontransreg)),
            "PhosphositeCount_MappedProteome" : np.concatenate((counts_comb_mappedminus_ccd, counts_comb_mappedminus_nonccd, 
                                                 counts_comb_mappedminus_transreg, counts_comb_mappedminus_nontransreg)),
            "FractionPhosphositesDownstreamOfKinase_MappedProteome" : np.concatenate((values_comb_mappedminus_ccd, values_comb_mappedminus_nonccd, 
                                               values_comb_mappedminus_transreg, values_comb_mappedminus_nontransreg)),
            "PhosphositeCount" : np.concatenate((counts_comb_ccd, counts_comb_nonccd, 
                                                      counts_comb_transreg, counts_comb_nontransreg)),
            "FractionPhosphositesDownstreamOfKinase" : np.concatenate((values_comb_ccd, values_comb_nonccd, 
                                                values_comb_transreg, values_comb_nontransreg)),
            "FisherPValue" : np.concatenate((fisher_comb_ccd[:,1], fisher_comb_nonccd[:,1], 
                                            fisher_comb_transreg[:,1], fisher_comb_nontransreg[:,1])),
            "FisherPValue_BonfCorrected" : pvals_corrected_bonf,
            "FisherPValue_BonfAlpha0.05Pass" : reject_bonf,
            }).to_csv("output/upstreamKinaseResults.csv", index = False)
        
