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
    
def read_variants(vcffilepath):
    '''Reads VCF annotated by SnpEff'''
    anncomma = re.compile(b',(?!\()')
    transcript_variant = {}
    with gzip.open(vcffilepath, 'rb') as file:
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
    # maybe check that it's the same geneId too?
    gene_info = pd.read_csv(f"input/processed/transcriptId_geneName.txt", sep="\t", index_col=False)
    transcriptId_proteinName = dict((info[1]["Transcript stable ID"], info[1]["Gene stable ID"]) for info in gene_info.iterrows())
    geneId_variant = {}
    for item in transcript_variant.items():
        baseProteinName = transcriptId_proteinName[item[0].decode('ascii')]#.split('-')[0]
        if baseProteinName not in geneId_variant: geneId_variant[baseProteinName] = item[1]
        else: geneId_variant[baseProteinName].extend(item[1])
    return geneId_variant

class ProteinProperties:
    '''Analyzes protein properties for evaluating differences in melting temperatures'''
    def __init__(self, wp_ensg, ensg_ccdprotein, 
            ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated, 
            bioccd, ensg_nonccdprotein, ensg_ccdtranscript,
            names_bioccd, names_ccdprotein, 
            names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, 
            names_nonccdprotein, names_ccdtranscript):
        '''Reads in information about protein properties for evaluating differences in melting temps'''
        self.proteinDisorder = dict([(r[1][1], r[1][2:]) for r in pd.read_csv("input/processed/ProteinDisorderProperties.csv.gz").iterrows()])
        self.aebersoldNumbers = pd.read_csv("C:/Users/antho/Dropbox/Projects/Nucleoli/AebersoldNumbers.csv", index_col=False)
        self.geneId_variant = read_variants("C:/Users/antho/Dropbox/ProjectData/CellCycle/combined.spritz.snpeff.vcf.gz")
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
        self.names_hpamapped = pd.read_csv("input/raw/proteinatlas.tsv.gz", sep="\t")["Gene"]
        self.test_results = {}

    def analyze2(self, aa, bb, aa_lab, bb_lab, property_label, testZerosSeparately=False):
        if testZerosSeparately:
            self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}_NonzeroDistribution")] = ("Median", np.median(aa[aa!=0]), np.median(bb[bb!=0]), 
                    "Kruskal", stats.kruskal(aa[aa!=0], bb[bb!=0])[1])
            self.test_results[(property_label, f"{aa_lab}_eqvar_{bb_lab}_NonzeroDistribution")] = ("Std", np.std(aa[aa!=0]), np.std(bb[bb!=0]), 
                   "EqualVarianceLevene", stats.levene(aa[aa!=0], bb[bb!=0])[1])
            self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}_BinaryFeature")] = ("Count|CountZero", f"{sum(aa != 0)}|{sum(aa == 0)}", f"{sum(bb != 0)}|{sum(bb == 0)}",
                    "FisherExact", stats.fisher_exact([[sum(aa != 0), sum(aa == 0)], [sum(bb != 0), sum(bb == 0)]])[1])
        else:
            self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}")] = ("Median", np.median(aa), np.median(bb), 
                    "Kruskal", stats.kruskal(aa, bb)[1])
            self.test_results[(property_label, f"{aa_lab}_eqvar_{bb_lab}")] = ("Std", np.std(aa), np.std(bb), 
                   "EqualVarianceLevene", stats.levene(aa, bb)[1])

    def analyze3(self, aa, bb, cc, aa_lab, bb_lab, cc_lab, property_label, testZerosSeparately=False):#, lab_lab, showfliers, fileprefix):
        '''General boxplot and kruskal-wallis test'''
        # utils.general_boxplot((aa, bb, cc), (aa_lab, bb_lab, cc_lab), "", lab_lab, "", True, f"{fileprefix}_outliers.png")
        # utils.general_boxplot((aa, bb, cc), (aa_lab, bb_lab, cc_lab), "", lab_lab, "", False, f"{fileprefix}.png")
        self.analyze2(aa, bb, aa_lab, bb_lab, property_label, testZerosSeparately)
        self.analyze2(bb, cc, bb_lab, cc_lab, property_label, testZerosSeparately)
        self.analyze2(aa, cc, aa_lab, cc_lab, property_label, testZerosSeparately)
        
    def analyzeAll(self, property_label, *groups, testZerosSeparately=False):
        self.test_results[(property_label), "allgroups"] = ("N/A", "N/A", "N/A", 
                        "Kruskal", stats.kruskal(*[group[group != 0] if testZerosSeparately else group for group in groups])[1])
    
    def analyze_melting_points(self):
        '''Gather measurements of melting points for each protein across 10 cell lines; get the median per protein and evaluate for differences between CCD groups'''
        # Load melting points for each protein across 10 samples
        meltingDf = pd.read_csv("input/raw/ProteinStability/human.csv.gz")
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

    def analyzeSequenceFraction(self, idx, property_label, isLength=False, testZerosSeparately=False):
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
                      testZerosSeparately=testZerosSeparately)
        self.analyze3(fractPropertyCcd, fractPropertyNonCcd, fractPropertyAllMapped, 
                      "ccd", "nonccd", "allmapped", property_label, 
                      testZerosSeparately=testZerosSeparately)
        self.analyze2(fractPropertyVariable, fractPropertyAllMapped, "variable", "allmapped", property_label, testZerosSeparately=testZerosSeparately)
        if property_label == "Disorder":
            self.analyze3(fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyAllMapped, 
                          "mitoticCCD", "interphaseCCD", "allmapped", property_label, 
                          testZerosSeparately=testZerosSeparately)
            self.analyze2(fractPropertyCcdPseudotime, fractPropertyNonCcd, 
                          "interphaseCCD", "nonccd", property_label, testZerosSeparately=testZerosSeparately)
            self.analyze2(fractPropertyBioccd, fractPropertyNonCcd, 
                          "mitoticCCD", "nonccd", property_label, testZerosSeparately=testZerosSeparately)
        return fractPropertyAll, fractPropertyAllMapped, fractPropertyVariable, fractPropertyTransreg, fractPropertyNontransreg, fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyNonCcd, fractPropertyCcd, fractPropertyTranscriptCCD

    def analyze_disorder(self):
        disorderedMitotic = sum([self.proteinDisorder[nn][0] for nn in self.names_bioccd if nn in self.proteinDisorder])
        totalMitotic = sum([nn in self.proteinDisorder for nn in self.names_bioccd])
        print(f"{disorderedMitotic / totalMitotic}: fraction of disordered mitotic proteins")
        disorderedCcdProtein = sum([self.proteinDisorder[nn][0] for nn in self.names_ccdprotein if nn in self.proteinDisorder and not nn in self.names_bioccd])
        totalCcdProtein = sum([nn in self.proteinDisorder for nn in self.names_ccdprotein if not nn in self.names_bioccd])
        print(f"{disorderedCcdProtein / totalCcdProtein}: fraction of disordered CCD pseuduotime-only proteins")
        results = self.analyzeSequenceFraction(2, "Disorder", testZerosSeparately=True)
        self.fractDisorderAll, self.fractDisorderAllMapped, self.fractDisorderVariable, self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderBioccd, self.fractDisorderCcdPseudotime, self.fractDisorderNonCcd, self.fractDisorderCcd, self.fractDisorderTranscriptCCD = results
        
    def analyze_cysteines(self):
        results = self.analyzeSequenceFraction(3, "Cysteines", testZerosSeparately=True)
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

    def analyze_abundances(self):
        all_intensities = self.aebersoldNumbers["Copies Per Cell"]
        self.aebGenes = list(self.aebersoldNumbers["GeneName"])
        self.aebGenesSet = set(self.aebGenes)
        self.abundanceTransreg = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_ccdprotein_transcript_regulated if nn in self.aebGenesSet]
        self.abundanceNontransreg = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_ccdprotein_nontranscript_regulated if nn in self.aebGenesSet]
        self.abundanceBioccd = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_bioccd if nn in self.aebGenesSet]
        self.abundanceCcdPseudotime = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_ccdprotein if nn in self.aebGenesSet]
        self.abundanceNonccd = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_nonccdprotein if nn in self.aebGenesSet]
        self.abundanceCcd = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_ccdprotein if nn in self.aebGenesSet]
        self.abundanceCcdTranscript = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_ccdtranscript if nn in self.aebGenesSet]
        self.abundanceAllMapped = [all_intensities[self.aebGenes.index(nn)] for nn in self.names_hpamapped if nn in self.aebGenesSet]
        self.abundanceVariable = np.concatenate((self.abundanceCcd, self.abundanceNonccd))
        self.abundanceAll = all_intensities
        property_label = "Abundance"
        self.analyze3(self.abundanceTransreg, self.abundanceNontransreg, self.abundanceAllMapped, "transreg", "nontransreg", "allmapped", property_label)
        self.analyze3(self.abundanceCcd, self.abundanceNonccd, self.abundanceAllMapped, "ccd", "nonccd", "allmapped", property_label)
        self.analyze2(self.abundanceVariable, self.abundanceAllMapped, "variable", "allmapped", property_label)

    def test_variants(self, idsA, idsB, aa_lab, bb_lab, property_label):
        countWithVariantsA = sum([ensg in self.geneId_variant for ensg in idsA])
        countWithoutVariantsA = len(idsA) - countWithVariantsA
        countWithVariantsB = sum([ensg in self.geneId_variant for ensg in idsB])
        countWithoutVariantsB = len(idsB) - countWithVariantsB
        _, fisherPvalue = stats.fisher_exact([[countWithVariantsA, countWithoutVariantsA],
                                              [countWithVariantsB, countWithoutVariantsB]])
        #self.test_results # (comparisonName, propertyName),  (valuename, valueA, valueB, testName, testPvalue)
        self.test_results[(property_label, f"{aa_lab}_vs_{bb_lab}")] = ("CountWith|CountWithout", 
                    f"{countWithVariantsA}|{countWithoutVariantsA}", f"{countWithVariantsB}|{countWithoutVariantsB}", 
                    "FisherExact", fisherPvalue)
        
    def test_variants_vsallmapped(self, ensg_mapped, idsA, aa_lab, property_label):
        ensg_allminusA = ensg_mapped[~np.isin(ensg_mapped, idsA)]
        self.test_variants(idsA, ensg_allminusA, aa_lab, "allmapped", property_label)

    def analyze_variants(self):
        ensg_mapped = self.wp_ensg[np.isin(utils.ccd_gene_names_gapped(self.wp_ensg, utils.getGeneNameDict()), self.names_hpamapped)]
        self.test_variants_vsallmapped(ensg_mapped, self.ensg_ccdprotein, "ccd", "Variants")
        self.test_variants_vsallmapped(ensg_mapped, self.ensg_nonccdprotein, "ccd", "Variants")
        self.test_variants_vsallmapped(ensg_mapped, self.ensg_ccdprotein_transcript_regulated, "transreg", "Variants")
        self.test_variants_vsallmapped(ensg_mapped, self.ensg_ccdprotein_nontranscript_regulated, "nontransreg", "Variants")
        self.test_variants(self.ensg_ccdprotein, self.ensg_nonccdprotein, "ccd", "nonccd", "Variants")
        self.test_variants(self.ensg_ccdprotein_transcript_regulated, self.ensg_ccdprotein_nontranscript_regulated, "transreg", "nontransreg", "Variants")
        self.test_variants_vsallmapped(ensg_mapped, np.concatenate((self.ensg_ccdprotein, self.ensg_nonccdprotein)), "variable", "Variants")
    
    def get_property_and_temps(self, propertyValues, names_filter, isAbundance=False):
        '''Filters temperatures given a list of names'''
        self.all_protnames_set = set(self.all_protnames)
        self.all_protnames_list = list(self.all_protnames)
        if isAbundance:
            filterNames = [nn in self.all_protnames_set and nn in self.proteinDisorder for nn in self.aebGenes] 
            filteredPropertyValues = np.log10(propertyValues)[filterNames]
            filteredTempValues = [self.all_temps[self.all_protnames_list.index(nn)] for nn in self.aebGenes if nn in self.all_protnames_set  and nn in self.proteinDisorder]
        else:
            filterNames = [nn in self.all_protnames_set for nn in names_filter if nn in self.proteinDisorder] 
            filteredPropertyValues = np.asarray(propertyValues)[filterNames]
            filteredTempValues = [self.all_temps[self.all_protnames_list.index(nn)] for nn in names_filter if nn in self.all_protnames_set and nn in self.proteinDisorder]
        return filteredPropertyValues, filteredTempValues
        
    def temp_property_scatter(self, transregProperty, nontransregProperty, nonccdProperty, propertyLabel):
        # plt.scatter(*self.get_property_and_temps(allProperty, self.proteinDisorder.keys()), label="Non-CCD")
        # plt.xlabel(f"Fraction {propertyLabel} Residues")
        # plt.ylabel("Melting Point (°C)")
        # plt.legend()
        # plt.savefig(f"figures/{propertyLabel}VsTm_all.png")
        # plt.show(); plt.close()
        
        plt.scatter(*self.get_property_and_temps(transregProperty, self.names_ccdprotein_transcript_regulated), label="Trans Reg")
        plt.scatter(*self.get_property_and_temps(nontransregProperty, self.names_ccdprotein_nontranscript_regulated), label="Non-Trans Reg")
        plt.scatter(*self.get_property_and_temps(nonccdProperty, self.names_nonccdprotein), label="Non-CCD")
        plt.xlabel(f"Fraction {propertyLabel} Residues")
        plt.ylabel("Melting Point (°C)")
        plt.legend()
        plt.savefig(f"figures/{propertyLabel}VsTm.png")
        plt.show(); plt.close()
    
    def apply_linregress(self, allProperty, propertyLabel, isAbundance = False):
        names = self.proteinDisorder.keys() if not isAbundance else self.aebGenes
        propertyValues, propertyTemps = self.get_property_and_temps(allProperty, names, isAbundance)
        linModel = stats.linregress(propertyValues, propertyTemps)
        plt.scatter(propertyValues, propertyTemps, alpha=0.1)
        xfit = np.min(propertyValues) + np.arange(100) / 100 * (np.max(propertyValues) - np.min(propertyValues))
        yfit = linModel.intercept + xfit * linModel.slope
        plt.scatter(xfit, yfit)
        plt.xlabel(f"Fraction {propertyLabel} Residues")
        plt.ylabel("Melting Point (°C)")
        plt.savefig(f"figures/All{propertyLabel}VsTm.png")
        plt.show(); plt.close()
        print(f"{linModel.slope}: slope for all melting temperatures vs fract {propertyLabel} residues")
        print(f"{linModel.rvalue**2}: r-squared for all melting temperatures vs fract {propertyLabel} residues")
        print(f"{linModel.pvalue}: p-value for nonzero slope for all melting temperatures vs fract {propertyLabel} residues")
    
    def tm_scatters(self):
        self.temp_property_scatter(self.fractHydrophobicTransreg, self.fractHydrophobicNontransreg, self.fractHydrophobicNonCcd, "Hydrophobic")
        self.temp_property_scatter(self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderNonCcd, "Disorder")
        self.apply_linregress(self.fractHydrophobicAll, "Hydrophobic")
        self.apply_linregress(self.fractDisorderAll, "Disorder")
        self.apply_linregress(self.fractCysteinesAll, "Cysteines")
        self.apply_linregress(self.lengthAll, "Length")
        self.apply_linregress(self.abundanceAll, "Log10 Abundance", True)
        
    def generate_properties_table(self):
        '''Make a table with all the values analyzed'''
        ensg_variant = np.array(list(self.geneId_variant.items()), dtype=object)
        variant_gene_names = np.array(utils.ccd_gene_names_gapped(self.geneId_variant.keys(), utils.getGeneNameDict()))
        variant_gene_name_set = set(variant_gene_names)
        protDisorderGeneList = list(self.proteinDisorder.keys())
        genes = np.sort(np.unique(np.concatenate((protDisorderGeneList, self.aebGenes, self.all_protnames_list, variant_gene_names))))
        lengths, disordered, hydrophob, cysteine, polar, variants = [],[],[],[],[],[]
        for nn in genes:
            isGene = nn in self.proteinDisorder
            idx = -1 if not isGene else protDisorderGeneList.index(nn)                
            lengths.append("" if not isGene else self.lengthAll[idx])
            disordered.append("" if not isGene else self.fractDisorderAll[idx])
            hydrophob.append("" if not isGene else self.fractHydrophobicAll[idx])
            cysteine.append("" if not isGene else self.fractCysteinesAll[idx])
            polar.append("" if not isGene else self.fractPolarAll[idx])
            variants.append("" if not nn in variant_gene_name_set else ";".join(
                ["_".join((xx[0], "_".join([vv.decode('ascii') for vv in xx[1]]))) for xx in ensg_variant[variant_gene_names == nn]]))
        pd.DataFrame({
            "Protein" : genes,
            "IsCCD" : np.isin(genes, self.names_ccdprotein),
            "IsNonCCD" : np.isin(genes, self.names_nonccdprotein),
            "IsInterphaseCCD" : np.isin(genes, self.names_ccdprotein) & ~np.isin(genes, self.names_bioccd),
            "IsMitoticCCD" : np.isin(genes, self.names_bioccd),
            "IsTranscriptRegCCD" : np.isin(genes, self.names_ccdprotein_transcript_regulated),
            "IsNonTranscriptRegCCD" : np.isin(genes, self.names_ccdprotein_nontranscript_regulated),
            "Melting Temperature (deg C)" : [self.all_temps[self.all_protnames_list.index(nn)] if nn in self.all_protnames_set else "" for nn in genes],
            "Length" : lengths,
            "Abundance (Copies per Cell)" : [self.abundanceAll[self.aebGenes.index(nn)] if nn in self.aebGenesSet else "" for nn in genes],
            "Disordered Residue Fraction" : disordered,
            "Hydrophobic Residue Fraction" : hydrophob,
            "Cysteine Fraction" : cysteine,
            "Polar Residue Fraction" : polar,
            "Variants" : variants
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
        plt.show()
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
        
        # Abundance ()
        values = (np.log10(self.abundanceAllMapped), np.log10(self.abundanceVariable), 
                       np.log10(self.abundanceNonccd), np.log10(self.abundanceCcd), 
                       np.log10(self.abundanceTransreg), np.log10(self.abundanceNontransreg))
        self.boxplot_summary(*values, "Log10 Abundance", False)
        utils.barplot_annotate_brackets(0, 2, f"p={'%.1e' % self.test_results['Abundance','nonccd_vs_allmapped'][4]}", np.arange(len(values)), [max(vv) for vv in values]) # highlight CCD vs nonccd
        utils.barplot_annotate_brackets(2, 3, f"p={'%.1e' % self.test_results['Abundance','ccd_vs_nonccd'][4]}", np.arange(len(values)), [1+max(vv) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose("figures/ProteinAbundancePointBox.pdf")
        
        # Disorder ()
        values = (self.fractDisorderAllMapped, self.fractDisorderVariable, 
                       self.fractDisorderNonCcd, self.fractDisorderCcd, 
                       self.fractDisorderTransreg, self.fractDisorderNontransreg)
        self.boxplot_summary(*values, "Disorder", False)
        utils.barplot_annotate_brackets(0, 1, f"p={'%.1e' % self.test_results['Disorder','variable_vs_allmapped_NonzeroDistribution'][4]}", np.arange(len(values)), [max(vv) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/DisorderBox.pdf")
        
        # Hydrophobicity ()
        values = (self.fractHydrophobicAllMapped, self.fractHydrophobicVariable, 
                       self.fractHydrophobicNonCcd, self.fractHydrophobicCcd, 
                       self.fractHydrophobicTransreg, self.fractHydrophobicNontransreg)
        self.boxplot_summary(*values, "Fraction Hydrophobic", False)
        utils.barplot_annotate_brackets(0, 4, f"p={'%.1e' % self.test_results['Hydrophobic','ccd_vs_allmapped'][4]}", np.arange(len(values)), [max(vv) for vv in values]) # highlight CCD vs nonccd
        utils.barplot_annotate_brackets(3, 4, f"p={'%.1e' % self.test_results['Hydrophobic','ccd_vs_nonccd'][4]}", np.arange(len(values)), [0.2+max(vv) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/HydrophobicBox.pdf")
        
        # Length ()
        values = (self.lengthAllMapped, self.lengthVariable, 
                       self.lengthNonCcd, self.lengthCcd, 
                       self.lengthTransreg, self.lengthNontransreg)
        self.boxplot_summary(*values, "Log2 Length", False)
        utils.barplot_annotate_brackets(0, 1, f"p={'%.1e' % self.test_results['Log2 Length','variable_vs_allmapped'][4]}", np.arange(len(values)), [np.percentile(vv, 98) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/LengthBox.pdf")
        
        # Cysteines
        values = (self.fractCysteinesAllMapped, self.fractCysteinesVariable, 
                       self.fractCysteinesNonCcd, self.fractCysteinesCcd, 
                       self.fractCysteinesTransreg, self.fractCysteinesNontransreg)
        self.boxplot_summary(*values, "Cysteines", False)
        utils.barplot_annotate_brackets(0, 1, f"p={'%.1e' % self.test_results['Cysteines','variable_vs_allmapped_NonzeroDistribution'][4]}", np.arange(len(values)), [np.percentile(vv, 95) for vv in values]) # highlight CCD vs nonccd
        self.boxplot_saveclose(f"figures/CysteinesBox.pdf")
                
    
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
        kinaseFams = pd.read_csv("C:\\Users\\antho\\Dropbox\\ProjectData\\PhosphositePlus\\KinHubKinaseFamilies.csv")
        kinFamDict = dict([(row[7], row[4]) for idx, row in kinaseFams.iterrows()])
        phosphositeplus = pd.read_csv("C:\\Users\\antho\\Dropbox\\ProjectData\\PhosphositePlus\\Kinase_Substrate_Dataset_NoHeader.gz", sep="\t")
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
        
