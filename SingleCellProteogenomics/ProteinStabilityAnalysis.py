# -*- coding: utf-8 -*-
"""
Analysis of mass spectrometry (MS) proteomic profiling of protein melting points:
    - Allows evaluation of protein stability, since lower melting points indicate higher propensity for unfolding
    - The results from three different cell lines (A549, HEK293, HepG2) were averaged for each protein

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

import re, gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from SingleCellProteogenomics import utils

def avg_prot_temps(temps, proteinNames):
    '''Finds median of the measurements from the same protein'''
    aaa = []
    ppp = []
    for i in range(len(temps)):
        df = pd.DataFrame({"name" : proteinNames[i], "temps": temps[i]})
        med = df.groupby("name")["temps"].median()
        aaa.append(list(med))
        ppp.append(list(med.index))
    result = aaa
    result.extend(ppp)
    return result

def stat_tests(title, all_temps, transcript_reg, nontranscr_reg, nonccd_temps, allccdtranscript):
    '''Generate text describing statistical tests'''
    results = ""
    results += f"{title} statistical tests\n"
    results += "Testing whether transcript reg. CCD melting points are different than non-transcript reg. CCD\n"
    results += f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)\n"
    results += f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)\n"
    results += f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(transcript_reg, nontranscr_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(transcript_reg, nontranscr_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether non-transcript reg. CCD is different than all proteins\n"
    results += f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)\n"
    results += f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(all_temps)}: All Proteins (median)\n"
    results += f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(all_temps, nontranscr_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(all_temps, nontranscr_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether transcript reg. CCD is different than all proteins\n"
    results += f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)\n"
    results += f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(all_temps)}: All Proteins (median)\n"
    results += f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(all_temps, transcript_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(all_temps, transcript_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether transcript reg. CCD is different than all transcript CCD\n"
    results += f"{np.mean(allccdtranscript)} +/- {np.std(allccdtranscript)}: all transcript CCD Proteins (mean, std)\n"
    results += f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(allccdtranscript)}: all transcript CCD Proteins (median)\n"
    results += f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(allccdtranscript, transcript_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(allccdtranscript, transcript_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether non-transcript reg. CCD is different than all transcript CCD\n"
    results += f"{np.mean(allccdtranscript)} +/- {np.std(allccdtranscript)}: all transcript CCD Proteins (mean, std)\n"
    results += f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(allccdtranscript)}: all transcript CCD Proteins (median)\n"
    results += f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(allccdtranscript, nontranscr_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(allccdtranscript, nontranscr_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether transcript reg. CCD is different than non-CCD\n"
    results += f"{np.mean(nonccd_temps)} +/- {np.std(nonccd_temps)}: Non-CCD Proteins (mean, std)\n"
    results += f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(nonccd_temps)}: Non-CCD Proteins (median)\n"
    results += f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(nonccd_temps, transcript_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(nonccd_temps, transcript_reg)}: two sided kruskal\n"
    results += "\n"
    results += "Testing whether nontranscript reg. CCD is different than non-CCD\n"
    results += f"{np.mean(nonccd_temps)} +/- {np.std(nonccd_temps)}: Non-CCD Proteins (mean, std)\n"
    results += f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)\n"
    results += f"{np.median(nonccd_temps)}: Non-CCD Proteins (median)\n"
    results += f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)\n"
    results += f"{stats.ttest_ind(nonccd_temps, nontranscr_reg)}: two sided t-test\n"
    results += f"{stats.kruskal(nonccd_temps, nontranscr_reg)}: two sided kruskal\n"
    results += "\n"
    return results

def temp_hist(title, all_temps, transcript_reg, nontranscr_reg, nonccd_temps, allccdtranscript):
    '''Generates a histogram of melting points with bins normalized to 1'''
    bins=np.histogram(np.hstack((all_temps, transcript_reg, nontranscr_reg, nonccd_temps)), bins=40)[1] #get the bin edges
    plt.hist(all_temps, bins=bins, weights=utils.weights(all_temps), color="#3753A4", alpha=0.5, label="All Proteins")
    plt.hist(transcript_reg, bins=bins, weights=utils.weights(transcript_reg), color="#FF0000", alpha=0.5, label="Transcript Reg. CCD")
    plt.hist(nontranscr_reg, bins=bins, weights=utils.weights(nontranscr_reg), color="#2CE100", alpha=0.6, label="Non-Transcript Reg. CCD")
    plt.hist(nonccd_temps, bins=bins, weights=utils.weights(nonccd_temps), color="#2CE100", alpha=0.6, label="Non-CCD")
    plt.legend(loc="upper right")
    plt.xlabel("Melting Point (°C)")
    plt.title(title)
    return stat_tests(title, all_temps, transcript_reg, nontranscr_reg, nonccd_temps, allccdtranscript)

def melting_point_analysis(ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, nonccdprotein):
    '''Gather measurements of melting points for each protein across 10 cell lines; get the median per protein and evaluate for differences between CCD groups'''
    # Load melting points for each protein across 3 cell lines
    all_temps, allnonccdtranscript, allccdtranscript, transcript_reg, nontranscr_reg, nonccd_temps = [],[],[],[],[],[]
    atp, ant, att, trp, ntp, nnp = [],[], [], [], [],[]
    ccd_groups = [[True] * len(ccdtranscript), ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, nonccdprotein]
    temp_groups = [all_temps, allccdtranscript, allnonccdtranscript, transcript_reg, nontranscr_reg, nonccd_temps]
    name_groups = [atp, att, ant, trp, ntp, nnp]
    
    meltingDf = pd.read_csv("input/raw/ProteinStability/human.csv.gz")
    print(f"{','.join(np.unique(meltingDf['cell_line_or_type']))}: unique human cell samples\n")
    
    notna = pd.notna(meltingDf["quan_norm_meltPoint"])
    meltingDfNotNa = meltingDf[notna]
    med = meltingDfNotNa.groupby("gene_name")["quan_norm_meltPoint"].median()
    all_temps, atp = np.asarray(med), np.asarray(med.index)
    allccdtranscript_temps, att = np.asarray(med[np.isin(med.index, list(ccdtranscript))]), np.asarray(med.index[np.isin(med.index, list(ccdtranscript))])
    allnonccdtranscript_temps, ant = np.asarray(med[np.isin(med.index, list(nonccdtranscript))]), np.asarray(med.index[np.isin(med.index, list(nonccdtranscript))])
    ccdprotein_transcript_regulated_temps, trp = np.asarray(med[np.isin(med.index, list(ccdprotein_transcript_regulated))]), np.asarray(med.index[np.isin(med.index, list(ccdprotein_transcript_regulated))])
    ccdprotein_nontranscript_regulated_temps, ntp = np.asarray(med[np.isin(med.index, list(ccdprotein_nontranscript_regulated))]), np.asarray(med.index[np.isin(med.index, list(ccdprotein_nontranscript_regulated))])
    nonccdprotein_temps, nnp = np.asarray(med[np.isin(med.index, list(nonccdprotein))]), np.asarray(med.index[np.isin(med.index, list(nonccdprotein))])
    
    # Aggregate measurements per protein
    statsresults = temp_hist("Aggregated Melting Points", 
         all_temps, ccdprotein_transcript_regulated_temps, ccdprotein_nontranscript_regulated_temps, nonccdprotein_temps, allccdtranscript_temps)
    plt.savefig("figures/ProteinMeltingPointsMedianed.png")
    plt.show()
    plt.close()
    pd.DataFrame({
        "protein_name":atp,
        "median_melting_point":all_temps,
        "transcript_reg":np.isin(atp, trp),
        "nontranscr_reg":np.isin(atp, ntp)}).to_csv("output/MedianMeltingPoints.csv",index=False)

    # Boxplots
    utils.general_boxplot((all_temps, ccdprotein_transcript_regulated_temps, ccdprotein_nontranscript_regulated_temps, nonccdprotein_temps),
          ("All Proteins", "Transcript's\nReg\nCCD", "Non-Transcript\nReg\nCCD", "Non-CCD"), 
        "Protein Set", "Melting Point (°C)", "", True, "figures/ProteinMeltingPointBox.pdf")
    utils.boxplot_with_stripplot((all_temps, allccdtranscript_temps, ccdprotein_nontranscript_regulated_temps, ccdprotein_transcript_regulated_temps, ccdprotein_nontranscript_regulated_temps, nonccdprotein_temps), 
        ("All Proteins", "All\nTranscript\nCCD", "All\nTranscript\nNon-CCD", "Transcript\nReg\nCCD", "Non-Transcript\nReg\nCCD", "Non-CCD"), 
        "", "Melting Point (°C)", "", True, "figures/ProteinMeltingPointBoxSelect.pdf")
    utils.general_boxplot((all_temps, ccdprotein_transcript_regulated_temps, ccdprotein_nontranscript_regulated_temps), ("All Proteins", "Transcript\nRegulated\nCCD Proteins", "Non-Transcript\nRegulated\nCCD Proteins"), 
        "", "Melting Point (°C)", "", False, "figures/ProteinMeltingPointBoxSelect2.pdf")

    print(statsresults)
    with open("output/meltingpointresults.txt", 'w') as file:
        file.write(statsresults)

    # Pickle the results
    utils.np_save_overwriting("output/temperatures.all_temps.npy", all_temps)
    utils.np_save_overwriting("output/temperatures.all_temp_prot.npy", atp)
    utils.np_save_overwriting("output/temperatures.transcript_reg.npy", transcript_reg)
    utils.np_save_overwriting("output/temperatures.transcript_reg_prot.npy", trp)
    utils.np_save_overwriting("output/temperatures.nontranscr_reg.npy", nontranscr_reg)
    utils.np_save_overwriting("output/temperatures.nontranscript_reg_prot.npy", ntp)
    utils.np_save_overwriting("output/temperatures.nonccd_temps.npy", nonccd_temps)
    utils.np_save_overwriting("output/temperatures.nonccd_temps_prot.npy", nnp)
    
    return all_temps, atp
    
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
    def __init__(self, all_temps, all_protnames,
            wp_ensg, ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated,
            names_bioccd, names_ccdprotein,  names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_nonccdprotein):
        '''Reads in information about protein properties for evaluating differences in melting temps'''
        self.proteinDisorder = dict([(r[1][1], r[1][2:]) for r in pd.read_csv("input/processed/ProteinDisorderProperties.csv.gz").iterrows()])
        self.aebersoldNumbers = pd.read_csv("C:/Users/antho/Dropbox/Projects/Nucleoli/AebersoldNumbers.csv", index_col=False)
        self.geneId_variant = read_variants("input/processed/combined.snpeff.vcf.gz")
        self.all_temps, self.all_protnames = all_temps, all_protnames
        self.wp_ensg, self.ensg_ccdprotein_transcript_regulated, self.ensg_ccdprotein_nontranscript_regulated = wp_ensg, ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated
        self.names_bioccd, self.names_ccdprotein, self.names_ccdprotein_transcript_regulated, self.names_ccdprotein_nontranscript_regulated, self.names_nonccdprotein = names_bioccd, names_ccdprotein,  names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_nonccdprotein

    def analyze(self, a, b, c, a_lab, b_lab, c_lab, lab_lab, showfliers, fileprefix):
        '''General '''
        print(f"{np.median(a)}: median of {a_lab}; {np.median(b)}: median of {b_lab}; {np.median(c)}: median of {c_lab}; ")
        utils.general_boxplot((a, b, c), (a_lab, b_lab, c_lab), "", lab_lab, "", True, f"{fileprefix}_outliers.png")
        utils.general_boxplot((a, b, c), (a_lab, b_lab, c_lab), "", lab_lab, "", False, f"{fileprefix}.png")
        print(f"{stats.kruskal(a, b)[1]}: {a_lab} vs {b_lab}, {lab_lab} kruskal p-value")
        print(f"{stats.kruskal(b, c)[1]}: {b_lab} vs {c_lab}, {lab_lab} kruskal p-value")
        print(f"{stats.kruskal(a, c)[1]}: {a_lab} vs {c_lab}, {lab_lab} kruskal p-value")
        print(f"{stats.kruskal(a, b, c)[1]}: {a_lab}, {b_lab}, {c_lab}, {lab_lab} kruskal p-value")

    def analyze_property(self, idx, property_label):
        protLen_idx = 6
        fractPropertyAll = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.proteinDisorder.keys()]
        fractPropertyTransreg = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_ccdprotein_transcript_regulated if nn in self.proteinDisorder]
        fractPropertyNontransreg = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_ccdprotein_nontranscript_regulated if nn in self.proteinDisorder]
        self.analyze(fractPropertyTransreg, fractPropertyNontransreg, fractPropertyAll, "transreg", "nontransreg", "all", f"Fraction {property_label}", True, f"figures/Fract{property_label}.png")
        
        fractPropertyBioccd = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_bioccd if nn in self.proteinDisorder]
        fractPropertyCcdPseudotime = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_ccdprotein if nn in self.proteinDisorder and nn not in self.names_bioccd]
        fractPropertyNonCcd = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_nonccdprotein if nn in self.proteinDisorder]
        fractPropertyCcd = [self.proteinDisorder[nn][idx] / self.proteinDisorder[nn][protLen_idx] for nn in self.names_ccdprotein if nn in self.proteinDisorder]
        self.analyze(fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyAll, "mitotic", "interphase", "all", f"Fraction {property_label}", True, f"figures/Fract{property_label}Mitotic.png")
        self.analyze(fractPropertyCcd, fractPropertyNonCcd, fractPropertyAll, "ccd", "nonccd", "all", f"Fraction {property_label}", True, f"figures/Fract{property_label}Ccd.png")
        self.analyze(fractPropertyCcdPseudotime, fractPropertyNonCcd, fractPropertyAll, "ccdInterphase", "nonccd", "all", f"Fraction {property_label}", True, f"figures/Fract{property_label}CcdInterphase.png")
        self.analyze(fractPropertyBioccd, fractPropertyNonCcd, fractPropertyAll, "ccdMitotic", "nonccd", "all", f"Fraction {property_label}", True, f"figures/Fract{property_label}CcdMitotic.png")
        return fractPropertyAll, fractPropertyTransreg, fractPropertyNontransreg, fractPropertyBioccd, fractPropertyCcdPseudotime, fractPropertyNonCcd, fractPropertyCcd
    
    def analyze_disorder(self):
        disorderedMitotic = sum([self.proteinDisorder[nn][0] for nn in self.names_bioccd if nn in self.proteinDisorder])
        totalMitotic = sum([nn in self.proteinDisorder for nn in self.names_bioccd])
        print(f"{disorderedMitotic / totalMitotic}: fraction of disordered mitotic proteins")
        disorderedCcdProtein = sum([self.proteinDisorder[nn][0] for nn in self.names_ccdprotein if nn in self.proteinDisorder and not nn in self.names_bioccd])
        totalCcdProtein = sum([nn in self.proteinDisorder for nn in self.names_ccdprotein if not nn in self.names_bioccd])
        print(f"{disorderedCcdProtein / totalCcdProtein}: fraction of disordered CCD pseuduotime-only proteins")
        results = self.analyze_property(2, "Disorder")
        self.fractDisorderAll, self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderBioccd, self.fractDisorderCcdPseudotime, self.fractDisorderNonCcd, self.fractDisorderCcd = results
        
    def analyze_cysteines(self):
        results = self.analyze_property(3, "Cysteines")
        self.fractCysteinesAll, self.fractCysteinesTransreg, self.fractCysteinesNontransreg, self.fractCysteinesBioccd, self.fractCysteinesCcdPseudotime, self.fractCysteinesNonCcd, self.fractCysteinesCcd = results

    def analyze_hydrophilic(self):
        results = self.analyze_property(4, "Hydrophilic")
        self.fractHydrophilicAll, self.fractHydrophilicTransreg, self.fractHydrophilicNontransreg, self.fractHydrophilicBioccd, self.fractHydrophilicCcdPseudotime, self.fractHydrophilicNonCcd, self.fractHydrophilicCcd = results

    def analyze_hydrophobic(self):
        results = self.analyze_property(5, "Hydrophobic")
        self.fractHydrophobicAll, self.fractHydrophobicTransreg, self.fractHydrophobicNontransreg, self.fractHydrophobicBioccd, self.fractHydrophobicCcdPseudotime, self.fractHydrophobicNonCcd, self.fractHydrophobicCcd = results

    def analyze_abundances(self):
        all_intensities = self.aebersoldNumbers["Copies Per Cell"]
        aebGenes = list(self.aebersoldNumbers["GeneName"])
        aebGenesSet = set(aebGenes)
        abundanceTransreg = [all_intensities[aebGenes.index(nn)] for nn in self.names_ccdprotein_transcript_regulated if nn in aebGenesSet]
        abundanceNontransreg = [all_intensities[aebGenes.index(nn)] for nn in self.names_ccdprotein_nontranscript_regulated if nn in aebGenesSet]
        abundanceBioccd = [all_intensities[aebGenes.index(nn)] for nn in self.names_bioccd if nn in aebGenesSet]
        abundanceCcdPseudotime = [all_intensities[aebGenes.index(nn)] for nn in self.names_ccdprotein if nn in aebGenesSet]
        abundanceNonccd = [all_intensities[aebGenes.index(nn)] for nn in self.names_nonccdprotein if nn in aebGenesSet]
        abundanceCcd = [all_intensities[aebGenes.index(nn)] for nn in self.names_ccdprotein if nn in aebGenesSet]
        abundanceAll = all_intensities
        self.analyze(abundanceTransreg, abundanceNontransreg, abundanceAll, "transreg", "nontransreg", "all", "Abundance", False, "figures/fractDisordered.png")
        self.analyze(abundanceBioccd, abundanceCcdPseudotime, abundanceAll, "mitotic", "interphase", "all", "Abundance", False, "figures/fractDisorderedMitotic.png")
        self.analyze(abundanceCcd, abundanceNonccd, abundanceAll, "ccd", "nonccd", "all",  "Abundance", False, "figures/fractDisorderedCcd.png")

    def analyze_variants(self):
        fractWithVariant_all = sum([ensg in self.geneId_variant for ensg in self.wp_ensg]) / len(self.wp_ensg)
        fractWithVariant_nonTransReg = sum([ensg in self.geneId_variant for ensg in self.ensg_ccdprotein_nontranscript_regulated]) / len(self.ensg_ccdprotein_nontranscript_regulated)
        fractWithVariant_transReg = sum([ensg in self.geneId_variant for ensg in self.ensg_ccdprotein_transcript_regulated]) / len(self.ensg_ccdprotein_transcript_regulated)
        print(f"{fractWithVariant_all}: fraction of all genes with variant")
        print(f"{fractWithVariant_nonTransReg}: fraction of nontranscript regulated genes with variant")
        print(f"{fractWithVariant_transReg}: fraction of transcript regulated genes with variant")
    
    def get_property_and_temps(self, propertyValues, names_filter):
        '''Filters temperatures given a list of names'''
        all_protnames_set, all_protnames_list = set(self.all_protnames), list(self.all_protnames)
        filterNames = [nn in all_protnames_set for nn in names_filter if nn in self.proteinDisorder]
        filteredPropertyValues = np.asarray(propertyValues)[filterNames]
        filteredTempValues = [self.all_temps[all_protnames_list.index(nn)] for nn in names_filter if nn in all_protnames_set and nn in self.proteinDisorder]
        return filteredPropertyValues, filteredTempValues
        
    def temp_property_scatter(self, transregProperty, nontransregProperty, nonccdProperty, propertyLabel):
        plt.scatter(*self.get_property_and_temps(transregProperty, self.names_ccdprotein_transcript_regulated), label="Trans Reg")
        plt.scatter(*self.get_property_and_temps(nontransregProperty, self.names_ccdprotein_nontranscript_regulated), label="Non-Trans Reg")
        plt.scatter(*self.get_property_and_temps(nonccdProperty, self.names_nonccdprotein), label="Non-CCD")
        plt.xlabel(f"Fraction {propertyLabel} Residues")
        plt.ylabel("Melting Point (°C)")
        plt.legend()
        plt.savefig(f"figures/{propertyLabel}VsTm.png")
        plt.show(); plt.close()
    
    def temp_scatters(self):
        self.temp_property_scatter(self.fractHydrophilicTransreg, self.fractHydrophilicNontransreg, self.fractHydrophilicNonCcd, "Hydrophilic")
        self.temp_property_scatter(self.fractDisorderTransreg, self.fractDisorderNontransreg, self.fractDisorderNonCcd, "Disorder")
        
        linModel = stats.linregress(*self.get_property_and_temps(self.fractHydrophilicAll, self.proteinDisorder.keys()))
        plt.scatter(*self.get_property_and_temps(self.fractHydrophilicAll, self.proteinDisorder.keys()), alpha=0.1)
        xfit = 0.3 + np.arange(100) / 100 * (0.8 - 0.3)
        yfit = linModel.intercept + xfit * linModel.slope
        plt.scatter(xfit, yfit)
        plt.xlabel("Fraction Hydrophilic Residues")
        plt.ylabel("Melting Point (°C)")
        plt.savefig("figures/AllHydryophilicVsTm.png")
        plt.show(); plt.close()
        print(f"{linModel.rvalue**2}: r-squared for all melting temperatures vs fract hydrophilic residues")