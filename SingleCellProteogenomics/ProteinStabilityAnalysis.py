# -*- coding: utf-8 -*-
"""
Analysis of mass spectrometry (MS) proteomic profiling of protein melting points:
    - Allows evaluation of protein stability, since lower melting points indicate higher propensity for unfolding
    - The results from three different cell lines (A549, HEK293, HepG2) were averaged for each protein

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

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