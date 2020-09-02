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
from SingleCellProteogenomics import utils, Loaders, FucciCellCycle, RNADataPreparation, PTMAnalysisPSMLevel, PTMAnalysisIntensityLevel, ProteinStabilityAnalysis
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
ensg_ccdtranscript, ensg_nonccdtranscript, ensg_ccdprotein, ensg_nonccdprotein, ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated, genes_analyzed, ccd_regev_filtered, ccd_filtered = name_results[0]
names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_genes_analyzed, names_ccd_regev_filtered, names_ccd_filtered = name_results[1]
bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
names_bioccd = utils.ccd_gene_names(bioccd, utils.getGeneNameDict())

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins
ptmAnalysis = PTMAnalysisPSMLevel.PTMAnalysis(names_ccdtranscript, 
                  names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, 
                  names_genes_analyzed, names_ccdprotein, names_ccd_regev_filtered, 
                  names_nonccdtranscript, names_nonccdprotein)

# Read in MetaMorpheus results
genemodsBulk = ptmAnalysis.analyze_ptms("input/raw/v305ProteinResults/bulk/search/AllProteinGroups.tsv.gz")
genemodsPhospho = ptmAnalysis.analyze_ptms("input/raw/v305ProteinResults/phospho/search/AllProteinGroups.tsv.gz")

# Analyze the modification occupancy for each PTM site
dfBulk, occdfBulk = ptmAnalysis.process_genemods(genemodsBulk)
dfPhospho, occdfPhospho = ptmAnalysis.process_genemods(genemodsPhospho)

# Investigate whether PTM occupancies have correlations with cell division pseudotime (they don't with this dataset)
ptmAnalysis.temporal_ptm_regulation_not_observed(dfBulk, occdfBulk, "Bulk", wp_max_pol, wp_ensg)
ptmAnalysis.temporal_ptm_regulation_not_observed(dfPhospho, occdfPhospho, "Phospho", wp_max_pol, wp_ensg)

# Compare the PTM regulation using PTM site occupancies between transcript and non-transcript regulated groups
ptmAnalysis.compare_ptm_regulation(dfBulk, occdfBulk, "Bulk")
ptmAnalysis.compare_ptm_regulation(dfPhospho, occdfPhospho, "Phospho")

#%% Do this analysis with the intensities
ptmAnalysis = PTMAnalysis.IntensityPTMAnalysis(ensg_ccdtranscript, ensg_nonccdtranscript, ensg_ccdprotein, ensg_nonccdprotein, ensg_ccdprotein_transcript_regulated, 
               ensg_ccdprotein_nontranscript_regulated, genes_analyzed, ccd_regev_filtered, ccd_filtered,
               names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, names_ccdprotein_transcript_regulated,
               names_ccdprotein_nontranscript_regulated, names_genes_analyzed, names_ccd_regev_filtered, names_ccd_filtered,
               bioccd)
# look for any significant trends in the Olsen 2010 cell cycle phospho data
ptmAnalysis.get_olsen_stoichiometry()

# look at the bulk data with MS1 MetaMorpheus intensities
result_bulk = ptmAnalysis.analyze_mm_ms1_ratios_bulk()
modRatios_mm_ms1_bulk, modificationsPerPeptide_bulk, targetsPep, unambigTargetPeps, unambigStarts = result_bulk

# look at the phospho data, compare MM and MQ
result_phospho_mm = ptmAnalysis.get_mm_ratios_phospho()
modRatios_mm_ms1_phospho, modPeptides_ms1_phospho, targetsPep, unambigTargetPeps, unambigStarts, modRatios_mm_ms2, modificationsPerPeptide_mm_ms2, seqToIntensities_mqMs2 = result_phospho_mm

result_phospho_mq = ptmAnalysis.get_mq_ratios_phospho()
modRatios_mq_ms1, modRatios_mq_ms2, modificationsPerPeptide_mq, evidenceAllPep = result_phospho_mq

#%% Perform melting point analysis
# Idea: A lower melting point for a protein indicates a higher propensity to unfold, 
#    and so the melting points measured by MS proteomics serve as a useful way to study the protein stability of CCD proteins
# Output: Boxplots illustrating differences in stability for different classes of CCD proteins
all_temps, all_protnames = ProteinStabilityAnalysis.melting_point_analysis(names_ccdtranscript, names_nonccdtranscript, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_nonccdprotein)

#%% Analyze properties of the different groups relative to melting points
proteinProperties = ProteinStabilityAnalysis.ProteinProperties(all_temps, all_protnames,
                            wp_ensg, ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated,
                            names_bioccd, names_ccdprotein,  names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_nonccdprotein, )
proteinProperties.analyze_disorder()
proteinProperties.analyze_cysteines()
proteinProperties.analyze_hydrophilic()
proteinProperties.analyze_hydrophobic()
proteinProperties.analyze_abundances()
proteinProperties.analyze_variants()
proteinProperties.temp_scatters()