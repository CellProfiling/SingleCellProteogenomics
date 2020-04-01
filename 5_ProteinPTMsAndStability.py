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
plate, valuetype, use_spikeins, biotype_to_use = "All", "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(plate, valuetype, use_spikeins, biotype_to_use)
adata, phasesfilt = RNADataPreparation.qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)

import_dict = Loaders.load_ptm_and_stability(adata)
wp_ensg, wp_max_pol = import_dict["wp_ensg"], import_dict["wp_max_pol"]
ccdtranscript, nonccdtranscript = import_dict["ccdtranscript_names"], import_dict["nonccdtranscript_names"]
ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated = import_dict["ccdprotein_transcript_regulated_names"], import_dict["ccdprotein_nontranscript_regulated_names"]
genes_analyzed = import_dict["genes_analyzed_names"]
ccd_regev_filtered, ccd_filtered = import_dict["ccd_regev_filtered_names"], import_dict["ccd_filtered_names"] 
ccdprotein, nonccdprotein = import_dict["ccdprotein_names"], import_dict["nonccdprotein_names"]
utils.save_gene_names_by_category(adata) # We're using gene names for this analysis instead of ENSG gene identifiers, so save those gene names

#%% Analyze the PTMs from bulk U2OS data and see if they are more expressed
# one or the other
# Idea: PTM annotation counts are pretty biased; PTM data might be better
# Execution: 
#   Take in the results from analyzing U2OS data with MetaMorpheus.
#   Count mods for each protein, excluding oxidations
# Output: # of PTMs per protein in each class and for all proteins

# Read in MetaMorpheus results
genemodsBulk = PTMAnalysis.analyze_ptms("input/raw/U2OSBulkAllProteinGroups.tsv", ccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, genes_analyzed)
genemodsPhospho = PTMAnalysis.analyze_ptms("input/raw/U2OSPhosphoAllProteinGroups.tsv", ccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, genes_analyzed)

# Analyze the modification occupancy for each PTM site
dfBulk, occdfBulk = PTMAnalysis.process_genemods(genemodsBulk)
dfPhospho, occdfPhospho = PTMAnalysis.process_genemods(genemodsPhospho)

# Investigate whether PTM occupancies have correlations with cell division pseudotime (they don't with this dataset)
PTMAnalysis.temporal_ptm_regulation_not_observed(dfBulk, occdfBulk, "Bulk", wp_max_pol, wp_ensg, ccdprotein)
PTMAnalysis.temporal_ptm_regulation_not_observed(dfPhospho, occdfPhospho, "Phospho", wp_max_pol, wp_ensg, ccdprotein)

# Compare the PTM regulation using PTM site occupancies between transcript and non-transcript regulated groups
PTMAnalysis.compare_ptm_regulation(dfBulk, occdfBulk, "Bulk", genes_analyzed, ccd_regev_filtered, ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, ccdprotein, nonccdprotein)
PTMAnalysis.compare_ptm_regulation(dfPhospho, occdfPhospho, "Phospho", genes_analyzed, ccd_regev_filtered, ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, ccdprotein, nonccdprotein)

#%% Perform melting point analysis
# Idea: A lower melting point for a protein indicates a higher propensity to unfold, 
#    and so the melting points measured by MS proteomics serve as a useful way to study the protein stability of CCD proteins
# Output: Boxplots illustrating differences in stability for different classes of CCD proteins
ProteinStabilityAnalysis.melting_point_analysis(ccdtranscript, nonccdtranscript, ccdprotein_transcript_regulated, ccdprotein_nontranscript_regulated, nonccdprotein)


