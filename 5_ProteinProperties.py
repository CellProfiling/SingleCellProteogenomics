# -*- coding: utf-8 -*-
"""
Investigation of the properties of proteins with different cell cycle regulation using PTMs and stability measurements
-  PTMs are observed using previously published bulk and phospho-enriched mass spectrometry (MS) proteomic data
-  Differences in PTM regulation is inferred using PTM occupancy for each PTM site
-  Protein stability was measured by MS thermal profiling in an external study
-  Differences in thermal shifts indicate different stabilities and propensity for unfolding

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import (FucciCellCycle, Loaders,
                                      ProteinPropertyAnalysis,
                                      RNADataPreparation, utils)
from SingleCellProteogenomics.utils import *

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300
fucci = FucciCellCycle.FucciCellCycle()
u_rna_plates = ["355","356","357"]

#%% Import the genes names we're analyzing
# Read in RNA-Seq data again and the CCD gene lists
valuetype, use_spikeins, biotype_to_use = "Tpms", False, "protein_coding"
adata, phases = RNADataPreparation.read_counts_and_phases(
    valuetype, use_spikeins, biotype_to_use, u_rna_plates
)
adata, phasesfilt = RNADataPreparation.qc_filtering(
    adata, do_log_normalize=True, do_remove_blob=True
)
adata = RNADataPreparation.zero_center_fucci(adata)
import_dict = Loaders.load_ptm_and_stability(adata)
wp_ensg, ccd_comp, nonccd_comp, ccdtranscript, wp_max_pol = (
    import_dict["wp_ensg"],
    import_dict["ccd_comp"],
    import_dict["nonccd_comp"],
    import_dict["ccdtranscript"],
    import_dict["wp_max_pol"],
)
ensg_results, name_results = utils.save_gene_names_by_category(
    adata, wp_ensg, ccd_comp, nonccd_comp, ccdtranscript
)
ensg_ccdtranscript, ensg_nonccdtranscript, ensg_ccdprotein, ensg_nonccdprotein, ensg_ccdprotein_transcript_regulated, ensg_ccdprotein_nontranscript_regulated, genes_analyzed, ccd_regev_filtered, ccd_filtered = ensg_results
names_ccdtranscript, names_nonccdtranscript, names_ccdprotein, names_nonccdprotein, names_ccdprotein_transcript_regulated, names_ccdprotein_nontranscript_regulated, names_genes_analyzed, names_ccd_regev_filtered, names_ccd_filtered = name_results
bioccd = np.genfromtxt(
    "input/ProteinData/BiologicallyDefinedCCD.txt", dtype="str"
)  # from mitotic structures
names_bioccd = utils.ccd_gene_names(bioccd, utils.getGeneNameDict())

#%% Analyze properties of the different groups relative to melting points
proteinProperties = ProteinPropertyAnalysis.ProteinProperties(
    wp_ensg,
    ensg_ccdprotein,
    ensg_ccdprotein_transcript_regulated,
    ensg_ccdprotein_nontranscript_regulated,
    bioccd,
    ensg_nonccdprotein,
    ensg_ccdtranscript,
    names_bioccd,
    names_ccdprotein,
    names_ccdprotein_transcript_regulated,
    names_ccdprotein_nontranscript_regulated,
    names_nonccdprotein,
    names_ccdtranscript,
)
proteinProperties.analyze_melting_points()
proteinProperties.analyze_disorder()
proteinProperties.analyze_length()
proteinProperties.statistical_properties_table()
proteinProperties.generate_properties_table()
proteinProperties.generate_statistical_boxplots()
proteinProperties.tm_scatters()
proteinProperties.kinase_families()
