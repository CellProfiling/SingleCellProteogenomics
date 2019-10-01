#%% Imports
from imports import *

#%% Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
dd = "All"
counts_or_rpkms = "Tpms"
do_log_normalization = True
# adata, phases = read_counts_and_phases(dd, counts_or_rpkms, False, "protein_coding")
# phases_filt = qc_filtering(adata, do_log_normalization)
# ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%%
biotype_to_use="protein_coding"
percent_variance_tests = pd.read_csv(f"output/transcript_regulation{biotype_to_use}.csv")
percent_variance_tests[percent_variance_tests["diana_ccd"] == True]["name"].to_csv("genewalk/ccd_proteins.csv", index=False)
percent_variance_tests[percent_variance_tests["diana_nonccd"] == True]["name"].to_csv("genewalk/nonccd_proteins.csv", index=False)

# rejectBonf_unsorted = percent_variance_tests["reject_B"]
# ccd_transcript_regulated = np.array(adata.var_names)[(rejectBonf_unsorted) & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff)]
# dianaccd_transcript_regulated = np.array(adata.var_names)[(dianaccdgenes) & (rejectBonf_unsorted) & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff)]
# dianaccd_nontranscript_regulated = np.array(adata.var_names)[(dianaccdgenes) & ~((rejectBonf_unsorted) & (percent_ccd_variance > percent_var_cutoff) & (total_variance > total_var_cutoff))]
# ccd_transcript_regulated.to_csv("genewalk/ccd_genes.csv", index=False)
# dianaccd_transcript_regulated.to_csv("genewalk/ccd_transcriptreg_proteins.csv", index=False)
# dianaccd_nontranscript_regulated.to_csv("genewalk/ccd_nontranscriptreg_proteins.csv", index=False)

#%%
# genewalk --project proteinccd --genes ccd_proteins.csv --id_type hgnc_symbol --nproc 4 --base_folder ./ccd_proteins
# genewalk --project nonccdprotein --genes nonccd_proteins.csv --id_type hgnc_symbol --nproc 4 --base_folder ./nonccd_proteins

#%%
