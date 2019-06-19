#%% Imports
from imports import *

from pybiomart import Server
server = Server("http://www.ensembl.org")
ensembl = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

#%% Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
adata, phases_filt = read_counts_and_phases()
qc_filtering(adata)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Read in anova results
anova_tests = pd.read_csv("output/transcript_regulation.csv")

#%% How many of the variable genes of each biotype are there?
    # pd.DataFrame([x for x in list(countsT_acceptedish.columns) if x not in accepted_symbols]).to_csv("output/genes_filtered_no_accepted_symbol")
    # countsT_accepted = countsT_acceptedish.loc[:, countsT_acceptedish.columns.isin(accepted_symbols)]

gene_biotypes = ensembl.query(attributes=[ "external_gene_name", "gene_biotype"])#, filters={"chromosome_name": ["22"]})
anova_tests_typed = anova_tests.merge(gene_biotypes, left_on="gene", right_on="Gene name", how="left")
anova_tests_typed.groupby("gene")["Gene type"].agg(["Gene type", ", ".join])

#%% This should be much easier with the RSEM dataframes
# Can possibly use the script already written to transfer GFF info to the table