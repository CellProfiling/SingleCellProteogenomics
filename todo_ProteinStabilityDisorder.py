#%% Imports
from imports import *
from Bio import SeqIO

#%% Import the genes names we're analyzing
ccd_transcript_regulated = np.array(pd.read_csv("output/ccd_transcript_regulated.csv")["gene"])
ccd_nontranscript_regulated = np.array(pd.read_csv("output/ccd_nontranscript_regulated.csv")["gene"])
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])

#%% Look at protein disorder
# Idea: 
#      Is there more protein disorder among proteins that are transcript CCD
#      Is there more protein disorder among proteins that have lower melting points?
# Execution: 
#      Use public tool for analyzing disordered regions
#      Import proteins with disorered regions / count disordered regions per protein
#      Maybe look at whether the disordered regions are at the terminus or not
# Output:
#      Disordered regions of proteins
#      Significant difference between transcript CCD and non-transcript CCD?