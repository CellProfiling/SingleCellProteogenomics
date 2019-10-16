#%% Imports
from imports import *
from Bio import SeqIO
import sys
import re
import math

#%% Take the results out of the pickle jars
modocc = pd.read_pickle("output/modocc.pkl")
all_modocc = modocc["all_modocc"]
ccd_t_modocc = modocc["ccd_t_modocc"]
ccd_rt_modocc = modocc["ccd_rt_modocc"]
ccd_at_modocc = modocc["ccd_at_modocc"]
ccd_n_modocc = modocc["ccd_n_modocc"]
nonccd_modocc = modocc["nonccd_modocc"]

all_modocc_genes = modocc["all_modocc_genes"]
ccd_t_modocc_genes = modocc["ccd_t_modocc_genes"]
ccd_rt_modocc_genes = modocc["ccd_rt_modocc_genes"]
ccd_at_modocc_genes = modocc["ccd_at_modocc_genes"]
ccd_n_modocc_genes = modocc["ccd_n_modocc_genes"]
nonccd_modocc_genes = modocc["nonccd_modocc_genes"]

all_temps = np.load("output/temperatures.all_temps.npy", allow_pickle=True)
t_reg_temps = np.load("output/temperatures.transcript_reg.npy", allow_pickle=True)
nont_reg_temps = np.load("output/temperatures.nontranscr_reg.npy", allow_pickle=True)
nonccd_temps = np.load("output/temperatures.nonccd_temps.npy", allow_pickle=True)

all_temp_genes = np.load("output/temperatures.all_temp_prot.npy", allow_pickle=True)
t_reg_temp_genes = np.load("output/temperatures.transcript_reg_prot.npy", allow_pickle=True)
nont_reg_temp_genes = np.load("output/temperatures.nontranscript_reg_prot.npy", allow_pickle=True)
nonccd_temps_genes = np.load("output/temperatures.nonccd_temps_prot.npy", allow_pickle=True)

#%%
def temp_modocc_plot(modoccdf, modocc_genes, modocc_genename, modocc_valname, temp_genes, temp_vals, label):
    all_modocc_medians = modocc[pd.notna(modocc_genes)].groupby(modocc_genename)[modocc_valname].median()
    all_modocc_medians_genes = all_modocc_medians.index
    all_temp_both = np.isin(temp_genes, all_modocc_medians_genes)
    all_modocc_both = np.isin(all_modocc_medians_genes, temp_genes)
    all_temp_sort_inds = np.argsort(temp_genes)[all_temp_both]
    all_modocc_genes_sort_inds = np.argsort(all_modocc_medians_genes)[all_modocc_both]
    all_temp_sort = temp_vals[all_temp_sort_inds]
    all_modocc_sort = all_modocc_medians[all_modocc_genes_sort_inds]
    plt.scatter(all_modocc_sort, all_temp_sort, label=label)

temp_modocc_plot(modocc, all_modocc_genes, "all_modocc_genes", "all_modocc", all_temp_genes, all_temps, "All Proteins")
temp_modocc_plot(modocc, ccd_at_modocc_genes, "ccd_at_modocc_genes", "ccd_at_modocc", t_reg_temp_genes, t_reg_temps, "Transcript Regulated CCD")
temp_modocc_plot(modocc, ccd_n_modocc_genes, "ccd_n_modocc_genes", "ccd_n_modocc", nont_reg_temp_genes, nont_reg_temps, "Non-Transcript Regulated CCD")
plt.xlabel("Median Modification Occupancy Per Protein")
plt.ylabel("Median Melting Temperature Per Protein")
plt.legend()
plt.show()
plt.savefig("figures/TempModocc.png")

temp_modocc_plot(modocc, all_modocc_genes, "all_modocc_genes", "all_modocc", all_temp_genes, all_temps, "All Proteins")
temp_modocc_plot(modocc, ccd_at_modocc_genes, "ccd_at_modocc_genes", "ccd_at_modocc", t_reg_temp_genes, t_reg_temps, "Transcript Regulated CCD")
temp_modocc_plot(modocc, ccd_n_modocc_genes, "ccd_n_modocc_genes", "ccd_n_modocc", nont_reg_temp_genes, nont_reg_temps, "Non-Transcript Regulated CCD")
temp_modocc_plot(modocc, nonccd_modocc_genes, "nonccd_modocc_genes", "nonccd_modocc", nonccd_temps_genes, nonccd_temps, "Non-CCD")
plt.xlabel("Median Modification Occupancy Per Protein")
plt.ylabel("Median Melting Temperature Per Protein")
plt.legend()
plt.show()
plt.savefig("figures/TempModoccWithNonccd.png")
plt.close()

#%%
