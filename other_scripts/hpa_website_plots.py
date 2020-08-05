# -*- coding: utf-8 -*-
"""
Examples (using anilin) of the code required to make the plots displayed on the HPA website.

@author: Anthony J. Cesnik
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn

class FucciCellCycle:
    '''
    Object representing the length of the FUCCI cell cycle phase transitions, which
    were manually determined by Diana M.
    '''
    def __init__(self):
        # Length of the cell cycle observed for the FUCCI cell line
        self.G1_LEN = 10.833 #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
        self.G1_S_TRANS = 2.625 #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
        self.S_G2_LEN = 11.975 #hours (this should be from the G2/M cutoff above to the end)
        self.M_LEN = 0.5 # We excluded M-phase from this analysis

        self.TOT_LEN = self.G1_LEN+self.G1_S_TRANS+self.S_G2_LEN

        self.G1_PROP = self.G1_LEN / self.TOT_LEN
        self.G1_S_PROP = self.G1_S_TRANS / self.TOT_LEN + self.G1_PROP
        self.S_G2_PROP = self.S_G2_LEN / self.TOT_LEN + self.G1_S_PROP
    
def get_fileprefixes(wp_ensg):
    '''Generate the file prefixes for given genes'''
    return np.array([f"{ensg}_{sum(wp_ensg[:ei] == ensg)}" for ei, ensg in enumerate(wp_ensg)])

def temporal_mov_avg(fucci_time, cell_intensities, mvavg_xvals, mvavg_yvals, 
                     mvperc_10p, mvperc_90p, mvperc_25p, mvperc_75p, folder, fileprefix, isprotein):
    '''
    Generates a moving average plot for one protein (isprotein = True) or transcript (isprotein = False) 
    Input: Antibody intensity measurements for the current protein
    Output: Moving average plot with scatter of protein abundance measurement for each cell over the cell cycle
    '''
    outfile = os.path.join(folder,fileprefix+'_mvavg.pdf')
    plt.figure(figsize=(5,5))
    plt.fill_between(mvavg_xvals, mvperc_10p, mvperc_90p, color="lightsteelblue" if isprotein else "bisque", label="10th & 90th Percentiles")
    plt.fill_between(mvavg_xvals, mvperc_25p, mvperc_75p, color="steelblue" if isprotein else "orange", label="25th & 75th Percentiles")
    plt.plot(mvavg_xvals, mvavg_yvals, color="blue" if isprotein else "tab:orange", label="Mean Intensity")
    plt.scatter(fucci_time, cell_intensities, c='b' if isprotein else "darkorange", alpha=0.3)
    plt.xlabel('Cell Cycle Time, hrs')
    plt.ylabel(f"{fileprefix} {'Protein' if isprotein else 'RNA'} Expression")
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.ylim(0, 1)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)
    plt.close()
    
def boxplot_result(g1, s, g2, folder, fileprefix, isprotein):
    '''Make boxplots of RNA expression by phase'''
    mmmm = np.concatenate((g1, s, g2))
    cccc = (["G1"] * len(g1))
    cccc.extend(["S"] * len(s))
    cccc.extend(["G2"] * len(g2))
    boxplot = sbn.boxplot(x=cccc, y=mmmm, color='lightsteelblue' if isprotein else 'bisque', showfliers=True)
    boxplot.set_xlabel("", size=36,fontname='Arial')
    boxplot.set_ylabel("Normalized Mean Intensity" if isprotein else "Normalized Log10(TPM)", size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    plt.title("")
    plt.savefig(f"{folder}/{fileprefix}_boxplot.pdf")
    plt.close()
    
def plot_expression_facs(cell_intensities, cell_fred, cell_fgreen, hist2d, hist2d_extent, folder, fileprefix, isprotein):
    '''
    Displaying relative expression on the FUCCI plot (log-log intensity of FUCCI markers)
    Output: a folder of FUCCI plots for each gene in a list
    '''
    plt.scatter(cell_fgreen, cell_fred, c=cell_intensities, cmap="plasma")
    plt.tight_layout()
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel("Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
    if isprotein:
        mycmap = plt.cm.gray_r
        mycmap.set_under(color='w',alpha=None)
        plt.imshow(hist2d, extent=hist2d_extent, origin="lower", cmap=mycmap)
    plt.xlabel(("Zero-Centered " if isprotein else "") + r'$Log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(("Zero-Centered " if isprotein else "") + r'$Log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(f"{folder}/{fileprefix}_fucciplot.png")
    plt.close()
    
#%% Load plotting tables
fucci = FucciCellCycle()
proteins = pd.read_csv("output/ProteinPseudotimePlotting.csv.gz", sep="\t")
rna = pd.read_csv("output/RNAPseudotimePlotting.csv.gz", sep="\t")

#%% Make an example of the protein plots
protein_idx = list(proteins["ENSG"]).index("ENSG00000011426") # anilin as an example
ensg_fileprefixes = get_fileprefixes(proteins["ENSG"]) # there are some replicates with the same ENSG, so make each file unique

# Pseudotime plot, protein
print(f"{proteins['ENSG'][protein_idx]}: Gene information")
print(f"{proteins['Compartment'][protein_idx]}: Compartment information")
print(f"{proteins['CCD'][protein_idx]}: Protein CCD information")
print(f"{ensg_fileprefixes[protein_idx]}: saving example file to this fileprefix")
temporal_mov_avg(
    [float(xx) * fucci.TOT_LEN for xx in proteins["cell_pseudotime"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["cell_intensity"][protein_idx].split(",")], 
    [float(xx) * fucci.TOT_LEN for xx in proteins["mvavg_x"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["mvavg_y"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["mvavgs_10p"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["mvavgs_90p"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["mvavgs_25p"][protein_idx].split(",")], 
    [float(xx) for xx in proteins["mvavgs_75p"][protein_idx].split(",")], 
    "figures/HpaWebsite", ensg_fileprefixes[protein_idx], isprotein=True)

# Boxplots, protein
phases = np.array(proteins["phase"][protein_idx].split(","))
g1 = phases == "G1"
s = phases == "S"
g2 = phases == "G2"
boxplot_result([float(xx) for xx in np.array(proteins["cell_intensity"][protein_idx].split(","))[g1]], 
               [float(xx) for xx in np.array(proteins["cell_intensity"][protein_idx].split(","))[s]],
               [float(xx) for xx in np.array(proteins["cell_intensity"][protein_idx].split(","))[g2]],
               "figures/HpaWebsite", ensg_fileprefixes[protein_idx], isprotein=True)

# Fucci plot, protein
h2d, xedges, yedges, image = plt.hist2d(np.concatenate([[float(xx) for xx in fgrelist.split(",")] for fgrelist in proteins["cell_fgreen"]]), 
                                        np.concatenate([[float(xx) for xx in fredlist.split(",")] for fredlist in proteins["cell_fred"]]), 
                                        bins=200)
plt.close()
plot_expression_facs([float(xx) for xx in proteins["cell_intensity"][protein_idx].split(",")],
                     [float(xx) for xx in proteins["cell_fred"][protein_idx].split(",")],
                     [float(xx) for xx in proteins["cell_fgreen"][protein_idx].split(",")],
                     h2d.T, (np.min(xedges), np.max(xedges), np.min(yedges), np.max(yedges)),
                     "figures/HpaWebsite", ensg_fileprefixes[protein_idx], True)

#%% Make example plots for RNA
rna_idx = list(rna["ENSG"]).index("ENSG00000011426") # anilin as an example

# Pseudotime plot, RNA
print(f"{rna['ENSG'][rna_idx]}: Gene information")
print(f"{rna['CCD'][rna_idx]}: RNA CCD information")
temporal_mov_avg(
    [float(xx) * fucci.TOT_LEN for xx in rna["cell_pseudotime"][rna_idx].split(",")], 
    [float(xx) for xx in rna["cell_intensity"][rna_idx].split(",")], 
    [float(xx) * fucci.TOT_LEN for xx in rna["mvavg_x"][rna_idx].split(",")], 
    [float(xx) for xx in rna["mvavg_y"][rna_idx].split(",")], 
    [float(xx) for xx in rna["mvavgs_10p"][rna_idx].split(",")], 
    [float(xx) for xx in rna["mvavgs_90p"][rna_idx].split(",")], 
    [float(xx) for xx in rna["mvavgs_25p"][rna_idx].split(",")], 
    [float(xx) for xx in rna["mvavgs_75p"][rna_idx].split(",")], 
    "figures/HpaWebsite", rna["ENSG"][rna_idx], isprotein=False)

# Boxplots, RNA
phases = np.array(rna["phase"][rna_idx].split(","))
g1 = phases == "G1"
s = phases == "S-ph"
g2 = phases == "G2M"
boxplot_result([float(xx) for xx in np.array(rna["cell_intensity"][rna_idx].split(","))[g1]], 
               [float(xx) for xx in np.array(rna["cell_intensity"][rna_idx].split(","))[s]],
               [float(xx) for xx in np.array(rna["cell_intensity"][rna_idx].split(","))[g2]],
               "figures/HpaWebsite", rna["ENSG"][rna_idx], isprotein=False)

# Fucci plot, RNA
plot_expression_facs([float(xx) for xx in rna["cell_intensity"][rna_idx].split(",")],
                     [float(xx) for xx in rna["cell_fred"][rna_idx].split(",")],
                     [float(xx) for xx in rna["cell_fgreen"][rna_idx].split(",")],
                     [],[],
                     "figures/HpaWebsite", rna["ENSG"][rna_idx], False)
