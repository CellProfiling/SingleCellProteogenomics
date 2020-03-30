#%% imports
from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, stretch_time, FucciCellCycle, FucciPseudotime
from scipy.optimize import least_squares
import decimal
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

#%% Fucci plots based on FACS intensities
# Idea: make log-log fucci intensity plots for the cells analyzed by RNA-Seq
# Execution: matplotlib
# Output: scatters
tt = "All"
adata, phases_filt = read_counts_and_phases(tt, "Counts", False, "protein_coding") # no qc, yet
colormap = { "G1" : "blue", "G2M" : "orange", "S-ph" : "green" }
legendboxes = []
labels = []
for key in colormap:
    legendboxes.append(mpatches.Rectangle((0,0), 1, 1, fc=colormap[key]))
    labels.append(key)

# heatmap
phasesFiltintSeqCenter = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.Stage)]
plt.hist2d(phasesFiltintSeqCenter["Green530"], phasesFiltintSeqCenter["Red585"], bins=200)
plt.tight_layout()
plt.savefig(f"figures/FucciPlot{tt}Density.png")
plt.show()
plt.close()

# scatters
def fucci_scatter(phases_filtered, outfile):
    '''Generate a FUCCI plot with log intensities of the GFP and RFP tags'''
    plt.scatter(phases_filtered["Green530"], phases_filtered["Red585"], c = phases_filtered["Stage"].apply(lambda x: colormap[x]))
    plt.legend(legendboxes, labels)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

fucci_scatter(phasesFiltintSeqCenter, f"figures/FucciPlot{tt}ByPhase.png")

#%% Convert FACS intensities to pseudotime
# Idea: Use the polar coordinate pseudotime calculations to calculate the pseudotime for each cell
# Execution: Adapt Devin's code for the cells sorted for RNA-Seq
# Output: Plot of all fucci pseudotimes; table of pseudotimes for each cell


phasesFilt = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585)] # stage may be null
polar_coord_results = FucciPseudotime.fucci_polar_coords(phasesFilt["Green530"], phasesFilt["Red585"], "RNA")
pol_sort_norm_rev, centered_data, pol_sort_centered_data0, pol_sort_centered_data1, pol_sort_inds, pol_sort_inds_reorder, more_than_start, less_than_start, cart_data_ur, start_point = polar_coord_results

# Assign cells a pseudotime and visualize in fucci plot
pol_unsort = np.argsort(pol_sort_inds_reorder)
fucci_time = pol_sort_norm_rev[pol_unsort]
adata.obs["fucci_time"] = fucci_time
phasesFilt["fucci_time"] = fucci_time

plt.figure(figsize=(6,5))
plt.scatter(phasesFilt["Green530"], phasesFilt["Red585"], c = phasesFilt["fucci_time"], cmap="RdYlGn")
cbar = plt.colorbar()
cbar.set_label('Pseudotime',fontname='Arial',size=20)
cbar.ax.tick_params(labelsize=18)
plt.xlabel("log10(GMNN GFP Intensity)",fontname='Arial',size=20)
plt.ylabel("log10(CDT1 RFP Intensity)",fontname='Arial',size=20)
plt.tight_layout()
plt.savefig(f"figures/FucciAllFucciPseudotime.pdf")
plt.show()
plt.close()

# Save fucci times, so they can be used in other workbooks
pd.DataFrame({"fucci_time": fucci_time}).to_csv("output/fucci_time.csv")

#%% Visualize that pseudotime result
# Idea: Generate a plot of the centered data
# Execution: hist2d
# Output: 2d hist

def plot_annotate_time(fraction):
    pt = pol2cart(R_2,start_phi + (1 - fraction) * 2 * np.pi)
    plt.scatter(pt[0],pt[1],c='c',linewidths=4)
    plt.annotate(f"  {round(fraction * TOT_LEN, 2)} hrs", (pt[0], pt[1]))

def drange(x, y, jump):
  while x < y:
    yield float(x)
    x += decimal.Decimal(jump)

def fucci_hist2d(centered_data, cart_data_ur, start_pt, outfolder, nbins=200):
    fig, ax1 = plt.subplots(figsize=(10,10))
    mycmap = plt.cm.gray_r
    mycmap.set_under(color='w',alpha=None)
    ax1.hist2d(centered_data[:,0],centered_data[:,1],bins=nbins,alpha=1,cmap=mycmap)
    hist, xbins, ybins = np.histogram2d(cart_data_ur[0], cart_data_ur[1], bins=nbins, normed=True)
    extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
    im = ax1.imshow(
            np.ma.masked_where(hist == 0, hist).T,
            interpolation='nearest',
            origin='lower',
            extent=extent,
            cmap='plasma')
    plt.scatter(start_pt[0],start_pt[1],c='c',linewidths=4)
    plt.scatter(g1_end_pt[0],g1_end_pt[1],c='c',linewidths=4)
    plt.scatter(g1s_end_pt[0],g1s_end_pt[1],c='c',linewidths=4)
    plt.scatter(0,0,c='m',linewidths=4)
    plt.annotate(f"  0 hrs (start)", (start_pt[0],start_pt[1]))
    plt.annotate(f"  {G1_LEN} hrs (end of G1)", (g1_end_pt[0],g1_end_pt[1]))
    plt.annotate(f"  {G1_LEN + G1_S_TRANS} hrs (end of S)", (g1s_end_pt[0],g1s_end_pt[1]))

    for yeah in list(drange(decimal.Decimal(0.1), 0.9, '0.1')):
        plot_annotate_time(yeah)

    plt.xlabel(r'$\propto log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(r'$\propto log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'masked_polar_hist.pdf'),transparent=True)

fucci_hist2d(centered_data, cart_data_ur, start_pt, "figures", NBINS)


#%% Make boxplots
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists
def boxplot_result(g1, s, g2, outfolder, ensg):
    if not os.path.exists(f"{outfolder}_png"): os.mkdir(f"{outfolder}_png")
    if not os.path.exists(f"{outfolder}_pdf"): os.mkdir(f"{outfolder}_pdf")
    mmmm = np.concatenate((g1, s, g2))
    cccc = (["G1"] * len(g1))
    cccc.extend(["G1/S"] * len(s))
    cccc.extend(["G2"] * len(g2))
    boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=True)
    boxplot.set_xlabel("", size=36,fontname='Arial')
    boxplot.set_ylabel("Normalized Log10(TPM)", size=18,fontname='Arial')
    boxplot.tick_params(axis="both", which="major", labelsize=14)
    plt.title("")
    plt.savefig(f"{outfolder}_png/GaussianClusteringProtein_{ensg}.png")
    plt.savefig(f"{outfolder}_pdf/GaussianClusteringProtein_{ensg}.pdf")
    plt.close()
    
dd = "All"
biotype_to_use="protein_coding"
count_or_rpkm = "Tpms"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, False, biotype_to_use)
adata, phasesfilt = qc_filtering(adata, do_log_normalize= True, do_remove_blob=True)

g1, s, g2 = adata.obs["phase"] == "G1", adata.obs["phase"] == "S-ph", adata.obs["phase"] == "G2M"
# for iii, ensg in enumerate(adata.var_names):
#     maxtpm = np.max(np.concatenate((adata.X[g1,iii], adata.X[s,iii], adata.X[g2,iii])))
#     boxplot_result(adata.X[g1,iii] / maxtpm, adata.X[s,iii] / maxtpm, adata.X[g2,iii] / maxtpm, "figures/RNABoxplotByPhase", ensg)
