#%% [markdown]
# # Cell Cycle scRNA-Seq Analysis
# We've shown in single cell imaging data that there is variability that is correlated to the cell cycle, as well as a majority of proteins that vary outside of the cell cycle, which might be due to metabolic processes or other sources of variation.
# 
# Here, we've collected single-cell RNA-Seq (scRNA-Seq) data for these cells.
# * How many cells were analyzed?
# * How many reads per cell?
# * How many genes show variability in expression at the RNA level?
# * How many of the genes that are temporally regulated over the cell cycle, using the fucci colors to build the cell cycle trajectory?
# * How many of the genes that show variability not correlated to the cell cycle?

#%%
# get_ipython().system('pip install scanpy pandas')
# get_ipython().system('conda install -y -c vtraag python-igraph')
# get_ipython().system('conda install -y -c conda-forge statsmodel')

# Render our plots inline
# get_ipython().run_line_magic('matplotlib', 'inline')

# some setup
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcolors
import matplotlib.patches as mpatches
import scanpy as sc
import os
import shutil
import scipy
import scipy.stats

#%% least squares, continuous time setup
from scipy.optimize import least_squares
def calc_R(xc, yc, x, y):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)
def f_2(c,x,y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    print(c)
    Ri = calc_R(c[0],c[1],x,y)
    return Ri - Ri.mean()

import decimal
from stretch_time import stretch_time

#%% TIMING OF PHASE TRANSITIONS (MANUALLY DETERMINED BY DIANA)

#hours (for the G1/S cutoff)
G1_LEN = 10.833
#hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
G1_S_TRANS = 2.625
#hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
S_G2_LEN = 11.975
#hours (this should be from the G2/M cutoff above to the end)
#M_LEN = 0.5
#We are excluding Mphase from this analysis
TOT_LEN = G1_LEN + G1_S_TRANS + S_G2_LEN
G1_PROP = G1_LEN / TOT_LEN
G1_S_PROP = G1_S_TRANS / TOT_LEN + G1_PROP
S_G2_PROP = S_G2_LEN / TOT_LEN + G1_S_PROP

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts355 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_355\\counts.tab", delimiter="\t", index_col=0)
    counts356 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_356\\counts.tab", delimiter="\t", index_col=0)
    counts357 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_357\\counts.tab", delimiter="\t", index_col=0)

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts355.columns += "_355"
    counts356.columns += "_356"
    counts357.columns += "_357"
    counts355.columns

#%%
if not os.path.isfile("AllCountsForScanpy.csv") or not os.path.isfile("355CountsForScanpy.csv"):
    counts = pd.concat([counts355,counts356,counts357], axis=1, sort=False)
    counts.T.sort_index().to_csv("AllCountsForScanpy.csv")
    counts355.T.sort_index().to_csv("355CountsForScanpy.csv")
    counts356.T.sort_index().to_csv("356CountsForScanpy.csv")
    counts357.T.sort_index().to_csv("357CountsForScanpy.csv")

#%% Read data into scanpy
dd = "All"
adata = sc.read_csv(dd + "CountsForScanpy.csv")
adata.var_names_make_unique()
adata.raw = adata

#%% Read phases and FACS intensities
phases = pd.read_csv("WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")
phases_filt = phases[phases["Well_Plate"].isin(adata.obs_names) ]
phases_filt = phases_filt.reset_index(drop=True) # remove filtered indices

#%% Assign phases and log intensities; require log intensity
adata.obs["phase"] = np.array(phases_filt["Stage"])
adata.obs["phase_ajc"] = np.array(phases_filt["StageAJC"])
adata.obs["Green530"] = np.array(phases_filt["Green530"])
adata.obs["Red585"] = np.array(phases_filt["Red585"])
adata = adata[pd.notnull(adata.obs["Green530"]) & pd.notnull(adata.obs["Red585"])] # removes dark mitotic cells

#%% Fucci plots
colormap = { "G1" : "blue", "G2M" : "orange", "S-ph" : "green" }
legendboxes = []
labels = []
for key in colormap:
    legendboxes.append(mpatches.Rectangle((0,0), 1, 1, fc=colormap[key]))
    labels.append(key)
phasesFilt = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585)] # stage may be null

tt = "All"

# heatmap
plt.hist2d(phasesFiltintSeqCenter["Green530"], phasesFiltintSeqCenter["Red585"], bins=200)
plt.tight_layout()
plt.savefig(f"figures/FucciPlot{tt}Density.png")
plt.show()

# scatters
def fucci_scatter(phases_filtered, outfile):
    plt.scatter(phases_filtered["Green530"], phases_filtered["Red585"], c = phases_filtered["StageAJC"].apply(lambda x: colormap[x]))
    plt.legend(legendboxes, labels)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.show()

phases_filtIntAjc = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.StageAJC)]
fucci_scatter(phases_filtIntAjc, f"figures/FucciPlot{tt}ByPhaseAJC.png")

phasesFiltintSeqCenter = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.Stage)]
fucci_scatter(phasesFiltintSeqCenter, f"figures/FucciPlot{tt}ByPhase.png")

for tt in ["355", "356", "357"]:
    phasesfilt355 = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.Stage) & phases_filt.Well_Plate.str.endswith(tt)]
    phasesfilt355ajc = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.StageAJC)  & phases_filt.Well_Plate.str.endswith(tt)]
    fucci_scatter(phasesfilt355, f"figures/FucciPlot{tt}ByPhase.png")
    fucci_scatter(phasesfilt355ajc, f"figures/FucciPlot{tt}ByPhaseAJC.png")

#%% Convert FACS intensities to pseudotime
x = phasesFilt["Green530"]
y = phasesFilt["Red585"]
fucci_data = np.column_stack([x, y])
center_est_xy = np.mean(x), np.mean(y)
center_est2_xy = least_squares(f_2, center_est_xy, args=(x, y))
xc_2, yc_2 = center_est2_xy.x
Ri_2       = calc_R(*center_est2_xy.x,x,y)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)

# Center data
centered_data = fucci_data - center_est2_xy.x

# Convert data to polar
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

pol_data = cart2pol(centered_data[:,0],centered_data[:,1])
pol_sort_inds = np.argsort(pol_data[1])
pol_sort_rho = pol_data[0][pol_sort_inds]
pol_sort_phi = pol_data[1][pol_sort_inds]
centered_data_sort0 = centered_data[pol_sort_inds,0]
centered_data_sort1 = centered_data[pol_sort_inds,1]

# Rezero to minimum --resoning, cells disappear during mitosis, so we should have the fewest detected cells there
NBINS = 150 #number of bins, arbitrary choice for now
bins = plt.hist(pol_sort_phi,NBINS)
start_phi = bins[1][np.argmin(bins[0])]

# Move those points to the other side
more_than_start = np.greater(pol_sort_phi,start_phi)
less_than_start = np.less_equal(pol_sort_phi,start_phi)
pol_sort_rho_reorder = np.concatenate((pol_sort_rho[more_than_start],pol_sort_rho[less_than_start]))
pol_sort_inds_reorder = np.concatenate((pol_sort_inds[more_than_start],pol_sort_inds[less_than_start]))
pol_sort_phi_reorder = np.concatenate((pol_sort_phi[more_than_start],pol_sort_phi[less_than_start]+np.pi*2))
pol_sort_centered_data0 = np.concatenate((centered_data_sort0[more_than_start],centered_data_sort0[less_than_start]))
pol_sort_centered_data1 = np.concatenate((centered_data_sort1[more_than_start],centered_data_sort1[less_than_start]))
pol_sort_shift = pol_sort_phi_reorder+np.abs(np.min(pol_sort_phi_reorder))

# Shift and re-scale "time"
# reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
pol_sort_norm = pol_sort_shift/np.max(pol_sort_shift)
pol_sort_norm_rev = 1 - pol_sort_norm 
pol_sort_norm_rev = stretch_time(pol_sort_norm_rev)
plt.tight_layout()
plt.savefig(f"figures/Fucci{dd}PseudotimeHist.png")
plt.show()

# Apply uniform radius (rho) and convert back
cart_data_ur = pol2cart(np.repeat(R_2, len(centered_data)), pol_data[1])

#%% Assign cells a pseudotime and visualize in fucci plot
pol_unsort = np.argsort(pol_sort_inds_reorder)
fucci_time = pol_sort_norm_rev[pol_unsort]
adata.obs["fucci_time"] = fucci_time
phasesFilt["fucci_time"] = fucci_time

plt.scatter(phasesFilt["Green530"], phasesFilt["Red585"], c = phasesFilt["fucci_time"])
plt.tight_layout()
plt.colorbar()
plt.savefig(f"figures/Fucci{dd}FucciPseudotime.png")
plt.show()

#%% Visualize that pseudotime result
start_pt = pol2cart(R_2,start_phi)
g1_end_pt = pol2cart(R_2,start_phi + (1 - G1_PROP) * 2 * np.pi)
g1s_end_pt = pol2cart(R_2,start_phi + (1 - G1_S_PROP) * 2 * np.pi)

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

#%% QC and filtering
sc.pl.highest_expr_genes(adata, n_top=20, show=True, save=True)
shutil.move("figures/highest_expr_genes.pdf", f"figures/highest_expr_genes_{dd}Cells.pdf")

sc.pp.filter_cells(adata, min_genes=400)
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

#%% Post filtering QC
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=True, save=True)
shutil.move("figures/filter_genes_dispersion.pdf", f"figures/filter_genes_dispersion{dd}Cells.pdf")

#%% UMAP plots

# UMAP statistics first
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.umap(adata, color=["phase"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqCenterPhase.pdf")
sc.pl.umap(adata, color = ["phase_ajc"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhase.pdf")

#%% Read in the published CCD genes
ccd_regev=pd.read_csv("ccd_regev.txt")
ccd_regev_filtered = [gene for gene in ccd_regev["gene"] if gene in adata.var_names]
adata_ccdregev = adata[:, ccd_regev_filtered]
adata_ccdregev.var_names_make_unique()

#%% UMAP with just the Regev cell-cycle dependent genes (CCD) 
# https://science.sciencemag.org/content/352/6282/189
sc.pp.neighbors(adata_ccdregev, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_ccdregev)
sc.pl.umap(adata_ccdregev, color="phase_ajc", show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsAjcPhaseCcdRegev.pdf")

#%% Heatmaps with published CCD genes again
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, ccd_regev_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdRegev.pdf")
# sc.pl.heatmap(adata, ccd_regev_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdRegev.pdf")

#%%
ccd=pd.read_csv("ccd_genes.txt")
nonccd=pd.read_csv("nonccd_genes.txt")
ccd_filtered = [gene for gene in ccd["gene"] if gene in adata.var_names]
nonccd_filtered = [gene for gene in nonccd["gene"] if gene in adata.var_names]

#%% Heatmaps with Diana's first few CCD genes
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, ccd_filtered[:50], "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectCcdDiana.pdf")
# sc.pl.heatmap(adata, ccd_filtered[:50], "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectCcdDiana.pdf")

# # Heatmaps with Diana's CCD genes
# sc.pl.heatmap(adata, ccd_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseCcdDiana.pdf")
# sc.pl.heatmap(adata, ccd_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseCcdDiana.pdf")

#%% Heatmaps with Diana's first few Non-CCD genes
# kind of broken; need to set the colormap back to 'reds' somehow
# sc.pl.heatmap(adata, nonccd_filtered[:50], "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseSelectNonCcdDiana.pdf")
# sc.pl.heatmap(adata, nonccd_filtered[:50], "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseSelectNonCcdDiana.pdf")

# # Heatmaps with Diana's Non-CCD genes
# sc.pl.heatmap(adata, nonccd_filtered, "phase", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsSeqCenterPhaseNonCcdDiana.pdf")
# sc.pl.heatmap(adata, nonccd_filtered, "phase_ajc", show=True, save=True)
# shutil.move("figures/heatmap.pdf", f"figures/heatmap{dd}CellsAjcPhaseNonCcdDiana.pdf")

#%% Louvian clustering -- this doesn't work; they use C++ packages that don't install easily.
# I checked and I can use Seurat for this in R
# sc.tl.louvain(adata)

#%% partition-based graph abstraction (PAGA), and this doesn't work because of louvain either
# sc.tl.paga(adata)
# sc.pl.paga_compare(adata, edges=True, threshold=0.05)

#%% pseudotime reconstruction
adata.uns["iroot"] = 0
sc.tl.dpt(adata, n_branchings=0)
sc.pl.umap(adata, color='dpt_pseudotime', show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsPredictedPseudotime.pdf")
sc.pl.diffmap(adata, color='dpt_pseudotime', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsPredictedPseudotime3d.pdf")

#%% fucci pseudotime
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.diffmap(adata, color='fucci_time', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsFucciPseudotime3d.pdf")
sc.pl.umap(adata, color=["fucci_time"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqFucciPseudotime.pdf")

#%% Expression vs Pseudotime, uncomment to run again
def plot_expression_pseudotime(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        plt.scatter(adata.obs["fucci_time"], normalized_exp_data)
        plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
        plt.ylabel("RNA-Seq Counts, Normalized By Cell",size=36,fontname='Arial')
        plt.title(gene,size=36,fontname='Arial')
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

# plot_expression_pseudotime(ccd_regev_filtered, "figures/RegevGeneProfiles")
# plot_expression_pseudotime(ccd_filtered, "figures/DianaCcdGeneProfiles")
# plot_expression_pseudotime(nonccd_filtered, "figures/DianaNonCcdGeneProfiles")
plot_expression_pseudotime(["ENO1"], "figures/OtherGeneProfiles")

#%% Cell cycle regulated ANOVA
expression_data = np.exp(adata.X) - 1
normalized_exp_data = expression_data / np.max(expression_data)
stages = np.array(adata.obs["phase"])
g1_exp = np.take(normalized_exp_data, np.nonzero(stages == "G1")[0], axis=0)
s_exp = np.take(normalized_exp_data, np.nonzero(stages == "S-ph")[0], axis=0)
g2_exp = np.take(normalized_exp_data, np.nonzero(stages == "G2M")[0], axis=0)
tests_fp = [stats.f_oneway(g1_exp[:,geneidx], s_exp[:,geneidx], g2_exp[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
pvals = [p for (F, p) in tests_fp]

# benjimini-hochberg multiple testing correction
# source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
alpha = 0.01
def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection'''
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)

pvals_sortind = np.argsort(pvals)
pvals_sorted = np.take(pvals, pvals_sortind)
ecdffactor = _ecdf(pvals_sorted)
reject = pvals_sorted <= ecdffactor*alpha
reject = pvals_sorted <= ecdffactor*alpha
if reject.any():
    rejectmax = max(np.nonzero(reject)[0])
    reject[:rejectmax] = True
pvals_corrected_raw = pvals_sorted / ecdffactor
pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
del pvals_corrected_raw
pvals_corrected[pvals_corrected>1] = 1
pvals_corrected_ = np.empty_like(pvals_corrected)

# deal with sorting
pvals_corrected_[pvals_sortind] = pvals_corrected
del pvals_corrected
reject_ = np.empty_like(reject)
reject_[pvals_sortind] = reject

# BenjiHoch is actually pretty liberal for this dataset. What about bonferroni?
# bonferroni MTC
alphaBonf = alpha / float(len(pvals))
rejectBonf = pvals_sorted <= alphaBonf
pvals_correctedBonf = pvals_sorted * float(len(pvals))
pvals_correctedBonf_unsorted = np.empty_like(pvals_correctedBonf) 
pvals_correctedBonf_unsorted[pvals_sortind] = pvals_correctedBonf
rejectBonf_unsorted = np.empty_like(rejectBonf)
rejectBonf_unsorted[pvals_sortind] = rejectBonf

anova_tests = pd.DataFrame(
    {"gene" : adata.var_names, 
    "pvalue" : pvals, 
    "pvaladj_BH" : pvals_corrected_, 
    "reject_BH" : reject_, 
    "pvaladj_B" : pvals_correctedBonf_unsorted, 
    "reject_B" : rejectBonf_unsorted})
anova_tests.to_csv("output/transcript_regulation.csv")


#%% Plotting variances of gene expression
means = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
variances = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:]))]
mean_regev = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_regev_filtered]
variances_regev = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_regev_filtered]
mean_dianaccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
variances_dianaccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
mean_diananonccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]
variances_diananonccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x = means, y = variances, label="All Genes")
ax1.scatter(x = mean_regev, y = variances_regev, label="Regev CCD Genes")
ax1.scatter(x = mean_dianaccd, y = variances_dianaccd, label="Fucci CCD Genes")
ax1.scatter(x = mean_diananonccd, y = variances_diananonccd, label="Fucci Non-CCD Genes")
plt.legend(loc="upper left")
plt.xlabel("Mean Expression")
plt.ylabel("Stdev Expression")
plt.tight_layout()
plt.savefig("figures/stdev_expression.png")
plt.show()
plt.close()

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))

fig = plt.figure()
ax1 = fig.add_subplot(111)
bins=np.histogram(np.hstack((variances, variances_regev, variances_dianaccd, variances_diananonccd)), bins=40)[1] #get the bin edges
ax1.hist(variances_regev, bins=bins, weights=weights(variances_regev), 
    label="Regev CCD Genes")
ax1.hist(variances_dianaccd, bins=bins, weights=weights(variances_dianaccd),
    label="Fucci CCD Genes")
ax1.hist(variances_diananonccd, bins=bins, weights=weights(variances_diananonccd), 
    label="Fucci Non-CCD Genes")
ax1.hist(variances, bins=bins, weights=weights(variances), 
    label="All Genes")
plt.legend(loc="upper right")
plt.xlabel("Stdev Expression")
plt.ylabel("Count, Normalized to 1")
plt.tight_layout()
plt.savefig("figures/stdev_expression_hist.png")
plt.show()
plt.close()

#%% Expression boxplots
def format_p(p):
    '''3 decimal places, scientific notation'''
    return '{:0.3e}'.format(p)

def plot_expression_boxplots(genelist, outanova, outfolder):
    notfound = []
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        if gene not in adata.var_names:
            notfound.append(gene)
            continue
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        expression_stage = pd.DataFrame({gene : normalized_exp_data, "stage" : adata.obs["phase"]})
        exp_stage_filt = expression_stage[expression_stage.stage != "nan"].reset_index(drop=True)
        boxplot = exp_stage_filt.boxplot(gene, by="stage", figsize=(12, 8))
        boxplot.set_xlabel("Cell Cycle Stage", size=36,fontname='Arial')
        boxplot.set_ylabel("Normalized Transcript Expression", size=36,fontname='Arial')
        boxplot.tick_params(axis="both", which="major", labelsize=16)
        qqq = anova_tests[anova_tests.gene == gene].iloc[0]["pvaladj"]
        rej = anova_tests[anova_tests.gene == gene].iloc[0]["reject"]
        boxplot.set_title(f"{gene} (Q_bh = {format_p(qqq)})",size=36,fontname='Arial')
        boxplot.get_figure().savefig(f"{outfolder}/{gene}_boxplot.png")
        plt.close()
    pd.DataFrame({"gene" : notfound}).to_csv(f"{outfolder}/GenesFilteredInScRNASeqAnalysis.csv")

plot_expression_boxplots(ccd_regev_filtered, "ccd_regev_filtered", "figures/RegevGeneBoxplots")
plot_expression_boxplots(ccd_filtered, "ccd_filtered", "figures/DianaCcdGeneBoxplots")
plot_expression_boxplots(nonccd_filtered, "nonccd_filtered", "figures/DianaNonCcdGeneBoxplots")

#%% [markdown]
## Summary of ANOVA analysis
#
# I looked for transcript regulation across the cell cycle  on all 17,199 genes that made it through filtering 
# (using ANOVA with Benjimini-Hochberg multiple testing correction).  
# Of those, about a third have evidence of cell cycle regulation at the transcript level (at 1% FDR). 
# Of the gene set from the Regev study (single cell transcriptomics), all but one (99%) of the genes 
# have evidence of transcript regulation over the cell cycle. An example of one of the surprises is above, 
# where geminin is significant (rejecting the hypothesis that all stages have the same expression level), 
# even while this effect is subtle in the pseudotime analysis. It also seems to make sense that 
# G1 and S phases have higher expression, while G2M has lower expression of geminin.
#
# In Diana's 464 cell-cycle dependent genes, 197 had evidence of cell-cycle regulation at the 
# transcript level (48% of the 410 genes that weren't filtered in scRNA-Seq preprocessing). 
# Strangely, of the 890 non cell-cycle dependent genes, including 811 that made it through 
# scRNA-Seq preprocessing, 362 genes (45%) showed evidence for cell cycle dependence at the 
# transcript level.
#
# I tried a Bonferroni correction (alpha=0.01), which is a bit more conservative, 
# and 97% of the Regev genes, 27% of Diana's CCD proteins, and 20% of Diana's 
# non-CCD proteins showed cell cycle regulation for transcripts.

#%% Expression Fucci, uncomment to run again
def plot_expression_facs(genelist, pppp, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        plt.scatter(pppp["Green530"], pppp["Red585"], normalized_exp_data)
        plt.tight_layout()
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel("Log Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
        plt.title(gene,size=20,fontname='Arial')
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

phasesfiltfilt = phasesFilt[phasesFilt["Well_Plate"].isin(adata.obs_names)]

# plot_expression_facs(ccd_regev_filtered, phasesfiltfilt, "figures/RegevGeneFucci")
# plot_expression_facs(ccd_filtered, phasesfiltfilt, "figures/DianaCcdGeneFucci")
# plot_expression_facs(nonccd_filtered, phasesfiltfilt, "figures/DianaNonCcdGeneFucci")

#%% Expression-FUCCI facs of anillin by plate
def plot_expression_facs_plate(genelist, pppp, plate, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    plt.subplot(1, 3)
    for tt in ["355", "356", "357"]:
        for gene in genelist:
            pppp[gene] = adata.X[:,list(adata.var_names).index(gene)]
            pppplate = pppp[pd.notnull(pppp.Green530) & pd.notnull(pppp.Red585) & pd.notnull(pppp.Stage) & phasesfiltfilt.Well_Plate.str.endswith(plate)]
            plt.scatter(pppplate["Green530"], pppplate["Red585"], c = pppplate[gene])
            plt.tight_layout()
            cbar = plt.colorbar()
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel("Log Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
            plt.title(gene,size=20,fontname='Arial')
            plt.savefig(f"{outfolder}/{gene}Plate{plate}.png")
            plt.close()
            pppp.drop(gene, 1)

    plot_expression_facs_plate(
        ["ANLN"], phasesfiltfilt, tt, f"figures/FucciPlotByPlates")

#%% Expression UMAP, uncomment to run again
def plot_expression_umap(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        adata.obs[gene] = normalized_exp_data
        sc.pl.umap(adata, color=gene, show=False, save=True)
        shutil.move("figures/umap.pdf", f"{outfolder}/{gene}.pdf")
        plt.close()
        adata.obs.drop(gene, 1)

# plot_expression_umap(ccd_regev_filtered, "figures/RegevGeneUmap")
# plot_expression_umap(ccd_filtered, "figures/DianaCcdGeneUmap")
# plot_expression_umap(nonccd_filtered, "figures/DianaNonCcdGeneUmap")


#%%
