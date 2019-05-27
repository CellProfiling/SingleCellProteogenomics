#%%
from imports import *

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

# TIMING OF PHASE TRANSITIONS (MANUALLY DETERMINED BY DIANA)
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


#%% Read data into scanpy; Read phases and FACS intensities
from ac_RNASeqData import read_counts_and_phases
adata, phases_filt = read_counts_and_phases() # no qc, yet

#%% Fucci plots based on FACS intensities
colormap = { "G1" : "blue", "G2M" : "orange", "S-ph" : "green" }
legendboxes = []
labels = []
for key in colormap:
    legendboxes.append(mpatches.Rectangle((0,0), 1, 1, fc=colormap[key]))
    labels.append(key)

tt = "All"

# heatmap
phasesFiltintSeqCenter = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.Stage)]
plt.hist2d(phasesFiltintSeqCenter["Green530"], phasesFiltintSeqCenter["Red585"], bins=200)
plt.tight_layout()
plt.savefig(f"figures/FucciPlot{tt}Density.png")
plt.show()
plt.close()

# scatters
def fucci_scatter(phases_filtered, outfile):
    plt.scatter(phases_filtered["Green530"], phases_filtered["Red585"], c = phases_filtered["StageAJC"].apply(lambda x: colormap[x]))
    plt.legend(legendboxes, labels)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

phases_filtIntAjc = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.StageAJC)]
fucci_scatter(phases_filtIntAjc, f"figures/FucciPlot{tt}ByPhaseAJC.png")
fucci_scatter(phasesFiltintSeqCenter, f"figures/FucciPlot{tt}ByPhase.png")

for tt in ["355", "356", "357"]:
    phasesfilt355 = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.Stage) & phases_filt.Well_Plate.str.endswith(tt)]
    phasesfilt355ajc = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585) & pd.notnull(phases_filt.StageAJC)  & phases_filt.Well_Plate.str.endswith(tt)]
    fucci_scatter(phasesfilt355, f"figures/FucciPlot{tt}ByPhase.png")
    fucci_scatter(phasesfilt355ajc, f"figures/FucciPlot{tt}ByPhaseAJC.png")

#%% Convert FACS intensities to pseudotime
phasesFilt = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585)] # stage may be null
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

#%% Save fucci times, so they can be used in other workbooks
pd.DataFrame({"fucci_time": fucci_time}).to_csv("output/fucci_time.csv")

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


#%%
