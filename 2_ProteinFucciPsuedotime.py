#%% Imports
from imports import *
import numpy as np
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
import collections
import fucci_plotting
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
plate = np.load("output/plate.filterNegStain.npy", allow_pickle=True)
u_plate = np.load("output/u_plate.filterNegStain.npy", allow_pickle=True)
well_plate=np.load("output/well_plate.filterNegStain.npy", allow_pickle=True)
imgnb=np.load("output/imgnb.filterNegStain.npy", allow_pickle=True)
u_well_plates=np.load("output/u_well_plates.filterNegStain.npy", allow_pickle=True)
ab_objnum=np.load("output/ab_objnum.filterNegStain.npy", allow_pickle=True)
ab_nuc=np.load("output/ab_nuc.filterNegStain.npy", allow_pickle=True)
ab_cyto=np.load("output/ab_cyto.filterNegStain.npy", allow_pickle=True)
ab_cell=np.load("output/ab_cell.filterNegStain.npy", allow_pickle=True)
mt_cell=np.load("output/mt_cell.filterNegStain.npy", allow_pickle=True)
green_fucci=np.load("output/green_fucci.filterNegStain.npy", allow_pickle=True)
red_fucci=np.load("output/red_fucci.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc=log_green_fucci_zeroc=np.load("output/log_green_fucci_zeroc.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc=np.load("output/log_red_fucci_zeroc.filterNegStain.npy", allow_pickle=True)
log_green_fucci_zeroc_rescale=np.load("output/log_green_fucci_zeroc_rescale.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc_rescale=np.load("output/log_red_fucci_zeroc_rescale.filterNegStain.npy", allow_pickle=True)
wp_cell_kruskal_adj = np.load("output/wp_cell_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_nuc_kruskal_adj = np.load("output/wp_nuc_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_cyto_kruskal_adj = np.load("output/wp_cyto_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_cell = np.load("output/wp_pass_bh_cell.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_nuc = np.load("output/wp_pass_bh_nuc.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_cyto = np.load("output/wp_pass_bh_cyto.filterNegStain.npy", allow_pickle=True)
fucci_data = np.load("output/fucci_data.filterNegStain.npy", allow_pickle=True)
print("loaded")

#%% 
# Idea: Calculate the polar coordinates and other stuff
# Exec: Devin's calculations
# Output: fucci plot with polar coordinates

NBINS = 150 #number of bins, arbitrary choice for now

def calc_R(xc, yc, x, y):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c,x,y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    print(c)
    Ri = calc_R(c[0],c[1],x,y)
    return Ri - Ri.mean()

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

# Center data
x0 = np.ones(5)
x = fucci_data[:,0]
y = fucci_data[:,1]
center_estimate = np.mean(fucci_data[:,0]), np.mean(fucci_data[:,1])
center_2 = least_squares(f_2, center_estimate, args=(x, y))
xc_2, yc_2 = center_2.x
Ri_2       = calc_R(*center_2.x,x,y)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)
centered_data = fucci_data-center_2.x

# Convert data to polar
pol_data = cart2pol(centered_data[:,0],centered_data[:,1])
pol_sort_inds = np.argsort(pol_data[1])
pol_sort_rho = pol_data[0][pol_sort_inds]
pol_sort_phi = pol_data[1][pol_sort_inds]
centered_data_sort0 = centered_data[pol_sort_inds,0]
centered_data_sort1 = centered_data[pol_sort_inds,1]

# Sort data by polar coordinates
def pol_sort(inds, nuc, cyto, cell, mt):
    return nuc[inds], cyto[inds], cell[inds], mt[inds]
well_plate_sort = well_plate[pol_sort_inds]
ab_nuc_sort, ab_cyto_sort, ab_cell_sort, mt_cell_sort = pol_sort(pol_sort_inds,ab_nuc,ab_cyto,ab_cell,mt_cell)
# ab_nuc_sort_int, ab_cyto_sort_int, ab_cell_sort_int, mt_cell_sort_int = pol_sort(pol_sort_inds,ab_nuc_int,ab_cyto_int,ab_cell_int,mt_cell_int)
fred_sort = red_fucci[pol_sort_inds]
fgreen_sort = green_fucci[pol_sort_inds]

# Rezero to minimum --reasoning, cells disappear during mitosis, so we should have the fewest detected cells there
bins = plt.hist(pol_sort_phi,NBINS)
start_phi = bins[1][np.argmin(bins[0])]

# Move those points to the other side
more_than_start = np.greater(pol_sort_phi,start_phi)
less_than_start = np.less_equal(pol_sort_phi,start_phi)
def pol_reord(arr):
    return np.concatenate((arr[more_than_start],arr[less_than_start]))

pol_sort_well_plate = pol_reord(well_plate_sort)
# gene, antibody, Uniprot, ENSG
pol_sort_rho_reorder = pol_reord(pol_sort_rho)
pol_sort_inds_reorder = pol_reord(pol_sort_inds)
pol_sort_phi_reorder = np.concatenate((pol_sort_phi[more_than_start],pol_sort_phi[less_than_start]+np.pi*2))
pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell = pol_reord(ab_nuc_sort), pol_reord(ab_cyto_sort), pol_reord(ab_cell_sort), pol_reord(mt_cell_sort)
# pol_sort_ab_nuc_int, pol_sort_ab_cyto_int, pol_sort_ab_cell_int, pol_sort_mt_cell_int = pol_reord(ab_nuc_sort_int), pol_reord(ab_cyto_sort_int), pol_reord(ab_cell_sort_int), pol_reord(mt_cell_sort_int)
pol_sort_centered_data0, pol_sort_centered_data1 = pol_reord(centered_data_sort0), pol_reord(centered_data_sort1)
pol_sort_fred = pol_reord(fred_sort)
pol_sort_fgreen = pol_reord(fgreen_sort)

#shift and re-scale "time"; reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
pol_sort_shift = pol_sort_phi_reorder+np.abs(np.min(pol_sort_phi_reorder))
pol_sort_norm = pol_sort_shift/np.max(pol_sort_shift)
pol_sort_norm_rev = 1-pol_sort_norm
pol_sort_norm_rev = stretch_time(pol_sort_norm_rev)

#apply uniform radius (rho) and convert back
cart_data_ur = pol2cart(np.repeat(R_2,len(centered_data)), pol_data[1])

def fucci_hist2d(centered_data,cart_data_ur,start_pt,nbins=200):
    fig, ax1 = plt.subplots(figsize=(10,10))
    mycmap = plt.cm.gray_r
    mycmap.set_under(color='w',alpha=None)
    ax1.hist2d(centered_data[:,0],centered_data[:,1],bins=nbins,alpha=1,cmap=mycmap)
    hist, xbins, ybins = np.histogram2d(cart_data_ur[0],cart_data_ur[1], bins=nbins, normed=True)
    extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
    im = ax1.imshow(
            np.ma.masked_where(hist == 0, hist).T,
            interpolation='nearest',
            origin='lower',
            extent=extent,
            cmap='plasma')
    plt.scatter(start_pt[0],start_pt[1],c='c',linewidths=4)
    plt.scatter(0,0,c='m',linewidths=4)
    plt.xlabel(r'$\propto log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(r'$\propto log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join('figures/masked_polar_hist.pdf'),transparent=True)
    plt.close()

#visualize that result
start_pt = pol2cart(R_2,start_phi)
fucci_hist2d(centered_data,cart_data_ur,start_pt)

#%%
# Idea: process the well data
# Exec: use Devin's code
# Output: the moving average plots for each gene studied
PSIN_INIT = [np.nan,1,1,np.nan]
PSIN_BOUNDS = ((0, 1/6, 1/2, 0), (1, 100, 100, 1))
OUTLIER_NAMES = ['KIAA2026_HPA002109',
                'FEN1_HPA006581',
                'FASN_HPA006461',
                'EMP2_HPA014711',
                'CYTIP_HPA007191',
                'KREMEN2_HPA003223',
                'NXNL2_HPA045526',
                'CTU2_HPA041894',
                'GMNN_HPA054597']
PERCVAR_CUT = 0.1 #Min % of the variance explained by the cell cycle.
FDR_CUT = 0.05 #False discovery rate we will tolerate
WINDOW = 20 #Number of points for moving average window, arbitrary choice
DO_PLOTS = True #flag of whether to plot each well and save the plot
TPLOT_MODE = 'avg' #can have values of: 'avg', 'psin'
HIGHLIGHTS = ['ORC6','DUSP19','BUB1B','DPH2', 'FLI1']

def fun(p, x, y):
    #return x[0] * np.exp(-x[1] * t) * np.sin(x[2] * t) - y
    return np.abs(x/p[0])**p[1] + np.abs(y/p[2])**p[3] * (1+p[4]*x)/(1-p[4]*x)

def psin_fit(p,x,y):
    return p[0]*np.sin(np.pi*x**p[1])**p[2]+p[3]-y

def psin_eq(x,b,a,g,c):
    return b*np.sin(np.pi*x**a)**g+c

def psin_predict(p,x):
    return p[0]*np.sin(np.pi*x**p[1])**p[2]+p[3]

def get_val_compartment(var_vec,curr_compartment):
    #var vec is ordered [cell,nuc,cyto]
    if curr_compartment=='cell':
        return var_vec[0]
    elif curr_compartment=='nucleus':
        return var_vec[1]
    elif curr_compartment=='cytosol':
        return var_vec[2]
    else:
        return -1

def mvavg_perc_var(yvals,mv_window):
    yval_avg = np.convolve(yvals,np.ones((mv_window,))/mv_window, mode='valid')
    return np.var(yval_avg)/np.var(yvals),yval_avg

def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    #Written by: Olivia Guest github.com/oliviaguest/gini/blob/master/gini.py
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))

def fix_nans(binned_values):
    for i,val in enumerate(binned_values):
        if np.isnan(val):
            if i==(len(binned_values)-1):
                prevval = i-1
                nextval = i+1
            else:
                prevval = binned_values[i-1]
                nextval = binned_values[i+1]
            prevval = np.nan
            jind = i
            count = 1
            while np.isnan(prevval) and count<len(binned_values):
                if jind>0:
                    jind = jind-1
                    prevval = binned_values[jind]
                elif jind==0:
                    jind = len(binned_values)-1
                    prevval = binned_values[jind]
                count = count+1
            nextval = np.nan
            jind = i
            count = 1
            while np.isnan(nextval) and count<len(binned_values):
                if jind<(len(binned_values)-1):
                    jind = jind+1
                    nextval = binned_values[jind]
                elif jind==(len(binned_values)-1):
                    print('doing the end of list')
                    jind = 0
                    nextval = binned_values[jind]
                count = count+1
            binned_values[i] = np.mean([prevval,nextval])
    return binned_values

def temporal_mov_avg(curr_pol,curr_ab_norm,x_fit,y_fit,outname,outsuff):
    plt.close()
    outfile = os.path.join(outname,outsuff+'_mvavg.pdf')
    if os.path.exists(outfile): return
    #plot data
    bin_size = 25
    df = pd.DataFrame({"time" : curr_pol, "intensity" : curr_ab_norm})
    # plt.scatter(curr_pol,curr_ab_norm,c='c')
    plt.figure(figsize=(5,5))
    plt.plot(df["time"],df["intensity"].rolling(bin_size).mean(),color="blue")
    plt.fill_between(df["time"], 
        df["intensity"].rolling(bin_size).quantile(0.10),
        df["intensity"].rolling(bin_size).quantile(0.90),
        color="lightsteelblue",
        label="10th & 90th Percentiles")
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xlabel('Pseudotime')
    plt.ylabel(outsuff.split('_')[0] + ' Protein Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
    # plt.legend(fontsize=14)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)

class experiment:
    def __init__(self, outsuff, curr_compartment,
               perc_var_cell, perc_var_nuc, perc_var_cyto,
               max_val_cell, max_val_nuc, max_val_cyto,
               curr_pval, var_compartment, gini_compartment, curr_fdr,
               perc_var_fred, perc_var_fgreen):
        self.outsuff = outsuff
        self.compartment = curr_compartment
        self.perc_var_cell =  perc_var_cell
        self.perc_var_nuc = perc_var_nuc
        self.perc_var_cyto = perc_var_cyto
        self.max_val_cell = max_val_cell
        self.max_val_nuc = max_val_nuc
        self.max_val_cyto = max_val_cyto
        self.pval = curr_pval
        self.var_compartment = var_compartment
        self.gini_compartment = gini_compartment
        self.curr_fdr = curr_fdr
        self.perc_var_fred = perc_var_fred
        self.perc_var_fgreen = perc_var_fgreen
        self.perc_var_compartment = get_val_compartment(
            [perc_var_cell, perc_var_nuc, perc_var_cyto], curr_compartment)
        self.max_val_compartment = get_val_compartment(
            [max_val_cell, max_val_nuc, max_val_cyto], curr_compartment)
                
x_fit = np.linspace(0,1,num=200)
ccd_coeff_list = []
not_ccd_coeff_list = []
model_free_list = []
model_free_list_all = []
xvals = np.linspace(0,1,num=21)
ccd_pvals = []
not_ccd_pvals = []

var_fred, var_fgreen = [],[] # variance of mean FUCCI intensities
var_mt_unorm = []
var_cell, var_nuc, var_cyto, var_mt = [],[],[],[] # mean intensity variances per antibody
# var_cell_int, var_nuc_int, var_cyto_int, var_mt_int = [],[],[],[] # integrated intensity variances per antibody

# These are tuples of (perc_var, y_val_avgs) where perc_var is the value described in the comments below, and the y_val_avgs are the moving averages
perc_var_fred, perc_var_fgreen = [],[] # percent variance attributed to cell cycle (FUCCI colors)
perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = [],[],[],[] # percent variance attributed to cell cycle (mean POI intensities)
# perc_var_cell_int, perc_var_nuc_int, perc_var_cyto_int, perc_var_mt_int = [],[],[],[] # percent variance attributed to cell cycle (integrated POI intensities)
mvavg_xvals = [] # this is a 2D array of xvals (rows) for all the antibodies (columns) for the moving average calculation

for well in u_well_plates:
    plt.close('all')
#    well = 'H05_55405991'#GMNN well, used for testing
    curr_well_inds = pol_sort_well_plate==well
    curr_pol = pol_sort_norm_rev[curr_well_inds]
    curr_fred = pol_sort_fred[curr_well_inds]
    curr_fgreen = pol_sort_fgreen[curr_well_inds]

    print(well)

    # check if the protein is expressed in the cyto,nuc,or both (cell)
    # curr_locs = plate_well_to_ab[plate][well][4]
    # curr_pval = plate_well_to_ab[plate][well][5]
    # curr_compartment = plate_well_to_ab[plate][well][6].lower()
    # curr_ensg = plate_well_to_ab[plate][well][7]

    # Normalize FUCCI colors
    curr_fred_norm = curr_fred/np.max(curr_fred)
    curr_fgreen_norm = curr_fgreen/np.max(curr_fgreen)
    # Normalize the mean intensities
    curr_ab_cell, curr_ab_nuc, curr_ab_cyto, curr_mt_cell = pol_sort_ab_cell[curr_well_inds], pol_sort_ab_nuc[curr_well_inds],pol_sort_ab_cyto[curr_well_inds], pol_sort_mt_cell[curr_well_inds]
    curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm, curr_mt_cell_norm = curr_ab_cell/np.max(curr_ab_cell), curr_ab_nuc/np.max(curr_ab_nuc), curr_ab_cyto/np.max(curr_ab_cyto), curr_mt_cell/np.max(curr_mt_cell)
    # Normalize the integrated intensities
    # curr_ab_cell_int, curr_ab_nuc_int, curr_ab_cyto_int, curr_mt_cell_int = pol_sort_ab_cell_int[curr_well_inds], pol_sort_ab_nuc_int[curr_well_inds],pol_sort_ab_cyto_int[curr_well_inds], pol_sort_mt_cell_int[curr_well_inds]
    # curr_ab_cell_norm_int, curr_ab_nuc_norm_int, curr_ab_cyto_norm_int, curr_mt_cell_norm_int = curr_ab_cell_int/np.max(curr_ab_cell_int), curr_ab_nuc_int/np.max(curr_ab_nuc_int), curr_ab_cyto_int/np.max(curr_ab_cyto_int), curr_mt_cell_int/np.max(curr_mt_cell_int)

    # Variance calculation FUCCI
    var_fred.append(np.var(curr_fred_norm))
    var_fgreen.append(np.var(curr_fgreen_norm))
    # Variance calculation -- mean intensities
    var_cell.append(np.var(curr_ab_cell_norm))
    var_nuc.append(np.var(curr_ab_nuc_norm))
    var_cyto.append(np.var(curr_ab_cyto_norm))
    var_mt.append(np.var(curr_mt_cell_norm))
    var_mt_unorm.append(np.var(curr_mt_cell))
    # Variance calculation -- integrated intensities
    # var_cell_int.append(np.var(curr_ab_cell_norm_int))
    # var_nuc_int.append(np.var(curr_ab_nuc_norm_int))
    # var_cyto_int.append(np.var(curr_ab_cyto_norm_int))
    # var_mt_int.append(np.var(curr_mt_cell_norm_int))

    # Compute Percent var, # if WINDOW>0: # always true for this use case, so I removed the other case
    perc_var_fred.append(mvavg_perc_var(curr_fred_norm, WINDOW))
    perc_var_fgreen.append(mvavg_perc_var(curr_fgreen_norm, WINDOW))
    # Compute percent variance due to the cell cycle (mean intensities)
    perc_var_cell.append(mvavg_perc_var(curr_ab_cell_norm,WINDOW))
    perc_var_nuc.append(mvavg_perc_var(curr_ab_nuc_norm, WINDOW))
    perc_var_cyto.append(mvavg_perc_var(curr_ab_cyto_norm, WINDOW))
    perc_var_mt.append(mvavg_perc_var(curr_mt_cell, WINDOW))
    # Compute percent variance due to the cell cycle (integrated intensities)
    # perc_var_cell_int.append(mvavg_perc_var(curr_ab_cell_norm_int,WINDOW))
    # perc_var_nuc_int.append(mvavg_perc_var(curr_ab_nuc_norm_int, WINDOW))
    # perc_var_cyto_int.append(mvavg_perc_var(curr_ab_cyto_norm_int, WINDOW))
    # perc_var_mt_int.append(mvavg_perc_var(curr_mt_cell_int, WINDOW))
    # Get x values for the moving average
    mvavg_xvals.append(mvavg_perc_var(curr_pol, WINDOW))

#%% Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
# Idea: create cutoffs for percent variance and 
# Execution: create cutoffs for perc_var and total variance per compartment and for integrated intensity and mean intensity
# Output: Graphs that illustrate the cutoffs (integrated, mean)
# Output: Overlap of the total variance cutoffs with the original filtering done manually

total_var_cutoff = np.mean(var_mt) + 1 * np.std(var_mt)
# total_var_cutoff_int = np.mean(var_mt_int) + 2 * np.std(var_mt_int)
percent_var_cutoff = np.mean([a[0] for a in perc_var_mt if not np.isinf(a[0])]) + 2 * np.std([a[0] for a in perc_var_mt if not np.isinf(a[0])])
# percent_var_cutoff_int = np.mean([a[0] for a in perc_var_mt_int if not np.isinf(a[0])]) + 2 * np.std([a[0] for a in perc_var_mt_int if not np.isinf(a[0])])
print(f"{total_var_cutoff}: cutoff for total variance")
print(f"{percent_var_cutoff}: cutoff for percent of total variance due to cell cycle")

plt.scatter(var_cell, [a[0] for a in perc_var_cell], c=wp_cell_kruskal_adj)
plt.vlines(total_var_cutoff, 0, 0.9)
plt.hlines(percent_var_cutoff, 0, 0.1)
plt.xlabel("Total Variance of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Cell - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CellProteinFractionVariance.png")
plt.show()
plt.close()
plt.scatter(var_cyto, [a[0] for a in perc_var_cyto], c=wp_cyto_kruskal_adj)
plt.vlines(total_var_cutoff, 0, 0.9)
plt.hlines(percent_var_cutoff, 0, 0.1)
plt.xlabel("Total Variance of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Cytoplasm - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CytoProteinFractionVariance.png")
plt.show()
plt.close()
plt.scatter(var_nuc, [a[0] for a in perc_var_nuc], c=wp_nuc_kruskal_adj)
plt.vlines(total_var_cutoff, 0, 0.9)
plt.hlines(percent_var_cutoff, 0, 0.1)
plt.xlabel("Total Variance of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Nucleoplasm - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/NucProteinFractionVariance.png")
plt.show()
plt.close()

# Output a list of the genes that show variation
name_df = pd.read_csv("input\\Fucci_staining_summary_first_plates.csv")
wppp, ensggg, abbb = list(name_df["well_plate"]), list(name_df["ENSG"]), list(name_df["Antibody"])
name_df2 = pd.read_csv("input\\FucciWellPlateGene.csv")
wppp.extend(name_df2["well_plate"])
ensggg.extend(name_df2["ENSG"])
abbb.extend(name_df2["Antibody"])
ensg_dict = dict([(wppp[i], ensggg[i]) for i in range(len(wppp))])
ab_dict = dict([(wppp[i], abbb[i]) for i in range(len(wppp))])
ENSG = np.asarray([ensg_dict[wp] if wp in ensg_dict else "" for wp in well_plate])
antibody = np.asarray([ab_dict[wp] if wp in ab_dict else "" for wp in well_plate])
wp_ensg = np.asarray([ensg_dict[wp] if wp in ensg_dict else "" for wp in u_well_plates])

alpha = 0.05
does_tot_vary_cell = (np.array(var_cell) >= total_var_cutoff) #| (np.array(var_cell_int) > total_var_cutoff_int)
does_tot_vary_nuc = (np.array(var_nuc) >= total_var_cutoff) #| (np.array(var_nuc_int) > total_var_cutoff_int)
does_tot_vary_cyto = (np.array(var_cyto) >= total_var_cutoff) #| (np.array(var_cyto_int) > total_var_cutoff_int)
does_perc_vary_cell = (np.array([a[0] for a in perc_var_cell]) >= percent_var_cutoff) #| (np.array([a[0] for a in perc_var_cell_int]) > percent_var_cutoff_int)
does_perc_vary_nuc = (np.array([a[0] for a in perc_var_nuc]) >= percent_var_cutoff) #| (np.array([a[0] for a in perc_var_nuc_int]) > percent_var_cutoff_int)
does_perc_vary_cyto = (np.array([a[0] for a in perc_var_cyto]) >= percent_var_cutoff) #| (np.array([a[0] for a in perc_var_cyto_int]) > percent_var_cutoff_int)
ccd_cell = does_tot_vary_cell & does_perc_vary_cell & (wp_cell_kruskal_adj < alpha)
ccd_cyto = does_tot_vary_cyto & does_perc_vary_cyto & (wp_cyto_kruskal_adj < alpha)
ccd_nuc = does_tot_vary_nuc & does_perc_vary_nuc & (wp_nuc_kruskal_adj < alpha)
ccd = ccd_cell | ccd_cyto | ccd_nuc
nonccd_cell = does_tot_vary_cell & ~ccd_cell
nonccd_cyto = does_tot_vary_cyto & ~ccd_cyto
nonccd_nuc = does_tot_vary_nuc & ~ccd_nuc
nonccd = nonccd_cell | nonccd_cyto | nonccd_nuc

df = pd.DataFrame({
    "well_plate" : u_well_plates, 
    "ENSG": wp_ensg,
    "var_cell":var_cell,
    "var_cyto":var_cyto,
    "var_nuc":var_nuc,
    "perc_var_cell":[a[0] for a in perc_var_cell],
    "perc_var_cyto":[a[0] for a in perc_var_cyto],
    "perc_var_nuc":[a[0] for a in perc_var_nuc],
    "wp_cell_kruskal_adj":wp_cell_kruskal_adj,
    "wp_cyto_kruskal_adj":wp_cyto_kruskal_adj,
    "wp_nuc_kruskal_adj":wp_nuc_kruskal_adj,
    "ccd_cell":ccd_cell,
    "ccd_cyto":ccd_cyto,
    "ccd_nuc":ccd_nuc,
    "ccd":ccd,
    "nonccd_cell":nonccd_cell,
    "nonccd_cyto":nonccd_cyto,
    "nonccd_nuc":nonccd_nuc,
    "nonccd":nonccd,
    })
df.to_csv("output/CellCycleVariationSummary.csv")
print(f"{sum(ccd)}: CCD variable proteins")
print(f"{sum(nonccd)}: non-CCD variable proteins")

#%%
