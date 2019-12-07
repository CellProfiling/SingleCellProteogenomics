#%% Imports
from imports import *
import numpy as np
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import iqr, variation

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
#plate = np.load("output/pickles/plate.npy", allow_pickle=True)
u_plate = np.load("output/pickles/u_plate.npy", allow_pickle=True)
well_plate=np.load("output/pickles/well_plate.npy", allow_pickle=True)
well_plate_imgnb=np.load("output/pickles/well_plate_imgnb.npy", allow_pickle=True)
u_well_plates=np.load("output/pickles/u_well_plates.npy", allow_pickle=True)
ab_nuc=np.load("output/pickles/ab_nuc.npy", allow_pickle=True)
ab_cyto=np.load("output/pickles/ab_cyto.npy", allow_pickle=True)
ab_cell=np.load("output/pickles/ab_cell.npy", allow_pickle=True)
mt_cell=np.load("output/pickles/mt_cell.npy", allow_pickle=True)
green_fucci=np.load("output/pickles/green_fucci.npy", allow_pickle=True)
red_fucci=np.load("output/pickles/red_fucci.npy", allow_pickle=True)
log_green_fucci_zeroc=np.load("output/pickles/log_green_fucci_zeroc.npy", allow_pickle=True)
log_red_fucci_zeroc=np.load("output/pickles/log_red_fucci_zeroc.npy", allow_pickle=True)
log_green_fucci_zeroc_rescale=np.load("output/pickles/log_green_fucci_zeroc_rescale.npy", allow_pickle=True)
log_red_fucci_zeroc_rescale=np.load("output/pickles/log_red_fucci_zeroc_rescale.npy", allow_pickle=True)
fucci_data = np.load("output/pickles/fucci_data.npy", allow_pickle=True)
wp_iscell = np.load("output/pickles/wp_iscell.npy", allow_pickle=True)
wp_isnuc = np.load("output/pickles/wp_isnuc.npy", allow_pickle=True)
wp_iscyto = np.load("output/pickles/wp_iscyto.npy", allow_pickle=True)
area_cell=np.load("output/pickles/area_cell.npy", allow_pickle=True)
area_nuc=np.load("output/pickles/area_nuc.npy", allow_pickle=True)
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
well_plate_imgnb_sort = well_plate_imgnb[pol_sort_inds]
ab_nuc_sort, ab_cyto_sort, ab_cell_sort, mt_cell_sort = pol_sort(pol_sort_inds,ab_nuc,ab_cyto,ab_cell,mt_cell)

# Rezero to minimum --reasoning, cells disappear during mitosis, so we should have the fewest detected cells there
bins = plt.hist(pol_sort_phi,NBINS)
start_phi = bins[1][np.argmin(bins[0])]

# Move those points to the other side
more_than_start = np.greater(pol_sort_phi,start_phi)
less_than_start = np.less_equal(pol_sort_phi,start_phi)
def pol_reord(arr):
    return np.concatenate((arr[more_than_start],arr[less_than_start]))

pol_sort_well_plate = pol_reord(well_plate_sort)
pol_sort_well_plate_imgnb = pol_reord(well_plate_imgnb_sort)
pol_sort_rho_reorder = pol_reord(pol_sort_rho)
pol_sort_inds_reorder = pol_reord(pol_sort_inds)
pol_sort_phi_reorder = np.concatenate((pol_sort_phi[more_than_start],pol_sort_phi[less_than_start]+np.pi*2))
pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell = pol_reord(ab_nuc_sort), pol_reord(ab_cyto_sort), pol_reord(ab_cell_sort), pol_reord(mt_cell_sort)
pol_sort_centered_data0, pol_sort_centered_data1 = pol_reord(centered_data_sort0), pol_reord(centered_data_sort1)
pol_sort_area_cell, pol_sort_area_nuc = pol_reord(area_cell), pol_reord(area_nuc)

#shift and re-scale "time"; reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
pol_sort_shift = pol_sort_phi_reorder+np.abs(np.min(pol_sort_phi_reorder))
pol_sort_norm = pol_sort_shift/np.max(pol_sort_shift)
pol_sort_norm_rev = 1-pol_sort_norm
pol_sort_norm_rev = stretch_time(pol_sort_norm_rev)

#apply uniform radius (rho) and convert back
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

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

#%% Measures of variance
# Idea: Calculate measures of variance, and show them in plots
# Execution: Now that we already have the data filtered for variability, this is just descriptive.
# Output: scatters of antibody vs microtubule variances by different measures of variaibility

# benjimini-hochberg multiple testing correction
# source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html

def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    # Written by: Olivia Guest github.com/oliviaguest/gini/blob/master/gini.py
    array = array.flatten()
    if np.amin(array) < 0: 
        array -= np.amin(array) # Values cannot be negative
    array = np.sort(array + 0.0000001) # Values must be sorted and nonzero
    index = np.arange(1, array.shape[0] + 1) # Index per array element
    n = array.shape[0] # Number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) # Gini coefficient

def values_comp(values_cell, values_nuc, values_cyto, wp_iscell, wp_isnuc, wp_iscyto):
    '''Get the values for the annotated compartment'''
    values_comp = np.empty_like(values_cell)
    values_comp[wp_iscell] = np.array(values_cell)[wp_iscell]
    values_comp[wp_isnuc] = np.array(values_nuc)[wp_isnuc]
    values_comp[wp_iscyto] = np.array(values_cyto)[wp_iscyto]
    return np.array(values_comp)

# Calculate variation
use_log = False
var_cell, var_nuc, var_cyto, var_mt = [],[],[],[] # mean intensity variances per antibody
cv_cell, cv_nuc, cv_cyto, cv_mt = [], [], [], []
gini_cell, gini_nuc, gini_cyto, gini_mt = [],[],[],[] # mean intensity ginis per antibody
var_cell_test_p, var_nuc_test_p, var_cyto_test_p = [],[],[]
mean_mean_cell, mean_mean_nuc, mean_mean_cyto, mean_mean_mt = [],[],[],[] # mean mean-intensity
cell_counts = []

wpi_img = []
var_cell_img, var_nuc_img, var_cyto_img, var_mt_img = [],[],[],[] # mean intensity variances per field of view
cv_cell_img, cv_nuc_img, cv_cyto_img, cv_mt_img = [], [], [], []

for well in u_well_plates:
    curr_well_inds = pol_sort_well_plate==well
    curr_ab_cell = pol_sort_ab_cell[curr_well_inds] if not use_log else np.log10(pol_sort_ab_cell[curr_well_inds])
    curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds] if not use_log else np.log10(pol_sort_ab_nuc[curr_well_inds])
    curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds] if not use_log else np.log10(pol_sort_ab_cyto[curr_well_inds])
    curr_mt_cell = pol_sort_mt_cell[curr_well_inds] if not use_log else np.log10(pol_sort_mt_cell[curr_well_inds])
    curr_fuccigreen=0
    curr_fuccired=0
    
    cell_counts.append(len(curr_ab_cell))
    
    var_cell.append(np.var(curr_ab_cell))
    var_nuc.append(np.var(curr_ab_nuc))
    var_cyto.append(np.var(curr_ab_cyto))
    var_mt.append(np.var(curr_mt_cell))
    
    cv_cell.append(variation(curr_ab_cell))
    cv_nuc.append(variation(curr_ab_nuc))
    cv_cyto.append(variation(curr_ab_cyto))
    cv_mt.append(variation(curr_mt_cell))
    
    gini_cell.append(gini(curr_ab_cell))
    gini_nuc.append(gini(curr_ab_nuc))
    gini_cyto.append(gini(curr_ab_cyto))
    gini_mt.append(gini(curr_mt_cell))
    
    # Save the mean mean intensities
    mean_mean_cell.append(np.mean(curr_ab_cell))
    mean_mean_nuc.append(np.mean(curr_ab_nuc)) 
    mean_mean_cyto.append(np.mean(curr_ab_cyto))
    mean_mean_mt.append(np.mean(curr_mt_cell))
    
    curr_well_plate_imgnbs = pol_sort_well_plate_imgnb[curr_well_inds]
    curr_wpi_img = []
    curr_var_cell_img, curr_var_nuc_img, curr_var_cyto_img, curr_var_mt_img = [],[],[],[] # mean intensity variances per field of view
    curr_cv_cell_img, curr_cv_nuc_img, curr_cv_cyto_img, curr_cv_mt_img = [], [], [], []
    for wpi in np.unique(curr_well_plate_imgnbs):
        curr_wpis = pol_sort_well_plate_imgnb == wpi
        curr_ab_cell = pol_sort_ab_cell[curr_wpis] if not use_log else np.log10(pol_sort_ab_cell[curr_wpis])
        curr_ab_nuc = pol_sort_ab_nuc[curr_wpis] if not use_log else np.log10(pol_sort_ab_nuc[curr_wpis])
        curr_ab_cyto = pol_sort_ab_cyto[curr_wpis] if not use_log else np.log10(pol_sort_ab_cyto[curr_wpis])
        curr_mt_cell = pol_sort_mt_cell[curr_wpis] if not use_log else np.log10(pol_sort_mt_cell[curr_wpis])
        
        curr_wpi_img.append(wpi)
        
        curr_var_cell_img.append(np.var(curr_ab_cell))
        curr_var_nuc_img.append(np.var(curr_ab_nuc))
        curr_var_cyto_img.append(np.var(curr_ab_cyto))
        curr_var_mt_img.append(np.var(curr_mt_cell))
        
        curr_cv_cell_img.append(variation(curr_ab_cell))
        curr_cv_nuc_img.append(variation(curr_ab_nuc))
        curr_cv_cyto_img.append(variation(curr_ab_cyto))
        curr_cv_mt_img.append(variation(curr_mt_cell))
    
    wpi_img.append(curr_wpi_img)
    var_cell_img.append(curr_var_cell_img)
    var_nuc_img.append(curr_var_nuc_img)
    var_cyto_img.append(curr_var_cyto_img)
    var_mt_img.append(curr_var_mt_img)
    
    cv_cell_img.append(curr_cv_cell_img)
    cv_nuc_img.append(curr_cv_nuc_img)
    cv_cyto_img.append(curr_cv_cyto_img)
    cv_mt_img.append(curr_cv_mt_img)

# looking at the mean mean intensities
mean_mean_comp = values_comp(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, wp_iscell, wp_isnuc, wp_iscyto)
firstbatch = np.asarray([not str(p).split("_")[1].startswith("67") for p in u_well_plates])
plt.hist(np.array(mean_mean_comp)[firstbatch], bins=200, label="firstbatch")
plt.hist(np.array(mean_mean_comp)[~firstbatch], bins=200, label="secondbatch")
plt.legend()
plt.savefig("figures/MeanMeancComp.png")
plt.show()
plt.close()

firstbatch = np.asarray([not str(p).split("_")[1].startswith("67") for p in u_well_plates])
plt.hist(np.array(mean_mean_mt)[firstbatch], bins=200, label="firstbatch")
plt.hist(np.array(mean_mean_mt)[~firstbatch], bins=200, label="secondbatch")
plt.legend()
plt.savefig("figures/MeanMeanMt.png")
plt.show()
plt.close()

var_cell, var_nuc, var_cyto, var_mt = np.array(var_cell), np.array(var_nuc), np.array(var_cyto), np.array(var_mt)
gini_cell, gini_nuc, gini_cyto, gini_mt = np.array(gini_cell), np.array(gini_nuc), np.array(gini_cyto), np.array(gini_mt)
cv_cell, cv_nuc, cv_cyto, cv_mt = np.array(cv_cell), np.array(cv_nuc), np.array(cv_cyto), np.array(cv_mt)

mmmm = np.concatenate((var_cell, var_cyto, var_nuc, var_mt))
cccc = (["var_cell"] * len(var_cell))
cccc.extend(["var_cyto"] * len(var_cyto))
cccc.extend(["var_nuc"] * len(var_nuc))
cccc.extend(["var_mt"] * len(var_mt))
moddf = pd.DataFrame({"variance": mmmm, "category" : cccc})
boxplot = moddf.boxplot("variance", by="category", figsize=(12, 8), showfliers=True)
boxplot.set_xlabel("Metacompartment", size=36,fontname='Arial')
boxplot.set_ylabel(f"Variance using {'log' if use_log else 'natural'} intensity values", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/VarianceBoxplot.png")
plt.show()
plt.close()

mmmm = np.concatenate((gini_cell, gini_cyto, gini_nuc, gini_mt))
cccc = (["gini_cell"] * len(gini_cell))
cccc.extend(["gini_cyto"] * len(gini_cyto))
cccc.extend(["gini_nuc"] * len(gini_nuc))
cccc.extend(["gini_mt"] * len(gini_mt))
moddf = pd.DataFrame({"variance": mmmm, "category" : cccc})
boxplot = moddf.boxplot("variance", by="category", figsize=(12, 8), showfliers=True)
boxplot.set_xlabel("Metacompartment", size=36,fontname='Arial')
boxplot.set_ylabel(f"Gini using {'log' if use_log else 'natural'} intensity values", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/GiniBoxplot.png")
plt.show()
plt.close()

mean_mean_comp = values_comp(mean_mean_cell, mean_mean_nuc, mean_mean_cyto, wp_iscell, wp_isnuc, wp_iscyto)
cv_comp = values_comp(cv_cell, cv_nuc, cv_cyto, wp_iscell, wp_isnuc, wp_iscyto)
gini_comp = values_comp(gini_cell, gini_nuc, gini_cyto, wp_iscell, wp_isnuc, wp_iscyto)
var_comp = values_comp(var_cell, var_nuc, var_cyto, wp_iscell, wp_isnuc, wp_iscyto)

# var vs mt var
plt.figure(figsize=(10,10))
plt.scatter(var_mt, var_comp, label="all")
plt.ylabel(f"var Comp")
plt.xlabel(f"var microtubules")
plt.legend()
plt.savefig(f"figures/varCompVsvarMt_knowns.png")
plt.show()
plt.close()

# CV vs mt cv
plt.figure(figsize=(10,10))
plt.scatter(cv_mt, cv_comp,  label="all")
plt.ylabel(f"cv Comp")
plt.xlabel(f"cv microtubules")
plt.legend()
plt.savefig(f"figures/cvCompVsCvMt.png")
plt.show()
plt.close()

# Gini vs mt gini
plt.figure(figsize=(10,10))
plt.scatter(gini_comp, gini_mt, label="all")
plt.xlabel(f"gini Comp")
plt.ylabel(f"gini microtubules")
plt.legend()
plt.savefig(f"figures/GiniCompVsGiniMt.png")
plt.show()
plt.close()

# variance vs intensity
plt.figure(figsize=(10,10))
plt.scatter(var_comp, mean_mean_comp, label="all")
plt.ylabel(f"{'log10' if use_log else 'natural'} intensity")
plt.xlabel(f"variance Comp")
plt.legend()
plt.savefig(f"figures/VarianceVsIntensityComp.png")
plt.show()
plt.close()

#%% How do the variances of antibody intensities in images compare to the variance of the samples?
var_comp_img = values_comp(var_cell_img, var_nuc_img, var_cyto_img, wp_iscell, wp_isnuc, wp_iscyto)
cv_comp_img = values_comp(cv_cell_img, cv_nuc_img, cv_cyto_img, wp_iscell, wp_isnuc, wp_iscyto)

plt.scatter(np.concatenate([[var_comp[i]] * len(vvv) for i, vvv in enumerate(var_comp_img)]), np.concatenate(var_comp_img))
plt.xlabel("variance within compartment")
plt.ylabel("variance for each image")
plt.savefig("figures/VarianceByImage.png")
plt.show()
plt.close()

plt.scatter(np.concatenate([[cv_comp[i]] * len(vvv) for i, vvv in enumerate(cv_comp_img)]), np.concatenate(cv_comp_img))
plt.xlabel("CV within compartment")
plt.ylabel("CV for each image")
plt.savefig("figures/CVByImage.png")
plt.show()
plt.close()

print(np.concatenate(wpi_img)[np.argmax(np.concatenate(var_comp_img))] + ": the image with the max variance")

plt.hist(np.concatenate([vvv / var_comp[i] for i, vvv in enumerate(var_comp_img)]))
plt.show();plt.close();
high_var_img = np.concatenate(wpi_img)[np.concatenate([vvv > 4 * var_comp[i] for i, vvv in enumerate(var_comp_img)])]
print(f"{high_var_img}: the images with greater than 4x the variance of the whole sample")

norm_cv_img = np.concatenate([vvv / cv_comp[i] for i, vvv in enumerate(cv_comp_img)])
plt.hist(norm_cv_img)
plt.show();plt.close();
cutoff = np.mean(norm_cv_img) + 3 * np.std(norm_cv_img)
high_cv_img = np.concatenate(wpi_img)[norm_cv_img > cutoff]
print(f"{high_cv_img}: the images with greater than 4x the variance of the whole sample")

np.intersect1d(high_var_img, high_cv_img)

#%%
def np_save_overwriting(fn, arr):
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)
        
np_save_overwriting("output/pickles/pol_sort_well_plate.npy", pol_sort_well_plate)
np_save_overwriting("output/pickles/pol_sort_norm_rev.npy", pol_sort_norm_rev)
np_save_overwriting("output/pickles/pol_sort_ab_nuc.npy", pol_sort_ab_nuc)
np_save_overwriting("output/pickles/pol_sort_ab_cyto.npy", pol_sort_ab_cyto)
np_save_overwriting("output/pickles/pol_sort_ab_cell.npy", pol_sort_ab_cell)
np_save_overwriting("output/pickles/pol_sort_mt_cell.npy", pol_sort_mt_cell)
np_save_overwriting("output/pickles/pol_sort_area_cell.npy", pol_sort_area_cell)
np_save_overwriting("output/pickles/pol_sort_area_nuc.npy", pol_sort_area_nuc)

np_save_overwriting("output/pickles/mean_mean_comp.npy", mean_mean_comp)
np_save_overwriting("output/pickles/cv_comp.npy", cv_comp)
np_save_overwriting("output/pickles/gini_comp.npy", gini_comp)
np_save_overwriting("output/pickles/var_comp.npy", var_comp)