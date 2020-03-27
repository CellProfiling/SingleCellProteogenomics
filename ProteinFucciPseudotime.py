# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:20:24 2020

@author: antho
"""

from utils import *
import utils
from stretch_time import stretch_time
import numpy as np
import scipy.optimize
import scipy.stats
import sklearn.mixture
import seaborn as sbn
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

NBINS = 150 #number of bins, arbitrary choice for now
WINDOW_FUCCI_PSEUDOTIME = 100

G1_LEN = 10.833 #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
G1_S_TRANS = 2.625 #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
S_G2_LEN = 11.975 #hours (this should be from the G2/M cutoff above to the end)
#M_LEN = 0.5
#We are excluding Mphase from this analysis
TOT_LEN = G1_LEN+G1_S_TRANS+S_G2_LEN

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

def pol_sort(inds, nuc, cyto, cell, mt):
    '''Sort data by polar coordinates'''
    return nuc[inds], cyto[inds], cell[inds], mt[inds]

def pol_reord(arr, more_than_start, less_than_start):
    '''Reorder an array based on the start position of the polar coordinate model'''
    return np.concatenate((arr[more_than_start],arr[less_than_start]))

def pol2cart(rho, phi):
    '''Apply uniform radius (rho) and convert back'''
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

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
    plt.show()
    plt.close()
    
def visualize_gmnn_agreement(centered_data, pol_sort_well_plate, pol_sort_ab_nuc, pol_sort_centered_data0, pol_sort_centered_data1, cart_data_ur):
    '''GMNN was tagged in the FUCCI cells and stained with antibodies; visualize the agreement'''
    fig, ax1 = plt.subplots(figsize=(10,10))
    mycmap = plt.cm.gray_r
    mycmap.set_under(color='w',alpha=None)
    ax1.hist2d(centered_data[:,0],centered_data[:,1],bins=200,alpha=1,cmap=mycmap)
    hist, xbins, ybins = np.histogram2d(cart_data_ur[0],cart_data_ur[1], bins=200, normed=True)
    extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
    im = ax1.imshow(
            np.ma.masked_where(hist == 0, hist).T,
            interpolation='nearest',
            origin='lower',
            extent=extent,
            cmap='plasma')
    gmnn = "H05_55405991"
    gmnn_well_inds = pol_sort_well_plate==gmnn
    gmnn_ab_nuc = pol_sort_ab_nuc[gmnn_well_inds]
    im = ax1.scatter(pol_sort_centered_data0[gmnn_well_inds],pol_sort_centered_data1[gmnn_well_inds], c=gmnn_ab_nuc)
    fig.colorbar(im, ax=ax1)
    plt.xlabel(r'$log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(r'$log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.savefig("figures/GMNN_FUCCI_plot.pdf")
    plt.show()
    plt.close()

def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')
def mvpercentiles(yvals_binned):
    return np.percentile(yvals_binned, [10, 25, 50, 75, 90], axis=1)


def plot_fucci_intensities_on_pseudotime(pol_sort_norm_rev, pol_sort_centered_data1, pol_sort_centered_data0):
    '''visualize FUCCI intensities over pseudotime'''
    plt.figure(figsize=(5,5))
    WINDOW_FUCCI_PSEUDOTIMEs = np.asarray([np.arange(start, start + WINDOW_FUCCI_PSEUDOTIME) for start in np.arange(len(pol_sort_norm_rev) - WINDOW_FUCCI_PSEUDOTIME + 1)])
    mvperc_red = mvpercentiles(pol_sort_centered_data1[WINDOW_FUCCI_PSEUDOTIMEs])
    mvperc_green = mvpercentiles(pol_sort_centered_data0[WINDOW_FUCCI_PSEUDOTIMEs])
    mvavg_xvals = mvavg(pol_sort_norm_rev, WINDOW_FUCCI_PSEUDOTIME)
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc_green[1], mvperc_green[-2], color="lightgreen", label="25th & 75th Percentiles")
    plt.fill_between(mvavg_xvals * TOT_LEN, mvperc_red[1], mvperc_red[-2], color="lightcoral", label="25th & 75th Percentiles")
    
    mvavg_red = mvavg(pol_sort_centered_data1, WINDOW_FUCCI_PSEUDOTIME)
    mvavg_green = mvavg(pol_sort_centered_data0, WINDOW_FUCCI_PSEUDOTIME)
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_red, color="r", label="Mean Intensity")
    plt.plot(mvavg_xvals * TOT_LEN, mvavg_green, color="g", label="Mean Intensity")
    plt.xlabel('Cell Cycle Time, hrs')
    plt.ylabel('Log10 Tagged CDT1 & GMNN Intensity')
    plt.xticks(size=14)
    plt.yticks(size=14)
    # plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig("figures/FUCCIOverPseudotime.pdf")
    plt.savefig("figures/FUCCIOverPseudotime.png")
    plt.show()
    plt.close()
    
def fucci_polar_coordinate_calculations(fucci_data, ab_nuc,ab_cyto,ab_cell,mt_cell,area_cell, area_nuc,well_plate,well_plate_imgnb,
                        log_red_fucci_zeroc_rescale,log_green_fucci_zeroc_rescale):
    '''Generate a polar coordinate model of cell cycle progression based on the FUCCI intensities'''
    # Center data
    x0 = np.ones(5)
    x = fucci_data[:,0]
    y = fucci_data[:,1]
    center_estimate = np.mean(fucci_data[:,0]), np.mean(fucci_data[:,1])
    center_2 = scipy.optimize.least_squares(f_2, center_estimate, args=(x, y))
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
    
    well_plate_sort = well_plate[pol_sort_inds]
    well_plate_imgnb_sort = well_plate_imgnb[pol_sort_inds]
    ab_nuc_sort, ab_cyto_sort, ab_cell_sort, mt_cell_sort = pol_sort(pol_sort_inds,ab_nuc,ab_cyto,ab_cell,mt_cell)
    
    # Rezero to minimum --reasoning, cells disappear during mitosis, so we should have the fewest detected cells there
    bins = plt.hist(pol_sort_phi,NBINS)
    start_phi = bins[1][np.argmin(bins[0])]
    
    # Move those points to the other side
    more_than_start = np.greater(pol_sort_phi,start_phi)
    less_than_start = np.less_equal(pol_sort_phi,start_phi)
    
    pol_sort_well_plate = pol_reord(well_plate_sort, more_than_start, less_than_start)
    pol_sort_well_plate_imgnb = pol_reord(well_plate_imgnb_sort, more_than_start, less_than_start)
    pol_sort_rho_reorder = pol_reord(pol_sort_rho, more_than_start, less_than_start)
    pol_sort_inds_reorder = pol_reord(pol_sort_inds, more_than_start, less_than_start)
    pol_sort_phi_reorder = np.concatenate((pol_sort_phi[more_than_start],pol_sort_phi[less_than_start]+np.pi*2))
    pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_ab_cell, pol_sort_mt_cell = pol_reord(ab_nuc_sort, more_than_start, less_than_start), pol_reord(ab_cyto_sort, more_than_start, less_than_start), pol_reord(ab_cell_sort, more_than_start, less_than_start), pol_reord(mt_cell_sort, more_than_start, less_than_start)
    pol_sort_centered_data0, pol_sort_centered_data1 = pol_reord(centered_data_sort0, more_than_start, less_than_start), pol_reord(centered_data_sort1, more_than_start, less_than_start)
    pol_sort_area_cell, pol_sort_area_nuc = pol_reord(area_cell, more_than_start, less_than_start), pol_reord(area_nuc, more_than_start, less_than_start)
    pol_sort_fred, pol_sort_fgreen = pol_reord(log_red_fucci_zeroc_rescale, more_than_start, less_than_start), pol_reord(log_green_fucci_zeroc_rescale, more_than_start, less_than_start)
    
    # shift and re-scale "time"; reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
    pol_sort_shift = pol_sort_phi_reorder+np.abs(np.min(pol_sort_phi_reorder))
    pol_sort_norm = pol_sort_shift/np.max(pol_sort_shift)
    pol_sort_norm_rev = 1-pol_sort_norm
    pol_sort_norm_rev = stretch_time(pol_sort_norm_rev)
    
    cart_data_ur = pol2cart(np.repeat(R_2,len(centered_data)), pol_data[1])
    
    # visualize that result
    start_pt = pol2cart(R_2,start_phi)
    fucci_hist2d(centered_data,cart_data_ur,start_pt)
    visualize_gmnn_agreement(centered_data, pol_sort_well_plate, pol_sort_ab_nuc, pol_sort_centered_data0, pol_sort_centered_data1,cart_data_ur)
    plot_fucci_intensities_on_pseudotime(pol_sort_norm_rev, pol_sort_centered_data1, pol_sort_centered_data0)
    
    # pickle the results
    utils.np_save_overwriting("output/pickles/pol_sort_well_plate.npy", pol_sort_well_plate)
    utils.np_save_overwriting("output/pickles/pol_sort_norm_rev.npy", pol_sort_norm_rev)
    utils.np_save_overwriting("output/pickles/pol_sort_ab_nuc.npy", pol_sort_ab_nuc)
    utils.np_save_overwriting("output/pickles/pol_sort_ab_cyto.npy", pol_sort_ab_cyto)
    utils.np_save_overwriting("output/pickles/pol_sort_ab_cell.npy", pol_sort_ab_cell)
    utils.np_save_overwriting("output/pickles/pol_sort_mt_cell.npy", pol_sort_mt_cell)
    utils.np_save_overwriting("output/pickles/pol_sort_area_cell.npy", pol_sort_area_cell)
    utils.np_save_overwriting("output/pickles/pol_sort_area_nuc.npy", pol_sort_area_nuc)
    utils.np_save_overwriting("output/pickles/pol_sort_fred.npy", pol_sort_fred)
    utils.np_save_overwriting("output/pickles/pol_sort_fgreen.npy", pol_sort_fgreen)
    
