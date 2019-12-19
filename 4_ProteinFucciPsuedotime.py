#%% Imports
from imports import *
import numpy as np
import os
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import iqr, variation
from scipy.stats import linregress
from sklearn.mixture import GaussianMixture

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
wp_comp_kruskal_gaussccd_adj = np.load("output/pickles/wp_comp_kruskal_gaussccd_adj.npy", allow_pickle=True)
wp_pass_kruskal_gaussccd_bh_comp = np.load("output/pickles/wp_pass_kruskal_gaussccd_bh_comp.npy", allow_pickle=True)
fucci_data = np.load("output/pickles/fucci_data.npy", allow_pickle=True)

wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True) 
wp_ab = np.load("output/pickles/wp_ab.npy", allow_pickle=True) 
wp_prev_ccd = np.load("output/pickles/wp_prev_ccd.npy", allow_pickle=True) 
wp_prev_notccd = np.load("output/pickles/wp_prev_notccd.npy", allow_pickle=True) 
wp_prev_negative = np.load("output/pickles/wp_prev_negative.npy", allow_pickle=True) 
prev_ccd_ensg = np.load("output/pickles/prev_ccd_ensg.npy", allow_pickle=True) 
prev_notccd_ensg = np.load("output/pickles/prev_notccd_ensg.npy", allow_pickle=True) 
prev_negative_ensg = np.load("output/pickles/prev_negative_ensg.npy", allow_pickle=True)

u_well_plates_old = np.load("output/pickles/u_well_plates.devin.npy", allow_pickle=True)
perc_var_compartment_old = np.load("output/pickles/perc_var_compartment.devin.npy", allow_pickle=True)
perc_var_cell_old = np.load("output/pickles/perc_var_cell.devin.npy", allow_pickle=True)
perc_var_nuc_old = np.load("output/pickles/perc_var_nuc.devin.npy", allow_pickle=True)
perc_var_cyto_old = np.load("output/pickles/perc_var_cyto.devin.npy", allow_pickle=True)

mean_mean_comp = np.load("output/pickles/mean_mean_comp.npy", allow_pickle=True)
cv_comp = np.load("output/pickles/cv_comp.npy", allow_pickle=True)
gini_comp = np.load("output/pickles/gini_comp.npy", allow_pickle=True)
var_comp = np.load("output/pickles/var_comp.npy", allow_pickle=True)

pol_sort_well_plate = np.load("output/pickles/pol_sort_well_plate.npy", allow_pickle=True)
pol_sort_norm_rev = np.load("output/pickles/pol_sort_norm_rev.npy", allow_pickle=True)
pol_sort_ab_nuc = np.load("output/pickles/pol_sort_ab_nuc.npy", allow_pickle=True)
pol_sort_ab_cyto = np.load("output/pickles/pol_sort_ab_cyto.npy", allow_pickle=True)
pol_sort_ab_cell = np.load("output/pickles/pol_sort_ab_cell.npy", allow_pickle=True)
pol_sort_mt_cell = np.load("output/pickles/pol_sort_mt_cell.npy", allow_pickle=True)
pol_sort_area_cell=np.load("output/pickles/pol_sort_area_cell.npy", allow_pickle=True)
pol_sort_area_nuc=np.load("output/pickles/pol_sort_area_nuc.npy", allow_pickle=True)

wp_iscell = np.load("output/pickles/wp_iscell.npy", allow_pickle=True)
wp_isnuc = np.load("output/pickles/wp_isnuc.npy", allow_pickle=True)
wp_iscyto = np.load("output/pickles/wp_iscyto.npy", allow_pickle=True)
print("loaded")

#%% Gaussian clustering to identify biomodal intensity distributions
def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection'''
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)

def benji_hoch(alpha, pvals):
    pvals = np.nan_to_num(pvals, nan=1) # fail the ones with not enough data
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
    pvals_corrected_BH = np.empty_like(pvals_corrected)

    # deal with sorting
    pvals_corrected_BH[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_BH = np.empty_like(reject)
    reject_BH[pvals_sortind] = reject
    return pvals_corrected_BH, reject_BH

wp_bimodal_cluster_idxs = []
wp_bimodal_diffmeans = []
wp_bimodal_fcmeans = []
wp_bimodal_clusterlabels = []
wp_isbimodal_p = []
wp_timebimodal_p = []

analysis = "MeanRng"
folder = f"figures/TemporalMovingAverages{analysis}191205"
if not os.path.exists(folder): os.mkdir(folder)
fileprefixes = np.array([f"{ensg}_{sum(wp_ensg[:ei] == ensg)}" for ei, ensg in enumerate(wp_ensg)])

gaussian = GaussianMixture(n_components=2, random_state=1, max_iter=500)
for i, well in enumerate(u_well_plates):
    curr_well_inds = pol_sort_well_plate==well # the reversal isn't really helpful here
    curr_pol = pol_sort_norm_rev[curr_well_inds]
    curr_ab_cell = pol_sort_ab_cell[curr_well_inds]
    curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds]
    curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds]
    curr_mt_cell = pol_sort_mt_cell[curr_well_inds]

    # Normalize mean intensities, normalized for display
    curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
    curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
    curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
    curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
    curr_comp_norm = np.asarray(curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm)

    cluster_labels = gaussian.fit_predict(curr_comp_norm.reshape(1, -1).T)
#    cluster_labels = gaussian.fit_predict(np.array([curr_pol, curr_comp_norm]).T)
    wp_bimodal_clusterlabels.append(cluster_labels)
    c1 = cluster_labels == 0
    c2 = cluster_labels == 1
    wp_bimodal_cluster_idxs.append([c1, c2])
    wp_bimodal_diffmeans.append(np.mean(curr_comp_norm[c2]) - np.mean(curr_comp_norm[c1]))
    wp_bimodal_fcmeans.append(np.mean(curr_comp_norm[c2]) / np.mean(curr_comp_norm[c1]))
    
    k, p = scipy.stats.kruskal(curr_comp_norm[c1], curr_comp_norm[c2])
    wp_isbimodal_p.append(p)
    k, p = scipy.stats.kruskal(curr_pol[c1], curr_pol[c2])
    wp_timebimodal_p.append(p)
    
wp_isbimodal_padj, wp_isbimodal_pass = benji_hoch(0.01, wp_isbimodal_p)
wp_timebimodal_padj, wp_timebimodal_pass = benji_hoch(0.01, wp_timebimodal_p)

wp_enoughcellsinbothclusters = np.array([sum(c1[0]) > 50 and sum(c1[1]) > 50 for c1 in wp_bimodal_cluster_idxs])
wp_isbimodal_fcpadj_pass = (np.abs(np.log(wp_bimodal_fcmeans) / np.log(2)) > 1) & wp_isbimodal_pass & ~wp_timebimodal_pass & wp_enoughcellsinbothclusters

plt.scatter(np.log(wp_bimodal_fcmeans) / np.log(2), -np.log10(wp_isbimodal_padj), c=wp_isbimodal_fcpadj_pass)
plt.show();plt.close()
plt.scatter([sum(c1[0]) for c1 in wp_bimodal_cluster_idxs], [sum(c1[1]) for c1 in wp_bimodal_cluster_idxs], c=wp_enoughcellsinbothclusters)
plt.show();plt.close()

bimodal = "figures/bimodal"
if not os.path.exists(bimodal): os.mkdir(bimodal)
for ensg in wp_ensg[wp_isbimodal_fcpadj_pass]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(bimodal, ensg+'_mvavg.png'))
        
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
WINDOW = 10 #Number of points for moving average window, arbitrary choice
DO_PLOTS = True #flag of whether to plot each well and save the plot
TPLOT_MODE = 'avg' #can have values of: 'avg', 'psin'
HIGHLIGHTS = ['ORC6','DUSP19','BUB1B','DPH2', 'FLI1']

def mvavg_perc_var(yvals,mv_window):
    yval_avg = np.convolve(yvals,np.ones((mv_window,))/mv_window, mode='valid')
    return np.var(yval_avg)/np.var(yvals), yval_avg

def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')

def mvmed(yvals_binned):
    return np.median(yvals_binned, axis=1)

def mvpercentiles(yvals_binned):
    return np.percentile(yvals_binned, [10, 25, 50, 75, 90], axis=1)

def mvmed_perc_var(yvals, windows):
    yval_avg = mvmed(yvals[windows])
    return np.var(yval_avg) / np.var(yvals), yval_avg

def temporal_mov_avg(curr_pol, curr_ab_norm, mvavg_xvals, mvavg_yvals, windows, folder, fileprefix):
    plt.close()
    outfile = os.path.join(folder,fileprefix+'_mvavg.png')
    if os.path.exists(outfile): return
    plt.figure(figsize=(5,5))
    mvperc = mvpercentiles(curr_ab_norm[windows])
    plt.fill_between(mvavg_xvals, mvperc[0], mvperc[-1], color="lightsteelblue", label="10th & 90th Percentiles")
    plt.fill_between(mvavg_xvals, mvperc[1], mvperc[-2], color="steelblue", label="25th & 75th Percentiles")
    plt.plot(mvavg_xvals, mvavg_yvals, color="blue", label="Mean Intensity")
    plt.scatter(curr_pol, curr_ab_norm, c='c')
    plt.xlabel('Pseudotime')
    plt.ylabel(fileprefix + ' Protein Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
#    plt.legend(fontsize=14)
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)
    plt.close()

# bonferroni MTC
def bonf(alpha, pvals):
    pvals = np.nan_to_num(pvals, nan=1) # fail the ones with not enough data
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)
    alphaBonf = alpha / float(len(pvals))
    rejectBonf = pvals_sorted <= alphaBonf
    pvals_correctedBonf = pvals_sorted * float(len(pvals))
    pvals_correctedBonf_unsorted = np.empty_like(pvals_correctedBonf) 
    pvals_correctedBonf_unsorted[pvals_sortind] = pvals_correctedBonf
    rejectBonf_unsorted = np.empty_like(rejectBonf)
    rejectBonf_unsorted[pvals_sortind] = rejectBonf
    return pvals_correctedBonf_unsorted, rejectBonf_unsorted

def values_comp(values_cell, values_nuc, values_cyto, wp_iscell, wp_isnuc, wp_iscyto):
    '''Get the values for the annotated compartment'''
    values_comp = np.empty_like(values_cell)
    values_comp[wp_iscell] = np.array(values_cell)[wp_iscell]
    values_comp[wp_isnuc] = np.array(values_nuc)[wp_isnuc]
    values_comp[wp_iscyto] = np.array(values_cyto)[wp_iscyto]
    return np.array(values_comp)

x_fit = np.linspace(0,1,num=200)
ccd_coeff_list = []
not_ccd_coeff_list = []
xvals = np.linspace(0,1,num=21)
ccd_pvals = []
not_ccd_pvals = []

use_log_ccd = False
perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = [],[],[],[] # percent variance attributed to cell cycle (mean POI intensities)
perc_var_comp_rng, perc_var_mt_rng = [],[] # randomized in pseudotime; percent variances
perc_var_comp_clust1, perc_var_comp_clust2, mvavgs_comp_clust1, mvavgs_comp_clust2, perc_var_comp_clust1_rng, perc_var_comp_clust2_rng, mvavgs_x_clust1, mvavgs_x_clust2 = [],[],[],[],[],[],[],[] # percent variances for bimodal
mvavgs_cell, mvavgs_nuc, mvavgs_cyto, mvavgs_mt, mvavgs_x = [],[],[],[],[] # moving average y values & x value
mvavgs_cell_rng, mvavgs_nuc_rng, mvavgs_cyto_rng, mvavgs_mt_rng = [],[],[],[] # moving average y values from randomization
ccd_var_cell_levenep, ccd_var_nuc_levenep, ccd_var_cyto_levenep = [],[],[] # two-tailed p values for equal variance of mvavg and raw values
cell_counts = []

for i, well in enumerate(u_well_plates):
#    print(well)
    plt.close('all')
    if i % 100 == 0: print(f"well {i} of {len(u_well_plates)}")
#    well = 'H05_55405991'#GMNN well, used for testing
    curr_well_inds = pol_sort_well_plate==well # the reversal isn't really helpful here
    curr_pol = pol_sort_norm_rev[curr_well_inds]
    curr_ab_cell = pol_sort_ab_cell[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_cell[curr_well_inds])
    curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_nuc[curr_well_inds])
    curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_ab_cyto[curr_well_inds])
    curr_mt_cell = pol_sort_mt_cell[curr_well_inds] if not use_log_ccd else np.log10(pol_sort_mt_cell[curr_well_inds])

    # Normalize mean intensities, normalized for display
    curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
    curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
    curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
    curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
    curr_area_cell_norm = pol_sort_area_cell[curr_well_inds] / np.max(pol_sort_area_cell[curr_well_inds])
    curr_area_nuc_norm = pol_sort_area_nuc[curr_well_inds] / np.max(pol_sort_area_nuc[curr_well_inds])
        
    # Original method from Devin's work
    perc_var_cell_val, mvavg_cell = mvavg_perc_var(curr_ab_cell_norm, WINDOW)
    perc_var_nuc_val, mvavg_nuc = mvavg_perc_var(curr_ab_nuc_norm, WINDOW)
    perc_var_cyto_val, mvavg_cyto = mvavg_perc_var(curr_ab_cyto_norm, WINDOW)
    perc_var_mt_val, mvavg_mt = mvavg_perc_var(curr_mt_cell_norm, WINDOW)
    mvavg_xvals = mvavg(curr_pol, WINDOW)
    
    perms = np.asarray([np.random.permutation(len(curr_pol)) for nnn in np.arange(1000)])
    curr_comp_norm = np.asarray(curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm)
    curr_comp_perm = np.asarray([curr_comp_norm[perm] for perm in perms])
    curr_mt_perm = np.asarray([curr_mt_cell_norm[perm] for perm in perms])
    curr_mvavg_rng_comp = np.apply_along_axis(mvavg, 1, curr_comp_perm, WINDOW)
    curr_mvavg_rng_mt = np.apply_along_axis(mvavg, 1, curr_mt_perm, WINDOW)
    curr_percvar_rng_comp = np.var(curr_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm, axis=1)
    curr_percvar_rng_mt = np.var(curr_mvavg_rng_mt, axis=1) / np.var(curr_mt_perm, axis=1)
    perc_var_comp_rng.append(curr_percvar_rng_comp)
    perc_var_mt_rng.append(curr_percvar_rng_mt)
    
    if wp_isbimodal_fcpadj_pass[i]:
        clust1_idx, clust2_idx = wp_bimodal_cluster_idxs[i]
        perc_var_comp_clust1_val, mvavg_clust1 = mvavg_perc_var(curr_comp_norm[clust1_idx], WINDOW)
        perc_var_comp_clust2_val, mvavg_clust2 = mvavg_perc_var(curr_comp_norm[clust2_idx], WINDOW)
        mvavgs_x_clust1.append(mvavg(curr_pol[clust1_idx], WINDOW))
        mvavgs_x_clust2.append(mvavg(curr_pol[clust2_idx], WINDOW))
        perc_var_comp_clust1.append(perc_var_comp_clust1_val)
        perc_var_comp_clust2.append(perc_var_comp_clust2_val)
        mvavgs_comp_clust1.append(mvavg_clust1)
        mvavgs_comp_clust2.append(mvavg_clust2)

        perms1 = np.asarray([np.random.permutation(sum(clust1_idx)) for nnn in np.arange(1000)])
        perms2 = np.asarray([np.random.permutation(sum(clust2_idx)) for nnn in np.arange(1000)])
        curr_comp_perm1 = np.asarray([curr_comp_norm[clust1_idx][perm] for perm in perms1])
        curr_comp_perm2 = np.asarray([curr_comp_norm[clust2_idx][perm] for perm in perms2])
        curr_clust1_mvavg_rng_comp = np.apply_along_axis(mvavg, 1, curr_comp_perm1, WINDOW)
        curr_clust2_mvavg_rng_comp = np.apply_along_axis(mvavg, 1, curr_comp_perm2, WINDOW)
        curr_clust1_percvar_rng_comp = np.var(curr_clust1_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm1, axis=1)
        curr_clust2_percvar_rng_comp = np.var(curr_clust2_mvavg_rng_comp, axis=1) / np.var(curr_comp_perm2, axis=1)
        perc_var_comp_clust1_rng.append(curr_clust1_percvar_rng_comp)
        perc_var_comp_clust2_rng.append(curr_clust2_percvar_rng_comp)
        
        windows1 = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(sum(clust1_idx) - WINDOW + 1)])
        windows2 = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(sum(clust2_idx) - WINDOW + 1)])
#        temporal_mov_avg(curr_pol[clust1_idx], curr_comp_norm[clust1_idx], mvavgs_x_clust1[-1], mvavg_clust1, windows1, folder, fileprefixes[i] + "_clust1")
#        temporal_mov_avg(curr_pol[clust2_idx], curr_comp_norm[clust2_idx], mvavgs_x_clust2[-1], mvavg_clust2, windows2, folder, fileprefixes[i] + "_clust2")

    # Test for equal variances of the moving averages and raw values
    perc_var_cell.append(perc_var_cell_val)
    perc_var_nuc.append(perc_var_nuc_val)
    perc_var_cyto.append(perc_var_cyto_val)
    perc_var_mt.append(perc_var_mt_val)
    
    mvavgs_cell.append(mvavg_cell)
    mvavgs_nuc.append(mvavg_nuc)
    mvavgs_cyto.append(mvavg_cyto)
    mvavgs_mt.append(mvavg_mt)
    mvavgs_x.append(mvavg_xvals)
    
    cell_counts.append(len(curr_pol))
    percvar = perc_var_cell_val if wp_iscell[i] else perc_var_nuc_val if wp_isnuc[i] else perc_var_cyto_val
    
    # Uncomment to make the plots (takes 10 mins)
    windows = np.asarray([np.arange(start, start + WINDOW) for start in np.arange(len(curr_pol) - WINDOW + 1)])
#    temporal_mov_avg(curr_pol, curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm, mvavg_xvals,
#         mvavg_cell if wp_iscell[i] else mvavg_nuc if wp_isnuc[i] else mvavg_cyto,
#         windows, folder, fileprefixes[i])
    
alpha_ccd = 0.01
perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = np.array(perc_var_cell),np.array(perc_var_nuc),np.array(perc_var_cyto),np.array(perc_var_mt) # percent variance attributed to cell cycle (mean POI intensities)
perc_var_mt_rng, perc_var_comp_rng = np.array(perc_var_mt_rng), np.array(perc_var_comp_rng) 

# Let's check out which percent variances are greater than the permuted values
perc_var_comp = values_comp(perc_var_cell, perc_var_nuc, perc_var_cyto, wp_iscell, wp_isnuc, wp_iscyto)
perc_var_comp_withbimodal = np.concatenate((perc_var_comp, perc_var_comp_clust1, perc_var_comp_clust2))
perc_var_comp_rng_withbimodal = np.concatenate((perc_var_comp_rng, perc_var_comp_clust1_rng, perc_var_comp_clust2_rng))
ccd_var_comp_rng_wilcoxp_withbimodal = np.apply_along_axis(scipy.stats.wilcoxon, 1, (perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T, None, "wilcox", False, "greater").T[1].T
ccd_var_mt_rng_wilcoxp = np.apply_along_axis(scipy.stats.wilcoxon, 1, (perc_var_mt - perc_var_mt_rng.T).T, None, "wilcox", False, "greater").T[1].T

# randomization tests, try being a bit more stringent, try drawing the cutoff based on microtubules per sample
wp_comp_eq_percvar_adj_withbimodal, wp_comp_pass_eq_percvar_adj_withbimodal = bonf(alpha_ccd, ccd_var_comp_rng_wilcoxp_withbimodal)
wp_comp_gtpass_eq_percvar_adj_withbimodal = wp_comp_pass_eq_percvar_adj_withbimodal & (perc_var_comp_withbimodal > np.median(perc_var_comp_rng_withbimodal, axis=1))

# median differences from random
MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM = 0.05
mean_diff_from_rng_withbimodal = np.mean((perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T, 1)
wp_comp_ccd_difffromrng_withbimodal = mean_diff_from_rng_withbimodal >= MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM

###### Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
alphaa = 0.05
perc_var_mt_valid = perc_var_mt[~np.isinf(perc_var_mt) & ~np.isnan(perc_var_mt)]
percent_var_cutoff = np.median(perc_var_mt_valid)
print(f"{percent_var_cutoff}: cutoff for percent of total variance due to cell cycle")

# separate unimodal from bimodal again
wp_comp_eq_percvar_adj = wp_comp_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
wp_comp_pass_eq_percvar_adj = wp_comp_pass_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
wp_comp_gtpass_eq_percvar_adj = wp_comp_gtpass_eq_percvar_adj_withbimodal[:len(perc_var_comp)]
mean_diff_from_rng = mean_diff_from_rng_withbimodal[:len(perc_var_comp)]
wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng_withbimodal[:len(perc_var_comp)]

def clust_to_wp(clust, clust_idx):
    wp_clust = np.array([False] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust

def clust_to_wp_doub(clust, clust_idx):
    wp_clust = np.array([0.0] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust

wp_comp_ccd_percvar_clust1 = clust_to_wp(perc_var_comp_clust1 >= percent_var_cutoff, wp_isbimodal_fcpadj_pass)
wp_comp_pass_eq_percvar_adj_clust1 = clust_to_wp(wp_comp_pass_eq_percvar_adj_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
mean_diff_from_rng_clust1 = clust_to_wp_doub(mean_diff_from_rng_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
wp_comp_ccd_difffromrng_clust1 = clust_to_wp(wp_comp_ccd_difffromrng_withbimodal[len(perc_var_comp):len(perc_var_comp) + len(perc_var_comp_clust1)], wp_isbimodal_fcpadj_pass)
wp_comp_ccd_percvar_clust2 = clust_to_wp(perc_var_comp_clust2 >= percent_var_cutoff, wp_isbimodal_fcpadj_pass)
wp_comp_pass_eq_percvar_adj_clust2 = clust_to_wp(wp_comp_pass_eq_percvar_adj_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
mean_diff_from_rng_clust2 = clust_to_wp_doub(mean_diff_from_rng_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)
wp_comp_ccd_difffromrng_clust2 = clust_to_wp(wp_comp_ccd_difffromrng_withbimodal[len(perc_var_comp) + len(perc_var_comp_clust1):], wp_isbimodal_fcpadj_pass)

wp_normal_randompercvar_p = np.apply_along_axis(scipy.stats.normaltest, 1, (perc_var_comp - perc_var_comp_rng.T).T).T[1].T
wp_randompercvarnorm_adj, wp_randompercvarnorm_pass = benji_hoch(0.05, wp_normal_randompercvar_p)
print(f"{sum(wp_randompercvarnorm_pass)}: number of genes with randomized percvars that form normal distributions")


wp_comp_ccd_percvar = perc_var_comp >= percent_var_cutoff
print(f"{sum(wp_comp_ccd_percvar)}: # proteins showing CCD variation, comp, percvar")
print(f"{sum(wp_comp_ccd_percvar) / len(wp_comp_ccd_percvar)}: fraction of variable proteins showing CCD variation, comp, percvar")
wp_comp_ccd_percvar_rng = wp_comp_pass_eq_percvar_adj
print(f"{sum(wp_comp_ccd_percvar_rng)}: # proteins showing CCD variation, comp, percvar rng")
print(f"{sum(wp_comp_ccd_percvar_rng) / len(wp_comp_ccd_percvar_rng)}: fraction of variable proteins showing CCD variation, comp, percvar rng")
wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng
print(f"{sum(wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, comp, percvar rng median diff")
print(f"{sum(wp_comp_ccd_difffromrng) / len(wp_comp_ccd_difffromrng)}: fraction of variable proteins showing CCD variation, comp, percvar rng median diff")
wp_comp_ccd_gauss = wp_comp_kruskal_gaussccd_adj <= alphaa
print(f"{sum(wp_comp_ccd_gauss)}: # proteins showing CCD variation, comp, gaussian analysis")
print(f"{sum(wp_comp_ccd_gauss) / len(wp_comp_ccd_gauss)}: fraction of variable proteins showing CCD variation, comp, gaussian analysis")
wp_comp_ccd_gausspercvar = wp_comp_ccd_percvar & wp_comp_ccd_gauss
print(f"{sum(wp_comp_ccd_gausspercvar)}: # proteins showing CCD variation, comp, percvar & gauss")
print(f"{sum(wp_comp_ccd_gausspercvar) / len(wp_comp_ccd_gausspercvar)}: # fraction of proteins showing CCD variation, comp, percvar & gauss")
wp_comp_ccd_gausspercvar_rng = wp_comp_ccd_percvar_rng & wp_comp_ccd_gauss
print(f"{sum(wp_comp_ccd_gausspercvar_rng)}: # proteins showing CCD variation, comp, percvar rng & gauss")
print(f"{sum(wp_comp_ccd_gausspercvar_rng) / len(wp_comp_ccd_gausspercvar_rng)}: # fraction of proteins showing CCD variation, comp, percvar rng & gauss")
wp_comp_ccd_gausspercvar_percvarrng = wp_comp_ccd_percvar_rng & wp_comp_ccd_percvar
print(f"{sum(wp_comp_ccd_gausspercvar_percvarrng)}: # proteins showing CCD variation, comp, percvar & percvar_rng")
print(f"{sum(wp_comp_ccd_gausspercvar_percvarrng) / len(wp_comp_ccd_gausspercvar_percvarrng)}: # fraction of proteins showing CCD variation, comp, percvar rng & gauss")
wp_comp_ccd_gausspercvar_percvarrng_mediandiff = wp_comp_ccd_percvar_rng & wp_comp_ccd_percvar & wp_comp_ccd_difffromrng
print(f"{sum(wp_comp_ccd_gausspercvar_percvarrng_mediandiff)}: # proteins showing CCD variation, comp, percvar > {percent_var_cutoff}, percvar_rng, median_diff")
print(f"{sum(wp_comp_ccd_gausspercvar_percvarrng_mediandiff) / len(wp_comp_ccd_gausspercvar_percvarrng_mediandiff)}: # fraction of proteins showing CCD variation, comp, percvar rng & gauss")

# Address bimodal ones
# 1) The number that are bimodal in one cluster (number of those that are also CCD as unimodal)
# 2) The number that are bimodal in both clusters (number of those that are also CCD as unimodal)
wp_comp_ccd_clust1 = wp_comp_ccd_percvar_clust1 & wp_comp_pass_eq_percvar_adj_clust1 & wp_comp_ccd_difffromrng_clust1
wp_comp_ccd_clust2 = wp_comp_ccd_percvar_clust2 & wp_comp_pass_eq_percvar_adj_clust2 & wp_comp_ccd_difffromrng_clust2
print(f"{sum(wp_isbimodal_fcpadj_pass)}: samples with bimodal antibody intensities")
print(f"{sum(wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2)}: bimodal samples with one CCD cluster ({sum((wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2) & wp_comp_ccd_gausspercvar_percvarrng_mediandiff)}: also CCD unimodally)")
print(f"{sum(wp_comp_ccd_clust1 & wp_comp_ccd_clust2)}: bimodal samples with two CCD clusters ({sum((wp_comp_ccd_clust1 & wp_comp_ccd_clust2) & wp_comp_ccd_gausspercvar_percvarrng_mediandiff)}: also CCD unimodally)")

wp_comp_ccd_use = wp_comp_ccd_gausspercvar_percvarrng_mediandiff # gauss & percvar randomization, like in original manuscript

plt.figure(figsize=(10,10))
plt.scatter(perc_var_comp, mean_diff_from_rng, c=wp_comp_ccd_use, cmap="bwr_r")
plt.vlines(percent_var_cutoff, np.min(mean_diff_from_rng), np.max(mean_diff_from_rng), color="gray")
plt.xlabel("Percent Variance Explained by Cell Cycle")
plt.ylabel("Mean Difference from Random")
plt.savefig("figures/MedianDiffFromRandom.png")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
plt.scatter(mean_diff_from_rng, -np.log10(wp_comp_eq_percvar_adj), c=wp_comp_ccd_use, cmap="bwr_r")
plt.vlines(percent_var_cutoff, np.min(mean_diff_from_rng), np.max(mean_diff_from_rng), color="gray")
plt.xlabel("Mean Difference from Random")
plt.ylabel("-log10 adj p-value from randomization")
plt.savefig("figures/MedianDiffFromRandomVolcano.png")
plt.show()
plt.close()

# Copy profiles to the right place:
# 1) Percent variance
# 2) Percent variance randomization
# 3) Gaussian analysis
# 1 & 2: CCD, 1 & ~2: CCDPercvar, 2 & ~1: CCDPercvarRandomization, 3 & ~1: CCDGauss
ccdbothfolder = f"figures/CCDTemporalMovingAverages{analysis}191211"
ccdpercvarfolder = f"figures/CCDPercvarTemporalMovingAverages{analysis}191211"
ccdpercvarfolderrng = f"figures/CCDPercvarRandomizatinoTemporalMovingAverages{analysis}191211"
ccdgaussfolder = f"figures/CCDGaussTemporalMovingAverages{analysis}191211"
nonccdfolder = f"figures/NonCCDTemporalMovingAverages{analysis}191211"
bimodalccd = f"figures/CCDBimodalTemporalMovingAverages{analysis}191211"
bimodalnonccd = f"figures/NonCCDBimodalTemporalMovingAverages{analysis}191211"
for f in [ccdbothfolder,ccdgaussfolder,ccdpercvarfolder,nonccdfolder,ccdpercvarfolderrng,bimodalccd,bimodalnonccd]:
    if not os.path.exists(f): os.mkdir(f)
for ensg in fileprefixes[wp_comp_ccd_percvar & wp_comp_ccd_percvar_rng & wp_comp_ccd_difffromrng]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(ccdbothfolder, ensg +'_mvavg.png'))
    
for ensg in fileprefixes[~wp_comp_ccd_percvar & ~wp_comp_ccd_difffromrng & ~wp_comp_ccd_percvar_rng]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(nonccdfolder, ensg+'_mvavg.png'))
    
for ensg in fileprefixes[wp_comp_ccd_percvar & (~wp_comp_ccd_percvar_rng | ~wp_comp_ccd_difffromrng)]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(ccdpercvarfolder, ensg+'_mvavg.png'))
    
for ensg in fileprefixes[~wp_comp_ccd_gausspercvar_percvarrng_mediandiff & wp_comp_ccd_gauss]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(ccdgaussfolder, ensg+'_mvavg.png'))
    
for ensg in fileprefixes[~wp_comp_ccd_percvar & wp_comp_ccd_percvar_rng & wp_comp_ccd_difffromrng]:
    shutil.copy(os.path.join(folder, ensg+'_mvavg.png'), os.path.join(ccdpercvarfolderrng, ensg+'_mvavg.png'))

for ensg in fileprefixes[~wp_comp_ccd_use & wp_comp_ccd_clust1]:
    shutil.copy(os.path.join(folder, ensg+'_clust1_mvavg.png'), os.path.join(bimodalccd, ensg+'_clust1_mvavg.png'))
for ensg in fileprefixes[~wp_comp_ccd_use & wp_comp_ccd_clust2]:
    shutil.copy(os.path.join(folder, ensg+'_clust2_mvavg.png'), os.path.join(bimodalccd, ensg+'_clust2_mvavg.png'))
    
for ensg in fileprefixes[~wp_comp_ccd_use & ~wp_comp_ccd_clust1]:
    shutil.copy(os.path.join(folder, ensg+'_clust1_mvavg.png'), os.path.join(bimodalnonccd, ensg+'_clust1_mvavg.png'))
for ensg in fileprefixes[~wp_comp_ccd_use & ~wp_comp_ccd_clust2]:
    shutil.copy(os.path.join(folder, ensg+'_clust2_mvavg.png'), os.path.join(bimodalnonccd, ensg+'_clust2_mvavg.png'))
    
#%% Do the percent variance values match up with what we measured before for the genes?
# Idea: take the perc var values from Devin's analysis and compare them to the ones now
# Execution: run old code and save the data
# Output: plot of percvar vs percvar
alphaa = 0.05

u_well_plates_list = list(u_well_plates)
u_plates_old_idx = np.array([u_well_plates_list.index(wp) for wp in u_well_plates_old if wp in u_well_plates])
old_notfiltered = np.isin(u_well_plates_old, u_well_plates)

plt.figure(figsize=(10,10))
plt.scatter(perc_var_compartment_old[old_notfiltered], perc_var_comp[u_plates_old_idx], c=-np.log10(wp_comp_kruskal_gaussccd_adj[u_plates_old_idx]))
plt.xlabel("percent variance old")
plt.ylabel("percent variance new")
cb = plt.colorbar()
cb.set_label("-log10 FDR for Cell Cycle Dependence")
plt.savefig(f"figures/PercVarAgreement_comp.png")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
plt.scatter(perc_var_comp, mean_mean_comp)
plt.xlabel("percent variance new")
plt.ylabel("mean mean intensity")
plt.savefig(f"figures/PercVarVsMeanMeanIntensity_comp.png")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
plt.scatter(perc_var_comp, -np.log10(wp_comp_eq_percvar_adj))
plt.xlabel("percent variance new")
plt.ylabel("-log10 FDR for CCD")
plt.hlines(-np.log10(alphaa), np.min(perc_var_comp), np.max(perc_var_comp))
plt.savefig(f"figures/PercVarVsLog10FdrCCD_comp.png")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
plt.scatter(mean_mean_comp, -np.log10(wp_comp_eq_percvar_adj))
plt.xlabel("mean mean intensity")
plt.ylabel("-log10 FDR for CCD")
plt.hlines(-np.log10(alphaa), np.min(mean_mean_comp), np.max(mean_mean_comp))
plt.savefig(f"figures/IntensityVsLog10FdrCCD_comp.png")
plt.show()
plt.close()
    
#%% Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
# Idea: create cutoffs for percent variance and 
# Execution: create cutoffs for perc_var and total variance per compartment and for integrated intensity and mean intensity
# Output: Graphs that illustrate the cutoffs (integrated, mean)
# Output: Overlap of the total variance cutoffs with the original filtering done manually

plt.figure(figsize=(10,10))
plt.scatter(gini_comp, perc_var_comp, c=-np.log10(wp_comp_kruskal_gaussccd_adj))
plt.hlines(percent_var_cutoff, np.min(gini_comp), np.max(gini_comp), color="gray")
plt.xlabel("Gini of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Compartment - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CompartmentProteinFractionVariance.png")
plt.show()
plt.close()
   
plt.figure(figsize=(10,10))
plt.scatter(cv_comp, perc_var_comp, c=-np.log10(wp_comp_kruskal_gaussccd_adj))
plt.hlines(percent_var_cutoff, np.min(cv_comp), np.max(cv_comp), color="gray")
plt.xlabel("CV of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Compartment - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CompartmentCVProteinFractionVariance.png")
plt.show()
plt.close()

plt.figure(figsize=(10,10))
plt.scatter(gini_comp, perc_var_comp, c=wp_comp_ccd_use, cmap="bwr_r", alpha=0.5)
plt.hlines(percent_var_cutoff, np.min(gini_comp), np.max(gini_comp), color="gray")
plt.xlabel("Gini of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
#cb = plt.colorbar()
#cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Compartment - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CompartmentProteinFractionVarianceTF.png")
plt.show()
plt.close()

# Output a list of the genes that show variation
name_df = pd.read_csv("input/processed/excel/Fucci_staining_summary_first_plates.csv")
wppp1, ensggg1, abbb1 = list(name_df["well_plate"]), list(name_df["ENSG"]), list(name_df["Antibody"])
name_df2 = pd.read_csv("input/processed/excel/Fucci_staining_review_variation_check.csv")
wppp2, ensggg2, abbb2 = list(name_df2["well_plate"]), list(name_df2["ENSG"]), list(name_df2["Antibody"])
wppp, ensggg, abbb = wppp1 + wppp2, ensggg1 + ensggg2, abbb1 +  abbb2

ensg_dict = dict([(wppp[i], ensggg[i]) for i in range(len(wppp))])
ab_dict = dict([(wppp[i], abbb[i]) for i in range(len(wppp))])
ENSG = np.asarray([ensg_dict[wp] if wp in ensg_dict else "" for wp in well_plate])
antibody = np.asarray([ab_dict[wp] if wp in ab_dict else "" for wp in well_plate])

alpha = 0.05
ccd_comp = wp_comp_ccd_use
nonccd_comp = ~ccd_comp

n_tot_variable = len(u_well_plates)

print(f"{n_tot_variable}: # total proteins showing variation")
print(f"{sum(wp_prev_ccd & wp_comp_ccd_percvar) / len(prev_ccd_ensg)}: fraction of previously annotated CCD genes called CCD (percvar)")
print(f"{sum(wp_prev_ccd & wp_comp_ccd_gausspercvar) / len(prev_ccd_ensg)}: fraction of previously annotated CCD genes called CCD (percvar,gauss)")
print(f"{sum(wp_prev_ccd & wp_comp_ccd_gausspercvar & wp_comp_ccd_percvar_rng) / len(prev_ccd_ensg)}: fraction of previously annotated CCD genes called CCD (percvar,percvarrng)")
print(f"{len([ensg for ensg in wp_ensg[wp_comp_ccd_percvar] if ensg in ensggg1 and ensg in prev_ccd_ensg]) / len(wp_ensg[wp_comp_ccd_percvar])}: fraction of CCD genes that are previously annotated CCD genes (mvavg & gauss)")
print(f"{sum(ccd_comp)}: CCD variable proteins (mvavg & gauss)")
print(f"{sum(wp_comp_ccd_percvar)}: CCD variable proteins (mvavg)")
print(f"{len(prev_ccd_ensg)}: CCD variable proteins (previous)")
print(f"{sum(nonccd_comp)}: non-CCD variable proteins")

### address gene redundancy
wp_ensg_counts = np.array([sum([eeee == ensg for eeee in wp_ensg]) for ensg in wp_ensg])
ensg_is_duplicated = wp_ensg_counts > 1
duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
duplicated_ensg_pairs = [u_well_plates[wp_ensg == ensg] for ensg in duplicated_ensg]
print(f"{sum(wp_comp_ccd_gausspercvar[~ensg_is_duplicated])}: number of CCD proteins (no replicate, percvar&gauss)")
duplicated_ensg_ccd = np.array([sum(wp_comp_ccd_use[wp_ensg == ensg]) for ensg in duplicated_ensg])
print(f"{sum(duplicated_ensg_ccd == 2)}: number of replicated stainings shown to be CCD in both replicates")
print(f"{sum(duplicated_ensg_ccd == 1)}: number of replicated stainings shown to be CCD in just one replicate")
print(f"{sum(duplicated_ensg_ccd == 0)}: number of replicated stainings shown to be non-CCD in both replicate")
pd.DataFrame({
        "ENSG":duplicated_ensg,
        "well_plate_pair":[",".join(pair) for pair in duplicated_ensg_pairs],
        "sum(ccd)":duplicated_ensg_ccd,
        "ccd_pair":[",".join([str(wp_comp_ccd_use[u_well_plates == wp][0]) for wp in pair]) for pair in duplicated_ensg_pairs]
    }).to_csv("output/ReplicatedCellCycleDependentProteins.csv")
###

examples=np.loadtxt("input/processed/manual/ensgexamplesfrompaper.txt", dtype=str)
print(f"{sum(~np.isin(examples,wp_ensg[wp_comp_ccd_use]))}: number of examples missing of {len(examples)}, which includes 1 negative control")
print(examples[~np.isin(examples,wp_ensg[wp_comp_ccd_use])])

plt.figure(figsize=(10,10))
plt.scatter(gini_comp, perc_var_comp, c=-np.log10(wp_comp_kruskal_gaussccd_adj))
for iii in np.nonzero(np.isin(wp_ensg, examples))[0]:
    print(iii)
    plt.annotate(wp_ensg[iii], (gini_comp[iii], perc_var_comp[iii]))
plt.hlines(percent_var_cutoff, np.min(gini_comp), np.max(gini_comp), color="gray")
plt.xlabel("Gini of Protein Expression")
plt.ylabel("Fraction of Variance Due to Cell Cycle")
cb = plt.colorbar()
cb.set_label("FDR for Cell Cycle Dependence")
plt.title("Compartment - Fraction of Variance Due to Cell Cycle")
plt.savefig("figures/CompartmentProteinFractionVarianceAnn.png")
plt.show()
plt.close()

pd.DataFrame({
        "ENSG":wp_ensg[np.nonzero(np.isin(wp_ensg, examples))[0]],
        "gini_comp":gini_comp[np.nonzero(np.isin(wp_ensg, examples))[0]],
        "perc_var_comp":perc_var_comp[np.nonzero(np.isin(wp_ensg, examples))[0]],
        "wp_comp_ccd":wp_comp_ccd_use[np.nonzero(np.isin(wp_ensg, examples))[0]],
        }).to_csv("output/CellCycleExamples.csv")

pd.DataFrame({
    "well_plate" : u_well_plates, 
    "ENSG": wp_ensg,
    "var_comp":var_comp,
    "perc_var_comp":perc_var_comp,
    "pass_percvar":wp_comp_ccd_percvar,
    "wp_comp_kruskal_gaussccd_adj":wp_comp_kruskal_gaussccd_adj,
    "pass_gauss":wp_comp_ccd_gauss,
    "perc_var_rng_wilcoxon_p_adj":wp_comp_eq_percvar_adj,
    "pass_percvar_rng":wp_comp_ccd_percvar_rng,
    "mean_percvar_diff_from_random":mean_diff_from_rng,
    "pass_median_diff":wp_comp_ccd_difffromrng,
    "ccd_comp":ccd_comp,
    "nonccd_comp":nonccd_comp,
    "wp_prev_ccd":wp_prev_ccd,
    }).to_csv("output/CellCycleVariationSummary.csv")
    
pd.DataFrame({"ENSG":wp_ensg[wp_prev_ccd & ~ccd_comp]}).to_csv("output/DianaCCDMissingProteins.csv")
pd.DataFrame({"ENSG":wp_ensg[~wp_prev_ccd & ccd_comp]}).to_csv("output/NewToDianaCCDProteins.csv")


#%% Pickle the results needed later
def np_save_overwriting(fn, arr):
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)

np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp)
np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp)
np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)

pd.DataFrame({"gene": wp_ensg[ccd_comp]}).to_csv("output/picklestxt/ccd_compartment_ensg.txt", index=False)
pd.DataFrame({"gene": wp_ensg[nonccd_comp]}).to_csv("output/picklestxt/nonccd_compartment_ensg.txt", index=False)

