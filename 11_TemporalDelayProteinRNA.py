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

from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms" # so that the gene-specific results scales match for cross-gene comparisons
print("reading scRNA-Seq data")
biotype_to_use="protein_coding"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=False, biotype_to_use=biotype_to_use)
qc_filtering(adata, do_log_normalize=False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
my_df = pd.read_csv("..\\CellCycleSingleCellRNASeq\\fucci_screen\\nuc_predicted_prob_phases.csv")
print("loaded")

green_fucci = np.asarray(my_df.Intensity_MeanIntensity_CorrResizedGreenFUCCI)
log_green_fucci = np.log10(green_fucci)
red_fucci = np.asarray(my_df.Intensity_MeanIntensity_CorrResizedRedFUCCI)
log_red_fucci = np.log10(red_fucci)
fucci_data = np.column_stack([log_green_fucci,log_red_fucci])

ab_nuc = np.asarray(my_df.Intensity_MeanIntensity_ResizedAb)
ab_cyto = np.asarray(my_df.Mean_ab_Cyto)
ab_cell = np.asarray(my_df.Mean_ab_cell)
well_plate = np.asarray(my_df.well_plate)
u_well_plates = np.unique(well_plate)
ab_objnum = np.asarray(my_df.ObjectNumber)

plt.hist2d(np.log10(green_fucci),np.log10(red_fucci),bins=200)
plt.savefig("figures/FucciPlotProteinIFData.png")
plt.show()
plt.close()

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

#Center data
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

#Convert data to polar
pol_data = cart2pol(centered_data[:,0],centered_data[:,1])
pol_sort_inds = np.argsort(pol_data[1])
pol_sort_rho = pol_data[0][pol_sort_inds]
pol_sort_phi = pol_data[1][pol_sort_inds]
ab_nuc_sort = ab_nuc[pol_sort_inds]
ab_cyto_sort = ab_cyto[pol_sort_inds]
ab_cell_sort = ab_cell[pol_sort_inds]
well_plate_sort = well_plate[pol_sort_inds]
centered_data_sort0 = centered_data[pol_sort_inds,0]
centered_data_sort1 = centered_data[pol_sort_inds,1]
fred_sort = red_fucci[pol_sort_inds]
fgreen_sort = green_fucci[pol_sort_inds]
#rezero to minimum --resoning, cells disappear during mitosis, so we should have the fewest detected cells there
bins = plt.hist(pol_sort_phi,NBINS)
start_phi = bins[1][np.argmin(bins[0])]

#move those points to the other side
more_than_start = np.greater(pol_sort_phi,start_phi)
less_than_start = np.less_equal(pol_sort_phi,start_phi)
pol_sort_rho_reorder = np.concatenate((pol_sort_rho[more_than_start],pol_sort_rho[less_than_start]))
pol_sort_inds_reorder = np.concatenate((pol_sort_inds[more_than_start],pol_sort_inds[less_than_start]))
pol_sort_phi_reorder = np.concatenate((pol_sort_phi[more_than_start],pol_sort_phi[less_than_start]+np.pi*2))
pol_sort_ab_nuc = np.concatenate((ab_nuc_sort[more_than_start],ab_nuc_sort[less_than_start]))
pol_sort_ab_cyto = np.concatenate((ab_cyto_sort[more_than_start],ab_cyto_sort[less_than_start]))
pol_sort_ab_cell = np.concatenate((ab_cell_sort[more_than_start],ab_cell_sort[less_than_start]))
pol_sort_well_plate = np.concatenate((well_plate_sort[more_than_start],well_plate_sort[less_than_start]))
pol_sort_centered_data0 = np.concatenate((centered_data_sort0[more_than_start],centered_data_sort0[less_than_start]))
pol_sort_centered_data1 = np.concatenate((centered_data_sort1[more_than_start],centered_data_sort1[less_than_start]))
pol_sort_fred = np.concatenate((fred_sort[more_than_start],fred_sort[less_than_start]))#+abs(np.min(fred_sort))
pol_sort_fgreen = np.concatenate((fgreen_sort[more_than_start],fgreen_sort[less_than_start]))#+abs(np.min(fgreen_sort))

#shift and re-scale "time"
pol_sort_shift = pol_sort_phi_reorder+np.abs(np.min(pol_sort_phi_reorder))
pol_sort_norm = pol_sort_shift/np.max(pol_sort_shift)
#reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
pol_sort_norm_rev = 1-pol_sort_norm
#stretch time so that each data point is
pol_sort_norm_rev = stretch_time(pol_sort_norm_rev)
#reverse the order of the rho as well
#pol_sort_rho_rev = pol_sort_rho_reorder[::-1]
#pol_sort_inds_rev = pol_sort_inds_reorder*-1+np.max(pol_sort_inds_reorder)


#apply uniform radius (rho) and convert back
#cart_data_ur = pol2cart(np.repeat(R_2,len(centered_data)), pol_data[1])
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

#%% Load the antibody information
ab_df = pd.read_excel("..\\CellCycleSingleCellRNASeq\\fucci_screen\\file_kruskal_by_compartment2.xlsx")
plate_well_to_ab = collections.defaultdict(dict)
plate_dict = collections.defaultdict(dict)
ab_list = ab_df.Antibody
pval_list = ab_df.X_p_value_compartment
gene_list = ab_df.Gene
ensg_list = ab_df.ENSG
plate_well_list = ab_df.well_plate
is_ccd_list = ab_df.correlation
loc_list = ab_df.IF_intensity_var_1
compartment_list = ab_df.Compartments
loc_set = {}
loc_num_pair = [[],[]]
for gene,ab,plate_well,is_ccd,locs,p_val,compartment,ensg in zip(
        gene_list,ab_list,plate_well_list,is_ccd_list,
        loc_list,pval_list,compartment_list,ensg_list):
    plate_num = int(plate_well.split('_')[-1])
    locs_sep = locs.split(';')
    locs_num = []
    for loc in locs_sep:
        if loc not in loc_set:
            loc_set[loc] = len(loc_set)
            loc_num_pair[0].append(loc)
            loc_num_pair[1].append(len(loc_set)-1)
        locs_num.append(loc_set[loc])

    plate_dict[plate_num][plate_well] = [
        gene,ab,is_ccd,locs_num,locs_sep,p_val,compartment,ensg]

plate_well_to_ab = plate_dict

#%%
# Idea: process the well data
# Exec: use Devin's code
# Output: we'll see
PSIN_INIT = [np.nan,1,1,np.nan]
PSIN_BOUNDS = ((0, 1/6, 1/2, 0),(1, 100, 100, 1))
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
DO_PLOTS = False #flag of whether to plot each well and save the plot
TPLOT_MODE = 'psin' #can have values of: 'avg', 'psin'

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
    plt.scatter(curr_pol,curr_ab_norm,c='c')
    plt.plot(x_fit,y_fit,'m')
    #label stuff
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xlabel('pseudo-time')
    plt.ylabel(outsuff+' expression')
    plt.tight_layout()
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.join(os.getcwd(), os.path.dirname(outfile)))
    plt.savefig(outfile)

def temporal_psin(curr_pol,curr_ab_norm,fit_coeffs,outname,outsuff,n_pts=200):
    x_fit = np.linspace(0,1,n_pts)
    plt.close()
    outfile = os.path.join(outname,outsuff+'_psin.pdf')
    if os.path.exists(outfile): return
    fig, ax = plt.subplots(figsize=(10, 10),facecolor='None')
    #plot data
    plt.scatter(curr_pol,curr_ab_norm,c='grey')
    #plot fit
    y_fit = psin_predict(fit_coeffs.x,x_fit)
    plt.plot(x_fit,y_fit,'m')
    #label stuff
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xticks(size=12,fontname='Arial')
    plt.yticks(size=12,fontname='Arial')
    plt.xlabel('pseudo-time',size=20,fontname='Arial')
    plt.ylabel(str.split(outsuff,'_')[0]+' normalized expression', size=20, fontname='Arial')
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
                
missed = []
x_fit = np.linspace(0,1,num=200)
ccd_coeff_list = []
not_ccd_coeff_list = []
model_free_list = []
xvals = np.linspace(0,1,num=21)
ccd_pvals = []
not_ccd_pvals = []
perc_green_list = []
perc_red_list = []
neigh = RadiusNeighborsRegressor(radius=0.1)
fdr_df = pd.read_excel("..\\CellCycleSingleCellRNASeq\\fucci_screen\\file_kruskal_by_compartment2.xlsx")
for well in u_well_plates:
    plt.close('all')
#    well = 'H05_55405991'#GMNN well, used for testing
    curr_well_inds = pol_sort_well_plate==well
    curr_pol = pol_sort_norm_rev[curr_well_inds]
    curr_fred = pol_sort_fred[curr_well_inds]
    curr_fgreen = pol_sort_fgreen[curr_well_inds]

    plate = int(well.split("_")[-1])
    if well in plate_well_to_ab[plate]:
        outsuff = plate_well_to_ab[plate][well][0]+'_'+plate_well_to_ab[plate][well][1]
    else:
        print(well+' was not found.')
        missed.append(well)
        continue

    print(outsuff)

    try:
        curr_fdr = list(fdr_df[fdr_df['well_plate']==well].fdr_pvalue_comp)[0]
    except:
        print(well)
        print(fdr_df['well_plate']==well)
        print(f'could not find fdr for {well}')

    #check if the protein is expressed in the cyto,nuc,or both (cell)
    curr_locs = plate_well_to_ab[plate][well][4]
    curr_pval = plate_well_to_ab[plate][well][5]
    curr_compartment = plate_well_to_ab[plate][well][6].lower()
    curr_ensg = plate_well_to_ab[plate][well][7]
    curr_ab_cell = pol_sort_ab_cell[curr_well_inds]
    curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds]
    curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds]
    curr_ab_cell_norm = curr_ab_cell/np.max(curr_ab_cell)
    curr_ab_nuc_norm = curr_ab_nuc/np.max(curr_ab_nuc)
    curr_ab_cyto_norm = curr_ab_cyto/np.max(curr_ab_cyto)
    curr_fred_norm = curr_fred/np.max(curr_fred)
    curr_fgreen_norm = curr_fgreen/np.max(curr_fgreen)

    #fit psin model
    ###THIS COULD ACTUALLY LITERALLY BE DONE WITH JUST A SMOOTHED MOVING AVERAGE
    #Init all the params
    PSIN_INIT_cell = PSIN_INIT
    PSIN_INIT_nuc = PSIN_INIT
    PSIN_INIT_cyto = PSIN_INIT
    PSIN_INIT_cell[0] = np.max(curr_ab_cell_norm)-np.min(curr_ab_cell_norm)
    PSIN_INIT_nuc[0] = np.max(curr_ab_nuc_norm)-np.min(curr_ab_nuc_norm)
    PSIN_INIT_cyto[0] = np.max(curr_ab_cyto_norm)-np.min(curr_ab_cyto_norm)
    PSIN_INIT_cell[-1] = np.mean(curr_ab_cell_norm)
    PSIN_INIT_nuc[-1] = np.mean(curr_ab_nuc_norm)
    PSIN_INIT_cyto[-1] = np.mean(curr_ab_cyto_norm)

    #Fit the params
    fit_coeffs_cell = least_squares(psin_fit, PSIN_INIT_cell, args=(
        curr_pol, curr_ab_cell_norm),bounds=PSIN_BOUNDS,method='trf')
    fit_coeffs_nuc = least_squares(psin_fit, PSIN_INIT_nuc, args=(
        curr_pol, curr_ab_nuc_norm),bounds=PSIN_BOUNDS,method='trf')
    fit_coeffs_cyto = least_squares(psin_fit, PSIN_INIT_cyto, args=(
        curr_pol, curr_ab_cyto_norm),bounds=PSIN_BOUNDS,method='trf')
    fit_coeffs_fred = least_squares(psin_fit, PSIN_INIT_cyto, args=(
        curr_pol, curr_fred_norm),bounds=PSIN_BOUNDS,method='trf')
    fit_coeffs_fgreen = least_squares(psin_fit, PSIN_INIT_cyto, args=(
        curr_pol, curr_fgreen_norm),bounds=PSIN_BOUNDS,method='trf')
    #Variance calculation
    var_cell = np.var(curr_ab_cell_norm)
    var_nuc = np.var(curr_ab_nuc_norm)
    var_cyto = np.var(curr_ab_cyto_norm)
    var_fred = np.var(curr_fred_norm)
    var_fgreen = np.var(curr_fgreen_norm)
    #Gini coeff calculation
    gini_cell = gini(curr_ab_cell_norm)
    gini_nuc = gini(curr_ab_nuc_norm)
    gini_cyto = gini(curr_ab_cyto_norm)
    #Compute Percent var
    if WINDOW>0:
        perc_var_cell,yval_cell = mvavg_perc_var(curr_ab_cell_norm,WINDOW)
        perc_var_nuc,yval_nuc = mvavg_perc_var(curr_ab_nuc_norm,WINDOW)
        perc_var_cyto,yval_cyto = mvavg_perc_var(curr_ab_cyto_norm,WINDOW)
        perc_var_fred, yval_fred = mvavg_perc_var(curr_fred_norm,WINDOW)
        perc_var_fgreen, yval_fgreen = mvavg_perc_var(curr_fgreen_norm,WINDOW)
        #Get X values
        _, curr_xvals = mvavg_perc_var(curr_pol,WINDOW)
        #get y for the specific compartment we are in
        yvals_compartment = get_val_compartment(
            [yval_cell,yval_nuc,yval_cyto],curr_compartment)
        #get the maximum values
        max_val_cell = np.max(yval_cell)
        max_val_nuc = np.max(yval_nuc)
        max_val_cyto = np.max(yval_cyto)


    else:
        perc_var_cell = 1-np.var(fit_coeffs_cell.fun)/var_cell
        perc_var_nuc = 1-np.var(fit_coeffs_nuc.fun)/var_nuc
        perc_var_cyto = 1-np.var(fit_coeffs_cyto.fun)/var_cyto
        perc_var_fred = 1-np.var(fit_coeffs_fred.fun)/var_fred
        perc_var_fgreen = 1-np.var(fit_coeffs_fgreen.fun)/var_fgreen

        #find the function maximum in our range
        max_val_cell = minimize_scalar(
            lambda x: -psin_predict(fit_coeffs_cell.x,x),
            bounds=[0,1], method='bounded')
        max_val_nuc = minimize_scalar(
            lambda x: -psin_predict(fit_coeffs_nuc.x,x),
            bounds=[0,1], method='bounded')
        max_val_cyto = minimize_scalar(
            lambda x: -psin_predict(fit_coeffs_cyto.x,x),
            bounds=[0,1], method='bounded')

    perc_green_list.append(perc_var_fgreen)
    perc_red_list.append(perc_var_fred)

    #Set values for the current compartment
    #where the protein of interest is expressed
    curr_ab_norm = get_val_compartment(
        [curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm], curr_compartment)
    fit_coeffs = get_val_compartment(
        [fit_coeffs_cell, fit_coeffs_nuc, fit_coeffs_cyto], curr_compartment)
    perc_var_compartment = get_val_compartment(
        [perc_var_cell,perc_var_nuc,perc_var_cyto],curr_compartment)
    var_compartment = get_val_compartment(
        [var_cell,var_nuc,var_cyto],curr_compartment)
    gini_compartment = get_val_compartment(
        [gini_cell,gini_nuc,gini_cyto],curr_compartment)
    max_val_compartment = get_val_compartment(
        [max_val_cell, max_val_nuc, max_val_cyto], curr_compartment)

    #Compute the mean of the data surrounding the fit. It should be near 0
    #If not near 0, throw an error!
    outname = os.path.join("output_devin",plate_well_to_ab[plate][well][2])
    if any(outsuff == i for i in OUTLIER_NAMES):
        fucci_plotting.hist2d_time(green_fucci,red_fucci,
            curr_fgreen,curr_fred,
            curr_pol,outname,outsuff,NBINS)
        #if we have the gmnn antibody, make the pseudo-time plots for
        #fucci-red(CDT1) and fucci-green(GMNN)
        if outsuff == 'GMNN_HPA054597':
            temporal_psin(
                    curr_pol,curr_fred_norm,
                    fit_coeffs_fred,outname,outsuff+'_CDT1_fucci')
            temporal_psin(
                    curr_pol,curr_fgreen_norm,
                    fit_coeffs_fgreen,outname,outsuff+'_GMNN_fucci')

    if perc_var_fred>0.8 and perc_var_fgreen>0.8:
        fucci_plotting.hist2d_time(green_fucci,red_fucci,
            curr_fgreen,curr_fred,
            curr_pol,outname,outsuff,NBINS)

    mean_err = np.mean(fit_coeffs.fun)
    if np.abs(mean_err)>1e-4: print('bad fit!')

    if perc_var_compartment>PERCVAR_CUT and curr_fdr<FDR_CUT:
        ccd_coeff_list.append(experiment(outsuff, curr_compartment,
            perc_var_cell, perc_var_nuc, perc_var_cyto,
            max_val_cell, max_val_nuc, max_val_cyto,
            curr_pval, var_compartment, gini_compartment, curr_fdr,
            perc_var_fred, perc_var_fgreen))
        ccd_pvals.append(curr_pval)

        #compute the model free max
        curr_ab_df = pd.DataFrame({outsuff:curr_ab_norm})
        curr_pol_df = pd.DataFrame({outsuff:curr_pol})
        winsize = np.min([len(curr_ab_norm),30])
        model_free_ab = curr_ab_df.rolling(window=winsize)
        model_free_ab = model_free_ab.mean()
        model_free_ab_np = np.asarray(model_free_ab)
        model_free_pol = curr_pol_df.rolling(window=winsize)
        model_free_pol = model_free_pol.mean()
        model_free_pol_np = np.asarray(model_free_pol)
        max_loc = np.argmax(model_free_ab_np[~np.isnan(model_free_ab_np)])
        max_pol = model_free_pol_np[max_loc]
        binned_values = []
        for xval in xvals:
            if xval==0:
                prev_xval = xval
                continue
            curr_less_inds = curr_pol<xval
            curr_greater_inds = curr_pol>=prev_xval

            binned_values.append(np.median(curr_ab_norm[curr_less_inds&curr_greater_inds]))
            prev_xval = xval
        #fix any nans because they are probably missing data. Use average of surrounding
        binned_values = binned_values/np.nanmax(binned_values)
        binned_values = fix_nans(binned_values)

        max_loc = np.nanargmax(binned_values)
        if np.isnan(xvals[max_loc]):
            print('what')
        print(perc_var_compartment)
        model_free_list.append((well,
            plate_well_to_ab[plate][well][0],
            binned_values[max_loc],
            xvals[max_loc],
            binned_values,
            [perc_var_compartment],
            plate_well_to_ab[plate][well][7]))
        if np.isnan(np.sum(binned_values)):
            print(binned_values)
            print('broken')
    else:
        not_ccd_coeff_list.append(experiment(outsuff, curr_compartment,
            perc_var_cell, perc_var_nuc, perc_var_cyto,
            max_val_cell, max_val_nuc, max_val_cyto,
            curr_pval, var_compartment, gini_compartment, curr_fdr,
            perc_var_fred, perc_var_fgreen))
        not_ccd_pvals.append(curr_pval)

    #visualize the correlation
    if DO_PLOTS:
        outname = os.path.join("output_devin",plate_well_to_ab[plate][well][2])
        if TPLOT_MODE == 'avg':
            temporal_mov_avg(
                    curr_pol,curr_ab_norm,
                    curr_xvals,yvals_compartment,outname,outsuff)
        elif TPLOT_MODE == 'psin':
            temporal_psin(
                    curr_pol,curr_ab_norm,fit_coeffs,outname,outsuff)
        else:
            os.error('Unreckognized TPLOT_MODE.')


#%% Make temporal heatmap and use those peak values to compare to RNA data
#rank by percent variance explained.

def compartment_num(c_name):
    #takes compartment name 'cell','nuc','cyto' and turns it into a number for color plotting
    if c_name=='cell': return 'b'#[255,0,0]
    elif c_name=='nucleus': return 'r'#[0,255,0]
    elif c_name=='cytosol': return 'g'#[0,0,255]
    else: return -1

def perc_var_v_pval(sorted_list):
    perc_var_mat,color_list,varlist,plist = [],[],[],[]
    for item in sorted_list:
        perc_var_mat.append(
            [item.perc_var_cell, item.perc_var_nuc, item.perc_var_cyto])
        color_list.append(compartment_num(item.compartment))
        varlist.append(item.perc_var_compartment)
        try:
            plist.append(np.float(item.pval.strip().split('<')[-1]))
        except:
            plist.append(np.float(item.pval))
    perc_var_mat = np.asarray(perc_var_mat)
    color_list = np.asarray(color_list)
    return plist,varlist,perc_var_mat,color_list

ccd_coeff_list.sort(key=operator.attrgetter('perc_var_compartment'), reverse=True)
not_ccd_coeff_list.sort(key=operator.attrgetter('perc_var_compartment'), reverse=True)

plist,perc_varlist,perc_var_mat,color_list = perc_var_v_pval(ccd_coeff_list)
nc_plist,nc_perc_varlist,nc_perc_var_mat,nc_color_list = perc_var_v_pval(not_ccd_coeff_list)

#Write the outputs to ranked lists
#write_csv_outputs(ccd_coeff_list,outfolder+'/ccd_ranked_outputs_v4.csv')
#write_csv_outputs(not_ccd_coeff_list,outfolder+'/not_ccd_ranked_outputs_v4.csv')


#perc_var_list = [item[2] for item in ccd_coeff_list]
var_list = [item.var_compartment for item in ccd_coeff_list]
gini_list = [item.gini_compartment for item in ccd_coeff_list]
fdr_list = [item.curr_fdr for item in ccd_coeff_list]
nc_var_list = [item.var_compartment for item in not_ccd_coeff_list]
nc_gini_list = [item.gini_compartment for item in not_ccd_coeff_list]
nc_fdr_list = [item.curr_fdr for item in not_ccd_coeff_list]

var_tot = np.concatenate([var_list,nc_var_list])
gini_tot = np.concatenate([gini_list,nc_gini_list])
fdr_tot = np.concatenate([fdr_list,nc_fdr_list])
perc_var_tot = np.concatenate([perc_varlist,nc_perc_varlist])
plist_tot = np.concatenate([plist,nc_plist])


###PLOTTING
# temporal heatmap
# def temporal_heatmap(model_free_list, plate_well_to_ab, loc_set,outfolder):
#PROTEINS SELECTED FOR HIGHLIGHTING
HIGHLIGHTS = ['ORC6','CCNE1','PHLDB1','DPH2']
HIGHLIGHTS_MINOR = [ 'MCM10', 'ZNF32', 'JUN',  
                'DUSP19', 'CCNB1', 'AURKB', 'BUB1B', 'PAPSS1',
                'N6AMT1', 'FLI1']
#HIGHLIGHTS['ORC6','MCM10','ZNF32','JUN','CCNE1','DUSP19',
#   'CCNB1','AURKB','BUB1B','PAPSS1','N6AMT1','PHLDB1','DPH2','FLI1']
#TIMING OF PHASE TRANSITIONS (MANUALLY DETERMINED BY DIANA)
#hours (for the G1/S cutoff)
G1_LEN = 10.833 #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
G1_S_TRANS = 2.625 #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
S_G2_LEN = 11.975 #hours (this should be from the G2/M cutoff above to the end)
#M_LEN = 0.5
#We are excluding Mphase from this analysis
TOT_LEN = G1_LEN+G1_S_TRANS+S_G2_LEN
G1_PROP = G1_LEN/TOT_LEN
G1_S_PROP = G1_S_TRANS/TOT_LEN+G1_PROP
S_G2_PROP = S_G2_LEN/TOT_LEN+G1_S_PROP

#create localization matrix
#sort the coefficients based on the max value point
loc_mat = np.zeros([len(model_free_list),len(loc_set)])
for i,response in enumerate(model_free_list):
    curr_well = response[0]
    plate = int(curr_well.split("_")[-1])
    print(curr_well)
    print(plate_well_to_ab[plate][curr_well])
    curr_locs = plate_well_to_ab[plate][curr_well][3]
    loc_mat[i,curr_locs] = 1

model_free_list.sort(key=operator.itemgetter(3))
prob_interact = np.zeros([len(model_free_list),len(model_free_list)])
sum_rows = []
prot_list = []
sorted_gene_array = []
sorted_ensg_array = []
sorted_maxpol_array = []
for i,response in enumerate(model_free_list):
    sorted_gene_array.append(response[4]) # this is the expression values
    sorted_maxpol_array.append(response[3])
    sorted_ensg_array.append(response[6])
    curr_well = response[0]
    #set all other rows that share a loc to 1
    curr_locs = loc_mat[i,:]
    match_rows = np.greater(np.sum(loc_mat*curr_locs,axis=1),0)*response[3]
    sum_rows.append(np.sum(match_rows))
    prob_interact[i,:] = match_rows
    prot_list.append(response[1])

fig, ax = plt.subplots(figsize=(10, 10))
sc = ax.imshow(sorted_gene_array, interpolation='nearest')

#Arange the x ticks
xtick_labels = [str(np.around(x * TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)]
xtick_labels = xtick_labels#+['G1/S','S/G2']
xphase_labels = ['G1/S','S/G2']
my_xticks = np.arange(-.5, 20, 2)
num_ticks = 20

phase_trans = np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])

ax.set_xticks(my_xticks,minor=True)
ax.set_xticklabels(xtick_labels,minor=True)
ax.set_xticks(phase_trans, minor=False)
ax.set_xticklabels(xphase_labels, minor=False)
ax.tick_params(length=12)
plt.xticks(size=12,fontname='Arial')
plt.yticks(size=10,fontname='Arial')

#Do the y ticks
ytick_locs = [prot_list.index(p) for p in HIGHLIGHTS]
ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
ax.set_yticks(ytick_locs,minor=False)
ax.set_yticklabels(HIGHLIGHTS,minor=False)
ax.set_yticks(ytick_locs_minor,minor=True)
ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)
ax.tick_params(direction='out', length=12, width=2, colors='k', 
                axis='x',which='major')
ax.tick_params(direction='out', length=56, width=1, colors='k', 
                axis='y',which='major')
ax.set_aspect('auto')
plt.xlabel('Division Cycle, hrs',size=20,fontname='Arial')
plt.ylabel('Gene',size=20,fontname='Arial')

divider1 = make_axes_locatable(ax)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(sc,cax = cax1)
cbar.set_label('Relative expression',fontname='Arial',size=20)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(os.path.join("output_devin",'sorted_heatmap21_sw30_take4.pdf'), transparent=True)
plt.show()
# return fig,ax,prot_list

#heatmap
#pval vs fracvar
#fucci_plotting.p_v_fracvar(plist_tot,perc_var_tot,fdr_tot,outfolder)
#gini vs fracvar
#fucci_plotting.gini_v_fracvar(gini_tot,perc_var_tot,fdr_tot,outfolder)
#variance vs fracvar
# fucci_plotting.var_v_fracvar(var_tot,perc_var_tot,fdr_tot,outfolder)
#perc var list
# fucci_plotting.fracvar_compartments(perc_varlist,nc_perc_varlist,outfolder)


#%%
# Idea: calculate the peak RNA expression and compare to the peak protein expression for each gene
# Execution: compare distribution of differences between peaks; grainger test, too
# Output: plot of dist of differences; grainger test results
def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')

# Get the peak RNA expression polar locations
expression_data = adata.X # log normalized
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, 100)
max_moving_avg_loc = np.argmax(moving_averages, 0)
max_moving_avg_pol = np.take(fucci_time_sort, max_moving_avg_loc)

# Get the gene names
gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
var_list = list(adata.var_names)
var_set = set(adata.var_names)
ensg_rna_loc = [var_list.index(gene_id) for gene_id in gene_info["gene_id"] if gene_id in var_set]
name_rna_list = list(np.take(np.array(gene_info["name"]), ensg_rna_loc))

# Are we getting the expected majority overlap in gene sets (protein & rna)
allccd_transcript_regulated = np.array(pd.read_csv("output/allccd_transcript_regulated.csv")["gene"])
allccd_transcript_regulated_names = set(ccd_gene_names(allccd_transcript_regulated))
prot_genes = set(prot_list)
rna_genes = set(name_rna_list)
print(f"length of prot genes: {len(prot_genes)}")
print(f"length of RNA-seq genes: {len(rna_genes)}")
print(f"length of CCD RNA genes: {len(allccd_transcript_regulated_names)}")
print(f"length of intersection betweeen CCD prot and RNA: {len(prot_genes.intersection(rna_genes))}")
print(f"length of intersection betweeen CCD prot and CCD RNA: {len(prot_genes.intersection(allccd_transcript_regulated_names))}")

# Sort them to the protein arrays
# prot_list # this one is the gene name (protein name) for the protein list
# sorted_ensg_array # this one is the sorted ensg for the proteins
# sorted_maxpol_array # this one is the protein max pol
name_prot_loc = [name_rna_list.index(name) for name in prot_list if name in name_rna_list and name in allccd_transcript_regulated_names]
name_prot_list = np.take(name_rna_list, name_prot_loc)
max_rna_avg_prot_pol = np.take(max_moving_avg_pol, name_prot_loc)
prot_list_filter_loc = np.isin(prot_list, name_prot_list)
prot_maxpol_filter_array = np.array(sorted_maxpol_array)[prot_list_filter_loc]
diff_max_pol = prot_maxpol_filter_array - max_rna_avg_prot_pol

plt.hist(diff_max_pol * TOT_LEN)
plt.xlabel("Delay in peak protein expression from peak RNA expression, hrs")
plt.ylabel("Count of CCD Proteins")
plt.tight_layout()
plt.savefig("figures/DelayPeakProteinRNA.png")
plt.show()
plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hist(prot_maxpol_filter_array * TOT_LEN, alpha=0.5, label="Peak Protein Expression Time, hrs")
ax1.hist(max_rna_avg_prot_pol * TOT_LEN, alpha=0.5, label="Peak RNA Expression Time, hrs")
plt.legend(loc="upper left")
plt.xlabel("Psueodtime")
plt.ylabel("Count of Cell Cycle Genes")
plt.tight_layout()
plt.savefig(f"figures/DelayPeakProteinRNA_separate.png")
plt.show()
plt.close()

mmmm = np.concatenate((prot_maxpol_filter_array * TOT_LEN, max_rna_avg_prot_pol * TOT_LEN))
cccc = (["Protein"] * len(prot_maxpol_filter_array))
cccc.extend(["RNA"] * len(max_rna_avg_prot_pol))
moddf = pd.DataFrame({"time": mmmm, "category" : cccc})
boxplot = moddf.boxplot("time", by="category", figsize=(12, 8), showfliers=True)
boxplot.set_xlabel("", size=36,fontname='Arial')
boxplot.set_ylabel("Peak Expression, hrs", size=36,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=16)
plt.title("")
plt.savefig("figures/DelayPeakProteinRNA_boxplot.png")
plt.show()
plt.close()

print(f"Median RNA expression time for CCD proteins: {np.median(max_rna_avg_prot_pol)}")
print(f"Median protein expression time for CCD proteins: {np.median(prot_maxpol_filter_array)}")
t, p = scipy.stats.kruskal(max_rna_avg_prot_pol, prot_maxpol_filter_array)
print(f"One-sided kruskal for median protein expression time higher than median RNA expression time: {2*p}")
t, p = scipy.stats.ttest_1samp(diff_max_pol, 0)
print(f"One-sided, one-sample t-test for mean delay in protein expression larger than zero: {2*p}")

# diff_max_pol_regev = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_regev_filtered]
# mean_dianaccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
# diff_max_pol_dianaccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in ccd_filtered]
# mean_diananonccd = [np.mean(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]
# variances_diananonccd = [np.std(expression_data[:,geneidx]) for geneidx in range(len(g1_exp[0,:])) if adata.var_names[geneidx] in nonccd_filtered]

# def weights(vals):
#     '''normalizes all histogram bins to sum to 1'''
#     return np.ones_like(vals)/float(len(vals))

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# bins=np.histogram(np.hstack((variances, variances_regev, variances_dianaccd, variances_diananonccd)), bins=40)[1] #get the bin edges
# ax1.hist(variances_regev, bins=bins, weights=weights(variances_regev), 
#     label="Regev CCD Genes")
# ax1.hist(variances_dianaccd, bins=bins, weights=weights(variances_dianaccd),
#     label="Fucci CCD Genes")
# ax1.hist(variances_diananonccd, bins=bins, weights=weights(variances_diananonccd), 
#     label="Fucci Non-CCD Genes")
# ax1.hist(variances, bins=bins, weights=weights(variances), 
#     label="All Genes")
# plt.legend(loc="upper right")
# plt.xlabel("Stdev Expression")
# plt.ylabel("Count, Normalized to 1")
# plt.tight_layout()
# plt.savefig(f"figures/stdev_expression_hist_{biotype_to_use}.png")
# plt.show()
# plt.close()


#%% Sanity checks 
# double check that the names line up
prot_names = np.array(prot_list)[prot_list_filter_loc]
rna_names = name_prot_list
print(f"The name arrays are the same: {all(prot_names == rna_names)}")

# What are the smallest, largest, and median genes and what do they look like?
smallest = np.argmin(diff_max_pol)
smallest_gene = name_prot_list[smallest]
median = np.argsort(diff_max_pol)[len(diff_max_pol)//2]
median_gene = name_prot_list[median]
largest = np.argmax(diff_max_pol)
largest_gene = name_prot_list[largest]
print(f"smallest delay {diff_max_pol[smallest] * TOT_LEN} hr for {name_prot_list[smallest]}")
print(f"median delay {diff_max_pol[median] * TOT_LEN} hr for {name_prot_list[median]}")
print(f"largest delay {diff_max_pol[largest] * TOT_LEN} hr for {name_prot_list[largest]}")

plt.rcParams['figure.figsize'] = (10, 10)
sorted_gene_array_array = np.array(sorted_gene_array).transpose()
def plot_avg_rna_and_prot(namelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for name in namelist:
        bin_size = 100
        mvgavg = moving_averages[:,name_rna_list.index(name)]
        mvgmax = mvgavg.max()
        plt.plot(
            fucci_time_sort[:-(bin_size-1)] * TOT_LEN, 
            mvgavg / mvgmax,  
            color="blue", 
            label=f"RNA Moving Average by {bin_size} Cells")
        prot_exp = sorted_gene_array_array[:,prot_list.index(name)]
        plt.plot(
            xvals[:-1] * TOT_LEN, 
            prot_exp, 
            color="red", 
            label="Protein Expression")
        plt.xlabel("Cell Division Time, hrs",size=36,fontname='Arial')
        plt.ylabel("Expression, Normalized by Cell",size=36,fontname='Arial')
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.title(name,size=36,fontname='Arial')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{name}.png")
        plt.close()

plot_avg_rna_and_prot(name_prot_list, "figures/RNAProteinCCDAvgs")


#%%
