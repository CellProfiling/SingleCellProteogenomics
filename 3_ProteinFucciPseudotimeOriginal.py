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
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
my_df = pd.read_csv("input/raw/nuc_predicted_prob_phases.csv")
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
ab_df = pd.read_excel("input/processed/excel/file_kruskal_by_compartment2.xlsx")
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

###
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
    
# def plot_expression_avg_pseudotime(genelist, outfolder):
#     if not os.path.exists(outfolder): os.mkdir(outfolder)
#     for gene in genelist:
#         nexp = norm_exp_sort[:,list(adata.var_names).index(gene)]
#         df = pd.DataFrame({"fucci_time" : fucci_time_sort, gene : nexp})
#         # plt.scatter(df["fucci_time"], df[gene], label="Normalized Expression")
#         bin_size = 100
#         plt.plot(df["fucci_time"], 
#             df[gene].rolling(bin_size).mean(), 
#             color="blue", 
#             label=f"Moving Average by {bin_size} Cells")
#         plt.fill_between(df["fucci_time"], 
#             df[gene].rolling(bin_size).quantile(0.10),
#             df[gene].rolling(bin_size).quantile(0.90), 
#             color="lightsteelblue", 
#             label="10th & 90th Percentiles")
#         # plt.plot(df["fucci_time"], color="orange", label="Normalized Expression, 10th Percentile")
#         # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() + 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
#         # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() - 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
#         plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
#         plt.ylabel("RNA-Seq Counts, Normalized By Cell",size=36,fontname='Arial')
#         plt.xticks(size=14)
#         plt.yticks(size=14)
#         plt.title(gene,size=36,fontname='Arial')
#         plt.legend(fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f"{outfolder}/{gene}.png")
#         plt.close()

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
    #label stuff
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xlabel('Pseudotime')
    # plt.ylabel(outsuff+' expression')
    plt.ylabel(outsuff.split('_')[0] + ' Protein Expression')
    plt.xticks(size=14)
    plt.yticks(size=14)
    # plt.legend(fontsize=14)
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
model_free_list_all = []
xvals = np.linspace(0,1,num=21)
ccd_pvals = []
not_ccd_pvals = []
perc_green_list = []
perc_red_list = []
neigh = RadiusNeighborsRegressor(radius=0.1)
fdr_df = pd.read_excel("input/processed/excel/file_kruskal_by_compartment2.xlsx")
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

    model_free_list_all.append((well,
            plate_well_to_ab[plate][well][0],
            binned_values[max_loc],
            xvals[max_loc],
            binned_values,
            [perc_var_compartment, perc_var_cell, perc_var_nuc, perc_var_cyto],
            plate_well_to_ab[plate][well][7]))

    if perc_var_compartment>PERCVAR_CUT and curr_fdr<FDR_CUT:
        ccd_coeff_list.append(experiment(outsuff, curr_compartment,
            perc_var_cell, perc_var_nuc, perc_var_cyto,
            max_val_cell, max_val_nuc, max_val_cyto,
            curr_pval, var_compartment, gini_compartment, curr_fdr,
            perc_var_fred, perc_var_fgreen))
        ccd_pvals.append(curr_pval)
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
    if DO_PLOTS:# and any([h in outsuff for h in HIGHLIGHTS]):
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


np.save("output/pickles/u_well_plates.devin.npy", [m[0] for m in model_free_list_all])
np.save("output/pickles/perc_var_compartment.devin.npy", [m[5][0] for m in model_free_list_all])
np.save("output/pickles/perc_var_cell.devin.npy", [m[5][1] for m in model_free_list_all])
np.save("output/pickles/perc_var_nuc.devin.npy", [m[5][2] for m in model_free_list_all])
np.save("output/pickles/perc_var_cyto.devin.npy", [m[5][3] for m in model_free_list_all])
