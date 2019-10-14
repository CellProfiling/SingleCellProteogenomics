#%% Imports
from imports import *
import numpy as np
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
from sklearn.mixture import GaussianMixture
import collections
import fucci_plotting
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists, ccd_gene_names

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms" # so that the gene-specific results scales match for cross-gene comparisons
print("Reading scRNA-Seq data")
biotype_to_use="protein_coding"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=False, biotype_to_use=biotype_to_use)
qc_filtering(adata, do_log_normalize=False)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
#       use the mean intensity and integrated intensity for different cutoffs (ab and microtubules)
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
my_df1 = pd.read_csv("..\\CellCycleSingleCellRNASeq\\fucci_screen\\nuc_predicted_prob_phases_mt_all_firstbatch_plates.csv")
my_df2 = pd.read_csv("..\\CellCycleSingleCellRNASeq\\fucci_screen\\nuc_predicted_prob_phases_190909.csv")
my_df = pd.concat((my_df1, my_df2), sort=True)
print("loaded")

# Sample information (FucciWellPlateGene is still missing some information)
plate = np.asarray(my_df.plate)
u_plate = np.unique(plate)
well_plate = np.asarray(my_df.well_plate)
imgnb = np.asarray(my_df.ImageNumber)
u_well_plates = np.unique(well_plate)
ab_objnum = np.asarray(my_df.ObjectNumber)
name_df = pd.read_csv("input\\Fucci_staining_summary_first_plates.csv")
wppp, ensggg, abbb, rrrr = list(name_df["well_plate"]), list(name_df["ENSG"]), list(name_df["Antibody"]), list(name_df["Results_final_update"])
# name_df = pd.read_csv("input\\FucciWellPlateGene.csv")
# wppp, ensggg, abbb = list(name_df["well_plate"]), list(name_df["ENSG"]), list(name_df["Antibody"])
ensg_dict = dict([(wppp[i], ensggg[i]) for i in range(len(wppp))])
ab_dict = dict([(wppp[i], abbb[i]) for i in range(len(wppp))])
result_dict = dict([(wppp[i], rrrr[i]) for i in range(len(wppp))])
ENSG = np.asarray([ensg_dict[wp] if wp in result_dict else "" for wp in well_plate])
antibody = np.asarray([ab_dict[wp] if wp in result_dict else "" for wp in well_plate])
result = np.asarray([result_dict[wp] if wp in result_dict else "" for wp in well_plate])
neg="neg"
print(f"{len([s for s in result if s.lower().startswith(neg)])}: number of negative stainings")

# Let's take a look at cell area to see if that's a good cutoff
# Looks fine. They already filtered on the area of the object; the area of the cell looks relatively normal.
# area_object = np.asarray(my_df.AreaShape_Area)
# bins = plt.hist(np.log10(area_object), bins=100)
# plt.vlines(np.mean(np.log10(area_object)), 0, np.max(bins[0]))
# area_object_mean = np.mean(np.log10(area_object))
# upper_area_cutoff = area_object_mean + 4 * np.std(np.log10(area_object))
# lower_area_cutoff = area_object_mean - 4 * np.std(np.log10(area_object))
# plt.vlines(upper_area_cutoff, 0, np.max(bins[0]))
# plt.vlines(lower_area_cutoff, 0, np.max(bins[0]))
# plt.show()
# plt.
# area_cell = np.asarray(my_df.Area_cell)
# area_cyto = np.asarray(my_df.Area_cyto)
# bins = plt.hist(np.log10(area_cell), bins=100)
# plt.vlines(np.mean(np.log10(area_cell)), 0, np.max(bins[0]))
# area_cell_mean = np.mean(np.log10(area_cell))
# upper_area_cutoff = area_cell_mean + 4 * np.std(np.log10(area_cell))
# lower_area_cutoff = area_cell_mean - 4 * np.std(np.log10(area_cell))
# plt.vlines(upper_area_cutoff, 0, np.max(bins[0]))
# plt.vlines(lower_area_cutoff, 0, np.max(bins[0]))
# plt.show()
# plt.hist(np.log10(area_cyto), bins=100)
# plt.show()

# Fucci data (mean intensity)
# The microscope used for the new plates has a wider dynamic range, and the gain is set for the highest expressing well on the whole plate.
# It appears this also makes the 
green_offset = dict([
    (6716,0),
    (6718,0.002),
    (6719,0),
    (6720,0.002),
    (6721,0.004),
    (6722,0.004),
    (6724,0.004),
    (6725,0.004),
    (6731,0),
    (6734,0),
    (6735,0.004),
    (6736,0),
    (6745,0.002)])
green_fucci = np.asarray(my_df.Intensity_MeanIntensity_CorrResizedGreenFUCCI)
green_fucci_new = [green_fucci[i] + green_offset[plate[i]] if plate[i] in green_offset else green_fucci[i] for i in range(len(green_fucci))]
red_fucci = np.asarray(my_df.Intensity_MeanIntensity_CorrResizedRedFUCCI)
log_red_fucci = np.log10(red_fucci)
fucci_data = np.column_stack([log_green_fucci,log_red_fucci])
plt.hist2d(np.log10(green_fucci_new),np.log10(red_fucci),bins=200)
plt.savefig("figures/FucciPlotProteinIFData_unfiltered.png")
plt.show()
plt.close()
for p in u_plate:
    plt.hist2d(np.log10(green_fucci[plate == p]),np.log10(red_fucci[plate == p]),bins=200)
    plt.savefig(f"figures/FucciPlotProteinIFData_filterNegativeStains{p}.png")
    plt.title(f"FUCCI mean intensities, plate {p}")
    plt.show()
    plt.close()

# Fucci data (integrated intensity)
green_fucci_int = np.asarray(my_df.Intensity_IntegratedIntensity_CorrResizedGreenFUCCI)
log_green_fucci_int = np.log10(green_fucci_int)
red_fucci_int = np.asarray(my_df.Intensity_IntegratedIntensity_CorrResizedRedFUCCI)
log_red_fucci_int = np.log10(red_fucci_int)
fucci_data_int = np.column_stack([log_green_fucci_int,log_red_fucci_int])
plt.hist2d(np.log10(green_fucci_int),np.log10(red_fucci_int),bins=200)
plt.savefig("figures/FucciPlotProteinIFData_unfiltered.png")
plt.show()
plt.close()

# Antibody data (mean intensity)
ab_nuc = np.asarray(my_df.Intensity_MeanIntensity_ResizedAb)
ab_cyto = np.asarray(my_df.Mean_ab_Cyto)
ab_cell = np.asarray(my_df.Mean_ab_cell)
mt_cell = np.asarray(my_df.Mean_mt_cell)
bins = plt.hist(np.log10(ab_cell), bins=100)
plt.vlines(np.mean(np.log10(ab_cell)), 0, np.max(bins[0]))
ab_cell_mean = np.mean(np.log10(ab_cell))
upper_ab_cell_cutoff = ab_cell_mean + 3 * np.std(np.log10(ab_cell))
lower_ab_cell_cutoff = ab_cell_mean - 3 * np.std(np.log10(ab_cell))
lower_ab_cell_heuristic_cutoff = -2.4
plt.vlines(upper_ab_cell_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_ab_cell_heuristic_cutoff, 0, np.max(bins[0]))
plt.title("All Mean Intensities in Antibody Channel including Neg Controls")
plt.show()

# Negative controls (mean intensity)
neg_control = np.asarray([s.startswith("H12") for s in well_plate])
neg_control_ab_cell = ab_cell[neg_control]
bins = plt.hist(np.log10(neg_control_ab_cell), bins=100)
plt.vlines(np.mean(np.log10(neg_control_ab_cell)), 0, np.max(bins[0]))
neg_control_ab_cell_mean = np.mean(np.log10(neg_control_ab_cell))
upper_neg_control_ab_cell_cutoff = neg_control_ab_cell_mean + 3 * np.std(np.log10(neg_control_ab_cell))
lower_neg_control_ab_cell_cutoff = neg_control_ab_cell_mean - 3 * np.std(np.log10(neg_control_ab_cell))
plt.vlines(upper_neg_control_ab_cell_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_neg_control_ab_cell_cutoff, 0, np.max(bins[0]))
plt.title("All Mean Intensities in Antibody Channel of Neg Controls")
plt.show()

# Antibody data (integrated intensity)
ab_nuc_int = np.asarray(my_df.Intensity_IntegratedIntensity_ResizedAb)
ab_cyto_int = np.asarray(my_df.Integrated_ab_cyto)
ab_cell_int = np.asarray(my_df.Integrated_ab_cell)
mt_cell_int = np.asarray(my_df.Integrated_mt_cell)
bins = plt.hist(np.log10(ab_cell_int), bins=100)
plt.vlines(np.mean(np.log10(ab_cell_int)), 0, np.max(bins[0]))
ab_cell_int_mean = np.mean(np.log10(ab_cell_int))
upper_ab_cell_int_cutoff = ab_cell_int_mean + 3 * np.std(np.log10(ab_cell_int))
lower_ab_cell_int_cutoff = ab_cell_int_mean - 3 * np.std(np.log10(ab_cell_int))
plt.vlines(upper_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.title("All Integrated Intensities in Antibody Channel including Neg Controls")
plt.show()

# Negative controls (integrated intensity)
neg_control_ab_cell_int = ab_cell_int[neg_control]
bins = plt.hist(np.log10(neg_control_ab_cell_int), bins=100)
plt.vlines(np.mean(np.log10(neg_control_ab_cell_int)), 0, np.max(bins[0]))
neg_control_ab_cell_int_mean = np.mean(np.log10(neg_control_ab_cell_int))
upper_neg_control_ab_cell_int_cutoff = neg_control_ab_cell_int_mean + 1 * np.std(np.log10(neg_control_ab_cell_int))
lower_neg_control_ab_cell_int_cutoff = neg_control_ab_cell_int_mean - 1 * np.std(np.log10(neg_control_ab_cell_int))
plt.vlines(upper_neg_control_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_neg_control_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.title("All Mean Intensities in Antibody Channel of Neg Controls")
plt.show()

# Negative staining (mean intensity)
ab_cell_neg = ab_cell[[s.lower().startswith("neg") for s in result]]
bins = plt.hist(np.log10(ab_cell_neg), bins=100)
plt.title("Negative staining (mean intensity)")
plt.show()
ab_cell_int_neg = ab_cell_int[[s.lower().startswith("neg") for s in result]]
bins = plt.hist(np.log10(ab_cell_int_neg), bins=100)
plt.title("Negative staining (integrated intensity)")
plt.show()

# Filter out negative staining
# pass_filter_ab_cell_mean = (np.log10(ab_cell) > lower_ab_cell_heuristic_cutoff) & (log_green_fucci > lower_ab_cell_heuristic_cutoff) & (log_red_fucci > lower_ab_cell_heuristic_cutoff)
pass_filter_ab_cell_mean = (np.log10(ab_cell_int) > upper_neg_control_ab_cell_int_cutoff) & (~neg_control)
green_fucci = green_fucci[pass_filter_ab_cell_mean]
log_green_fucci = log_green_fucci[pass_filter_ab_cell_mean]
red_fucci = red_fucci[pass_filter_ab_cell_mean]
log_red_fucci = log_red_fucci[pass_filter_ab_cell_mean]
fucci_data = fucci_data[pass_filter_ab_cell_mean]
ab_nuc = ab_nuc[pass_filter_ab_cell_mean]
ab_cyto = ab_cyto[pass_filter_ab_cell_mean]
ab_cell = ab_cell[pass_filter_ab_cell_mean]
mt_cell = mt_cell[pass_filter_ab_cell_mean]
ab_nuc_int = ab_nuc_int[pass_filter_ab_cell_mean]
ab_cyto_int = ab_cyto_int[pass_filter_ab_cell_mean]
ab_cell_int = ab_cell_int[pass_filter_ab_cell_mean]
mt_cell_int = mt_cell_int[pass_filter_ab_cell_mean]
ENSG = ENSG[pass_filter_ab_cell_mean]
antibody = antibody[pass_filter_ab_cell_mean]
result = result[pass_filter_ab_cell_mean]
neg="neg"
print(f"{len([s for s in result if s.lower().startswith(neg)])}: number of negative stainings")

# New antibody intensity plots with cutoffs
bins = plt.hist(np.log10(ab_cell), bins=100)
plt.vlines(np.mean(np.log10(ab_cell)), 0, np.max(bins[0]))
ab_cell_mean = np.mean(np.log10(ab_cell))
upper_ab_cell_cutoff = ab_cell_mean + 3 * np.std(np.log10(ab_cell))
lower_ab_cell_cutoff = ab_cell_mean - 3 * np.std(np.log10(ab_cell))
plt.vlines(upper_ab_cell_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_ab_cell_heuristic_cutoff, 0, np.max(bins[0]))
plt.show()
bins = plt.hist(np.log10(ab_cell_int), bins=100)
plt.vlines(np.mean(np.log10(ab_cell_int)), 0, np.max(bins[0]))
ab_cell_int_mean = np.mean(np.log10(ab_cell_int))
upper_ab_cell_int_cutoff = ab_cell_int_mean + 3 * np.std(np.log10(ab_cell_int))
lower_ab_cell_int_cutoff = ab_cell_int_mean - 3 * np.std(np.log10(ab_cell_int))
plt.vlines(upper_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.vlines(lower_ab_cell_int_cutoff, 0, np.max(bins[0]))
plt.show()
ab_cell_neg = ab_cell[[s.lower().startswith("neg") for s in result]]
bins = plt.hist(np.log10(ab_cell_neg), bins=100)
plt.title("Negative staining (mean intensity)")
plt.show()
ab_cell_int_neg = ab_cell_int[[s.lower().startswith("neg") for s in result]]
bins = plt.hist(np.log10(ab_cell_int_neg), bins=100)
plt.title("Negative staining (integrated intensity)")
plt.show()

# Plot FUCCI intensities
plt.hist2d(np.log10(green_fucci),np.log10(red_fucci),bins=200)
plt.savefig("figures/FucciPlotProteinIFData_filterNegativeStains.png")
plt.show()
plt.close()

#%%
plt.hist2d(np.log10(ab_cell), np.log10(ab_cell_int), bins=200)

#%%a


#%%
# Idea: Gaussian clustering per plate to identify G1/S/G2 and do kruskal test for variance
# Exec: sklearn.mixture.GaussianMixture & scipy.stats.kruskal
# Output: FDR for cell cycle variation per well

# out_of_focus_well_plate_imagenb = pd.read_csv("input/outoffocusimages.txt") # not doing this for now; use all the data, since there are surely more that are out of focus
for plate in u_plate:
    gaussian = GaussianMixture(n_components=3, random_state=1, max_iter=500)
    use_idx = my_df.plate == plate
    gaussian.fit_predict(log_green_fucci[my_df.plate == plate], log_red_fucci[my_df.plate == plate])

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
ab_nuc_sort_int, ab_cyto_sort_int, ab_cell_sort_int, mt_cell_sort_int = pol_sort(pol_sort_inds,ab_nuc_int,ab_cyto_int,ab_cell_int,mt_cell_int)
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
pol_sort_ab_nuc_int, pol_sort_ab_cyto_int, pol_sort_ab_cell_int, pol_sort_mt_cell_int = pol_reord(ab_nuc_sort_int), pol_reord(ab_cyto_sort_int), pol_reord(ab_cell_sort_int), pol_reord(mt_cell_sort_int)
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
fucci_hist2d(centered_data,cart_data_ur,start_pt)

#%% Load the antibody information
# Thinking of redoing this; this is pretty much just statistics for each sample
# well_plate unique
# decides ccd/nonccd
# some lookups from that File_Devin file...?
# # some statistics based on kruskal for G1/S/G2
# ab_df = pd.read_excel("..\\CellCycleSingleCellRNASeq\\fucci_screen\\file_kruskal_by_compartment2.xlsx")
# plate_well_to_ab = collections.defaultdict(dict)
# plate_dict = collections.defaultdict(dict)
# ab_list = ab_df.Antibody
# pval_list = ab_df.X_p_value_compartment
# gene_list = ab_df.Gene
# ensg_list = ab_df.ENSG
# plate_well_list = ab_df.well_plate
# is_ccd_list = ab_df.correlation
# loc_list = ab_df.IF_intensity_var_1
# compartment_list = ab_df.Compartments
# loc_set = {}
# loc_num_pair = [[],[]]
# for gene,ab,plate_well,is_ccd,locs,p_val,compartment,ensg in zip(
#         gene_list,ab_list,plate_well_list,is_ccd_list,
#         loc_list,pval_list,compartment_list,ensg_list):
#     plate_num = int(plate_well.split('_')[-1])
#     locs_sep = locs.split(';')
#     locs_num = []
#     for loc in locs_sep:
#         if loc not in loc_set:
#             loc_set[loc] = len(loc_set)
#             loc_num_pair[0].append(loc)
#             loc_num_pair[1].append(len(loc_set)-1)
#         locs_num.append(loc_set[loc])

#     plate_dict[plate_num][plate_well] = [
#         gene,ab,is_ccd,locs_num,locs_sep,p_val,compartment,ensg]

# plate_well_to_ab = plate_dict

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
var_cell, var_nuc, var_cyto, var_mt = [],[],[],[] # mean intensity variances per antibody
var_cell_int, var_nuc_int, var_cyto_int, var_mt_int = [],[],[],[] # integrated intensity variances per antibody

# These are tuples of (perc_var, y_val_avgs) where perc_var is the value described in the comments below, and the y_val_avgs are the moving averages
perc_var_fred, perc_var_fgreen = [],[] # percent variance attributed to cell cycle (FUCCI colors)
perc_var_cell, perc_var_nuc, perc_var_cyto, perc_var_mt = [],[],[],[] # percent variance attributed to cell cycle (mean POI intensities)
perc_var_cell_int, perc_var_nuc_int, perc_var_cyto_int, perc_var_mt_int = [],[],[],[] # percent variance attributed to cell cycle (integrated POI intensities)
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
    curr_ab_cell_int, curr_ab_nuc_int, curr_ab_cyto_int, curr_mt_cell_int = pol_sort_ab_cell_int[curr_well_inds], pol_sort_ab_nuc_int[curr_well_inds],pol_sort_ab_cyto_int[curr_well_inds], pol_sort_mt_cell_int[curr_well_inds]
    curr_ab_cell_norm_int, curr_ab_nuc_norm_int, curr_ab_cyto_norm_int, curr_mt_cell_norm_int = curr_ab_cell_int/np.max(curr_ab_cell_int), curr_ab_nuc_int/np.max(curr_ab_nuc_int), curr_ab_cyto_int/np.max(curr_ab_cyto_int), curr_mt_cell_int/np.max(curr_mt_cell_int)

    # Variance calculation FUCCI
    var_fred.append(np.var(curr_fred_norm))
    var_fgreen.append(np.var(curr_fgreen_norm))
    # Variance calculation -- mean intensities
    var_cell.append(np.var(curr_ab_cell_norm))
    var_nuc.append(np.var(curr_ab_nuc_norm))
    var_cyto.append(np.var(curr_ab_cyto_norm))
    var_mt.append(np.var(curr_mt_cell_norm))
    # Variance calculation -- integrated intensities
    var_cell_int.append(np.var(curr_ab_cell_norm_int))
    var_nuc_int.append(np.var(curr_ab_nuc_norm_int))
    var_cyto_int.append(np.var(curr_ab_cyto_norm_int))
    var_mt_int.append(np.var(curr_mt_cell_norm_int))

    # Compute Percent var, # if WINDOW>0: # always true for this use case, so I removed the other case
    perc_var_fred.append(mvavg_perc_var(curr_fred_norm, WINDOW))
    perc_var_fgreen.append(mvavg_perc_var(curr_fgreen_norm, WINDOW))
    # Compute percent variance due to the cell cycle (mean intensities)
    perc_var_cell.append(mvavg_perc_var(curr_ab_cell_norm,WINDOW))
    perc_var_nuc.append(mvavg_perc_var(curr_ab_nuc_norm, WINDOW))
    perc_var_cyto.append(mvavg_perc_var(curr_ab_cyto_norm, WINDOW))
    perc_var_mt.append(mvavg_perc_var(curr_mt_cell, WINDOW))
    # Compute percent variance due to the cell cycle (integrated intensities)
    perc_var_cell_int.append(mvavg_perc_var(curr_ab_cell_norm_int,WINDOW))
    perc_var_nuc_int.append(mvavg_perc_var(curr_ab_nuc_norm_int, WINDOW))
    perc_var_cyto_int.append(mvavg_perc_var(curr_ab_cyto_norm_int, WINDOW))
    perc_var_mt_int.append(mvavg_perc_var(curr_mt_cell_int, WINDOW))
    # Get x values for the moving average
    mvavg_xvals.append(mvavg_perc_var(curr_pol, WINDOW))

#%% Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
# Idea: create cutoffs for percent variance and 
# Execution: create cutoffs for perc_var and total variance per compartment and for integrated intensity and mean intensity
# Output: Graphs that illustrate the cutoffs (integrated, mean)
# Output: Overlap of the total variance cutoffs with the original filtering done manually
total_var_cutoff = np.mean(var_mt) + 2 * np.std(var_mt)
total_var_cutoff_int = np.mean(var_mt_int) + 2 * np.std(var_mt_int)
percent_var_cutoff = np.mean([a[0] for a in perc_var_mt if not np.isinf(a[0])]) + 2 * np.std([a[0] for a in perc_var_mt if not np.isinf(a[0])])
percent_var_cutoff_int = np.mean([a[0] for a in perc_var_mt_int if not np.isinf(a[0])]) + 2 * np.std([a[0] for a in perc_var_mt_int if not np.isinf(a[0])])

does_vary_cell = (np.array(var_cell) > total_var_cutoff) | (np.array(var_cell_int) > total_var_cutoff_int)
does_vary_nuc = (np.array(var_nuc) > total_var_cutoff) | (np.array(var_nuc_int) > total_var_cutoff_int)
does_vary_cyto = (np.array(var_cyto) > total_var_cutoff) | (np.array(var_cyto_int) > total_var_cutoff_int)

does_perc_vary_cell = (np.array([a[0] for a in perc_var_cell]) > percent_var_cutoff) | (np.array([a[0] for a in perc_var_cell_int]) > percent_var_cutoff_int)
does_perc_vary_nuc = (np.array([a[0] for a in perc_var_nuc]) > percent_var_cutoff) | (np.array([a[0] for a in perc_var_nuc_int]) > percent_var_cutoff_int)
does_perc_vary_cyto = (np.array([a[0] for a in perc_var_cyto]) > percent_var_cutoff) | (np.array([a[0] for a in perc_var_cyto_int]) > percent_var_cutoff_int)

df = pd.DataFrame({"well_plate" : u_well_plates, 
    "does_vary_cell":does_vary_cell,
    "does_vary_nuc":does_vary_nuc,
    "does_vary_cyto":does_vary_cyto,
    "does_perc_vary_cell":does_perc_vary_cell,
    "does_perc_vary_nuc":does_perc_vary_nuc,
    "does_perc_vary_cyto":does_perc_vary_cyto})
df.to_csv("output/wellplatevary.csv")

#%%
# Idea: generate a histogram of cell intensities
# Execution: pyplot
# Output: show histogram
plt.hist([np.std(a[1]) for a in perc_var_mt])
plt.show()

#%%

    # y for the specific compartment we are in, and the maximum y values
    yvals_compartment = get_val_compartment([yval_cell,yval_nuc,yval_cyto], curr_compartment)
    max_val_cell, max_val_nuc, max_val_cyto = np.max(yval_cell), np.max(yval_nuc), np.max(yval_cyto)


    #Set values for the current compartment
    #where the protein of interest is expressed
    curr_ab_norm = get_val_compartment([curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm], curr_compartment)
    perc_var_compartment = get_val_compartment([perc_var_cell,perc_var_nuc,perc_var_cyto],curr_compartment)
    var_compartment = get_val_compartment([var_cell,var_nuc,var_cyto],curr_compartment)
    gini_compartment = get_val_compartment([gini_cell,gini_nuc,gini_cyto],curr_compartment)
    max_val_compartment = get_val_compartment([max_val_cell, max_val_nuc, max_val_cyto], curr_compartment)

    #Compute the mean of the data surrounding the fit. 
    outname = os.path.join("output_devin", plate_well_to_ab[plate][well][2])
    if any(outsuff == i for i in OUTLIER_NAMES):
        fucci_plotting.hist2d_time(green_fucci,red_fucci,curr_fgreen,curr_fred,curr_pol,outname,outsuff,NBINS)

    if perc_var_fred>0.8 and perc_var_fgreen>0.8:
        fucci_plotting.hist2d_time(green_fucci,red_fucci,curr_fgreen,curr_fred,curr_pol,outname,outsuff,NBINS)

    #compute the model free max
    curr_ab_df = pd.DataFrame({outsuff:curr_ab_norm})
    curr_pol_df = pd.DataFrame({outsuff:curr_pol})
    winsize = np.min([len(curr_ab_norm),30])
    model_free_ab = curr_ab_df.rolling(window=winsize)
    model_free_ab = np.asarray(model_free_ab.mean())
    model_free_pol = curr_pol_df.rolling(window=winsize)
    model_free_pol = model_free_pol.mean()
    model_free_pol_np = np.asarray(model_free_pol)
    max_loc = np.argmax(model_free_ab[~np.isnan(model_free_ab)])
    max_pol = model_free_pol_np[max_loc]
    binned_values = []
    for xval in xvals: # x values for the plot
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
            [perc_var_compartment],
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
            temporal_mov_avg(curr_pol,curr_ab_norm,curr_xvals,yvals_compartment,outname,outsuff)
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
HIGHLIGHTS = ['ORC6','DUSP19','BUB1B','DPH2', 'FLI1']
# HIGHLIGHTS = ['ORC6','CCNE1','PHLDB1','DPH2']
# HIGHLIGHTS_MINOR = [ 'MCM10', 'ZNF32', 'JUN',  
#                 'DUSP19', 'CCNB1', 'AURKB', 'BUB1B', 'PAPSS1',
#                 'N6AMT1', 'FLI1']
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
sorted_expression_table_ccd = []
sorted_expression_table_all = []
for i,response in enumerate(model_free_list):
    sorted_gene_array.append(response[4]) # these are the expression values
    sorted_maxpol_array.append(response[3])
    sorted_ensg_array.append(response[6])
    for i, exp in enumerate(response[4]):
        sorted_expression_table_ccd.append([response[6], (i + 1) / len(response[4]), exp])
    curr_well = response[0]
    #set all other rows that share a loc to 1
    curr_locs = loc_mat[i,:]
    match_rows = np.greater(np.sum(loc_mat*curr_locs,axis=1),0)*response[3]
    sum_rows.append(np.sum(match_rows))
    prob_interact[i,:] = match_rows
    prot_list.append(response[1])

for i, response in enumerate(model_free_list_all):
    for i, exp in enumerate(response[4]):
        sorted_expression_table_all.append([response[6], (i + 1) / len(response[4]), exp])   

fig, ax = plt.subplots(figsize=(10, 10))
sc = ax.imshow(sorted_gene_array, interpolation='nearest')
np.savetxt("output/temporal_protein_expression_ccd.tsv", np.array(sorted_expression_table_ccd), fmt="%s", delimiter="\t")
np.savetxt("output/temporal_protein_expression_all.tsv", np.array(sorted_expression_table_all), fmt="%s", delimiter="\t")

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
# ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
ax.set_yticks(ytick_locs,minor=False)
ax.set_yticklabels(HIGHLIGHTS,minor=False)
# ax.set_yticks(ytick_locs_minor,minor=True)
# ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)
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
# (take only intersection of CCD proteins and CCD transcripts)
#
# prot_list # this one is the gene name (protein name) for the protein list
# sorted_ensg_array # this one is the sorted ensg for the proteins
# sorted_maxpol_array # this one is the protein max pol

name_ccdprot_ccdtrans_loc = [name_rna_list.index(name) for name in prot_list if name in name_rna_list and name in allccd_transcript_regulated_names]
name_prot_list = np.take(name_rna_list, name_ccdprot_ccdtrans_loc)
max_rna_avg_prot_pol = np.take(max_moving_avg_pol, name_ccdprot_ccdtrans_loc)
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
plt.xlabel("Division Cycle, hrs")
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

print(f"Median delay of RNA and protein expression time for CCD proteins: {TOT_LEN * np.median(diff_max_pol)}")
print(f"Median RNA expression time for CCD proteins: {TOT_LEN * np.median(max_rna_avg_prot_pol)}")
print(f"Median protein expression time for CCD proteins: {TOT_LEN * np.median(prot_maxpol_filter_array)}")
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

