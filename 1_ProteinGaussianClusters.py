#%% Imports
from imports import *
import numpy as np
from stretch_time import stretch_time
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar
from sklearn.neighbors import RadiusNeighborsRegressor
from sklearn.mixture import GaussianMixture
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
#       use the mean intensity and integrated intensity for different cutoffs (ab and microtubules)
# Exec: pandas
# Output: nothing

EMPTYWELLS = set(["B11_6745","C11_6745","D11_6745","E11_6745","F11_6745","G11_6745","H11_6745",
    "A12_6745","B12_6745","C12_6745","D12_6745","E12_6745","F12_6745","G12_6745"])

def read_raw_data():
    print("reading raw protein IF data")
    my_df1 = pd.read_csv("input/raw/nuc_predicted_prob_phases_mt_all_firstbatch_plates.csv")
    my_df2 = pd.read_csv("input/raw/nuc_predicted_prob_phases_190909.csv")
    my_df = pd.concat((my_df1, my_df2), sort=True)
    print("loaded raw data")
    return my_df

def read_sample_info(df):
    '''Get the metadata for all the samples'''
    plate = np.asarray(df.plate)
    u_plate = np.unique(plate)
    well_plate = np.asarray(df.well_plate)
    imgnb = np.asarray(df.ImageNumber)
    well_plate_imgnb = np.asarray([f"{wp}_{imgnb[i]}" for i,wp in enumerate(well_plate)])
    u_well_plates = np.unique(well_plate)
    ab_objnum = np.asarray(df.ObjectNumber)
    area_cell = np.asarray(df.Area_cell)
    area_nuc = np.asarray(df.AreaShape_Area)
    area_cyto = np.asarray(df.Area_cyto)
    name_df = pd.read_csv("input/processed/excel/Fucci_staining_summary_first_plates.csv")
    wppp1, ensggg1, abbb1, rrrr, cccc1 = list(name_df["well_plate"]), list(name_df["ENSG"]), list(name_df["Antibody"]), list(name_df["Results_final_update"]), list(name_df["Compartment"])
    name_df2 = pd.read_csv("input/processed/excel/Fucci_staining_review_variation_check.csv")
    wppp2, ensggg2, abbb2, cccc2 = list(name_df2["well_plate"]), list(name_df2["ENSG"]), list(name_df2["Antibody"]), list(name_df2["Compartment"])
    wppp, ensggg, abbb, cccc = wppp1 + wppp2, ensggg1 + ensggg2, abbb1 +  abbb2, cccc1 + cccc2
    ensg_dict = dict([(wppp[i], ensggg[i]) for i in range(len(wppp))])
    ab_dict = dict([(wppp[i], abbb[i]) for i in range(len(wppp))])
    result_dict = dict([(wppp[i], rrrr[i]) for i in range(len(wppp1))])
    compartment_dict = dict([(wppp[i], cccc[i]) for i in range(len(wppp))])
    ENSG = np.asarray([ensg_dict[wp] if wp in ensg_dict else "" for wp in well_plate])
    antibody = np.asarray([ab_dict[wp] if wp in ab_dict else "" for wp in well_plate])
    result = np.asarray([result_dict[wp] if wp in result_dict else "" for wp in well_plate])
    compartment = np.asarray([compartment_dict[wp] if wp in compartment_dict else "" for wp in well_plate])
    return plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment

def previous_results(u_well_plates, result_dict, ensg_dict, ab_dict):
    '''Process the results metadata into lists of previously annotated CCD proteins'''
    wp_ensg = np.asarray([ensg_dict[wp] if wp in ensg_dict else "" for wp in u_well_plates])
    wp_ab = np.asarray([ab_dict[wp] if wp in ab_dict else "" for wp in u_well_plates])
    wp_prev_ccd = np.asarray([wp in result_dict and result_dict[wp].startswith("ccd") for wp in u_well_plates])
    wp_prev_notccd = np.asarray([wp in result_dict and result_dict[wp].startswith("notccd") for wp in u_well_plates])
    wp_prev_negative = np.asarray([wp in result_dict and result_dict[wp].startswith("negative") for wp in u_well_plates])
    prev_ccd_ensg = wp_ensg[wp_prev_ccd]
    prev_notccd_ensg = wp_ensg[wp_prev_notccd]
    prev_negative_ensg = wp_ensg[wp_prev_negative]
    return wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg

# read raw data
my_df = read_raw_data()
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = read_sample_info(my_df)
wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% Idea: Filter the raw data
# Execution: Use manual annotations and nucleus size to filter samples and images
# Output: Filtered dataframe

def apply_manual_filtering(my_df, result_dict):
    '''Filter raw data based on manual annotations'''
    # filter some wells in the last plate didn't have anything.
    print(f"{len(my_df)}: number of cells before filtering empty wells")
    my_df = my_df[~my_df.well_plate.isin(EMPTYWELLS)]
    print(f"{len(my_df)}: number of cells after filtering empty wells")
    
    my_df_filtered = my_df
    print("filtering out of focus")
    oof = pd.read_csv("input/processed/excel/outoffocusimages.txt", header=None)[0]
    well_plate = np.asarray(my_df_filtered.well_plate)
    imgnb = np.asarray(my_df_filtered.ImageNumber)
    well_plate_imgnb = np.asarray([f"{wp}_{imgnb[i]}" for i,wp in enumerate(well_plate)])
    print(f"{len(my_df_filtered)}: number of cells before filtering out of focus images")
    my_df_filtered = my_df_filtered[~np.isin(well_plate_imgnb, oof)]
    print(f"{len(my_df_filtered)}: number of cells after filtering out of focus images")
    print("finished filtering")
    
    print("filtering negative staining")
    new_data_or_nonnegative_stain = [wp not in result_dict or (not result_dict[wp].lower().startswith("negative") and not wp.startswith("H12")) for wp in my_df_filtered.well_plate]
    print(f"{len(my_df_filtered)}: number of cells before filtering negative staining from first batch")
    my_df_filtered = my_df_filtered[new_data_or_nonnegative_stain]
    print(f"{len(my_df_filtered)}: number of cells after filtering negative staining from first batch")
    print("finished filtering")
    
    print("filtering bad fields of view (negative staining, unspecific, etc)")
    filterthese = pd.read_csv("input/processed/excel/FOV_ImgNum_Lookup.csv")
    badfov = filterthese["well_plate_imgnb"][(filterthese["UseImage"] == 0)]
    well_plate = np.asarray(my_df_filtered.well_plate)
    imgnb = np.asarray(my_df_filtered.ImageNumber)
    well_plate_imgnb = np.asarray([f"{wp}_{imgnb[i]}" for i,wp in enumerate(well_plate)])
    negative_controls = np.asarray([wp.startswith("H12") for wp in well_plate])
    print(f"{len(my_df_filtered)}: number of cells before filtering out of focus images")
    my_df_filtered = my_df_filtered[~np.isin(well_plate_imgnb, badfov) & ~negative_controls]
    print(f"{len(my_df_filtered)}: number of cells after filtering out of focus images")
    print("finished filtering")
    return my_df_filtered

def plot_areas(areas, title):
    '''histogram for areas of cell/nuc/cytoplasm'''
    bins = plt.hist(areas, bins=100, alpha=0.5)
    plt.vlines(np.mean(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) - 2 * np.std(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) + 2 * np.std(areas), 0, np.max(bins[0]))
    plt.title(title)
    plt.savefig(f"figures/areas{title}.png")
    plt.show()
    plt.close()
    
def apply_big_nucleus_filter(my_df):
    '''filter the super big nuclei'''
    area_cell, area_nuc, area_cyto = my_df.Area_cell, my_df.AreaShape_Area, my_df.Area_cyto
    plot_areas(area_cell, "area_cell")
    plot_areas(area_nuc, "area_nuc")
    plot_areas(area_cyto, "area_cyto")
    
    upper_nucleus_cutoff = np.mean(area_nuc) + 2 * np.std(area_nuc)

    my_df_filtered = my_df
    print("filtering super big nuclei")
    cell_passes_nucleus_filter = my_df_filtered.AreaShape_Area < upper_nucleus_cutoff
    print(f"{len(my_df_filtered)}: number of cells before filtering out super big nuclei")
    my_df_filtered = my_df_filtered[cell_passes_nucleus_filter]
    print(f"{len(my_df_filtered)}: number of cells after filtering out super big nuclei")
    print("finished filtering on nuclei")
    
    area_cell_filtered, area_nuc_filtered, area_cyto_filtered = my_df_filtered.Area_cell, my_df_filtered.AreaShape_Area, my_df_filtered.Area_cyto
    plot_areas(area_cell_filtered, "area_cell_filtered")
    plot_areas(area_nuc_filtered, "area_nuc_filtered")
    plot_areas(area_cyto_filtered, "area_cyto_filtered")
    return my_df_filtered

MIN_CELL_COUNT = 60
def apply_cell_count_filter(my_df):
    '''filter low cell counts per sample'''
    my_df_filtered = my_df
    well_plate = np.asarray(my_df_filtered.well_plate)
    cell_counts = [sum(my_df.well_plate == wp) for wp in well_plate]
    print("filtering low cell counts")
    my_df_filtered = my_df_filtered[cell_counts >= MIN_CELL_COUNT]
    print(f"{len(my_df)}: number of cells before filtering out samples with < {MIN_CELL_COUNT} cells")
    print(f"{len(my_df_filtered)}: number of cells after filtering out samples with < {MIN_CELL_COUNT} cells")
    print("finished filtering on cell count")
    return my_df_filtered

my_df_filtered = apply_manual_filtering(my_df, result_dict)
my_df_filtered = apply_big_nucleus_filter(my_df_filtered)
my_df_filtered = apply_cell_count_filter(my_df_filtered)
my_df_filtered.to_csv("input/processed/python/nuc_predicted_prob_phases_filtered.csv")

plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = read_sample_info(my_df_filtered)
wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% 
# Idea: Filter for variation and get compartments
# Execution: Use annotated variation and compartment information
# Output: none
def apply_variation_filter(my_df_filtered, result_dict):
    '''Separate the varying and nonvarying samples'''
    my_df_filtered_variation, my_df_filtered_novariation = my_df_filtered, my_df_filtered
    variable_firstbatch = np.asarray([wp in result_dict and not result_dict[wp].replace(" ","").startswith("novariation") for wp in my_df_filtered.well_plate])
    
    varann_secondbatch = pd.read_csv("input/processed/excel/SecondBatchVariableLookup.csv")
    variable_ann_secondbatch = np.asarray([str(vv).lower().startswith("yes") for vv in varann_secondbatch["IsVariable"]])
    variable_wp_secondbatch = np.asarray(varann_secondbatch["well_plate"][variable_ann_secondbatch])
    variable_secondbatch = np.isin(my_df_filtered.well_plate, variable_wp_secondbatch)
    
    my_df_filtered_variation = my_df_filtered[variable_firstbatch | variable_secondbatch]
    my_df_filtered_novariation = my_df_filtered[~(variable_firstbatch | variable_secondbatch)]
    print(f"{len(my_df)}: number of cells before filtering for variation")
    print(f"{len(my_df_filtered_variation)}: number of cells in samples with variation")
    print(f"{len(my_df_filtered_novariation)}: number of cells in samples without variation")
    return my_df_filtered_variation, my_df_filtered_novariation

def metacompartments(u_well_plates, compartment_dict, my_df_filtered_variation):
    '''Get the compartments for the unique wellplates'''
    wp_iscell = np.asarray([compartment_dict[wp].lower().startswith("cell") if wp in compartment_dict else False for wp in u_well_plates])
    wp_isnuc = np.asarray([compartment_dict[wp].lower().startswith("nuc") if wp in compartment_dict else False for wp in u_well_plates])
    wp_iscyto = np.asarray([compartment_dict[wp].lower().startswith("cyto") if wp in compartment_dict else False for wp in u_well_plates])
    
    wp_nocompartmentinfo = ~wp_iscell & ~wp_isnuc & ~wp_iscyto
    print(f"{sum(wp_nocompartmentinfo)}: samples without compartment information; to be filtered since they're biologically defined as CCD and not included in the analysis")
    print(f"{len(my_df_filtered_variation)}: number of cells before filtering for compartment information")
    my_df_filtered_compartmentvariation = my_df_filtered_variation[~np.isin(my_df_filtered_variation.well_plate, u_well_plates[wp_nocompartmentinfo])]
    print(f"{len(my_df_filtered_compartmentvariation)}: number of cells before filtering for compartment information")
    return wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation

my_df_filtered_variation, my_df_filtered_novariation = apply_variation_filter(my_df_filtered, result_dict)
#my_df_filtered_variation.to_csv("input/processed/python/nuc_predicted_prob_phases_filtered_variation.csv")
#my_df_filtered_novariation.to_csv("input/processed/python/nuc_predicted_prob_phases_filtered_novariation.csv")

plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = read_sample_info(my_df_filtered_variation)

wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation = metacompartments(u_well_plates, compartment_dict, my_df_filtered_variation)
plate, u_plate, well_plate, well_plate_imgnb, u_well_plates, ab_objnum, area_cell, area_nuc, area_cyto, ensg_dict, ab_dict, result_dict, compartment_dict, ENSG, antibody, result, compartment = read_sample_info(my_df_filtered_compartmentvariation)
wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation = metacompartments(u_well_plates, compartment_dict, my_df_filtered_compartmentvariation)

wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg = previous_results(u_well_plates, result_dict, ensg_dict, ab_dict)

#%% 
# Idea: Get and process intensities
# Execution: get intensities; zero center fucci intensities
# Output: Fucci plot

# 0: use mean, 
# 1: use integrated,
# 2: use integrated / nucleus area
INTENSITY_SWITCH = 0 # small cells are often brighter and this is just because they have rounded up and are thicker, so use integrated
def read_sample_data(df):
    # Antibody data (mean intensity)
    ab_nuc = np.asarray([df.Intensity_MeanIntensity_ResizedAb, 
                         df.Intensity_IntegratedIntensity_ResizedAb, 
                         df.Intensity_IntegratedIntensity_ResizedAb / df.AreaShape_Area][INTENSITY_SWITCH])
    ab_cyto = np.asarray([df.Mean_ab_Cyto, 
                          df.Integrated_ab_cyto, 
                          df.Integrated_ab_cyto / df.AreaShape_Area][INTENSITY_SWITCH])
    ab_cell = np.asarray([df.Mean_ab_cell, 
                          df.Integrated_ab_cell, 
                          df.Integrated_ab_cell / df.AreaShape_Area][INTENSITY_SWITCH])
    mt_cell = np.asarray([df.Mean_mt_cell, 
                          df.Integrated_mt_cell, 
                          df.Integrated_mt_cell / df.AreaShape_Area][INTENSITY_SWITCH])

    # Fucci data (mean intensity)
    green_fucci = np.asarray(df.Intensity_MeanIntensity_CorrResizedGreenFUCCI)
    red_fucci = np.asarray(df.Intensity_MeanIntensity_CorrResizedRedFUCCI)
    return ab_nuc, ab_cyto, ab_cell, mt_cell, green_fucci, red_fucci

ab_nuc, ab_cyto, ab_cell, mt_cell, green_fucci, red_fucci = read_sample_data(my_df_filtered_compartmentvariation)

# Zero center and rescale FUCCI data in the log space
log_green_fucci, log_red_fucci = np.log10(green_fucci), np.log10(red_fucci)
wp_p_dict = dict([(str(p), plate == p) for p in u_plate])
logmed_green_fucci_p = dict([(str(p), np.log10(np.median(green_fucci[wp_p_dict[str(p)]]))) for p in u_plate])
logmed_red_fucci_p = dict([(str(p), np.log10(np.median(red_fucci[wp_p_dict[str(p)]]))) for p in u_plate])
logmed_green_fucci = np.array([logmed_green_fucci_p[wp.split("_")[1]] for wp in well_plate])
logmed_red_fucci = np.array([logmed_red_fucci_p[wp.split("_")[1]] for wp in well_plate])
log_green_fucci_zeroc = np.array(log_green_fucci) - logmed_green_fucci
log_red_fucci_zeroc = np.array(log_red_fucci) - logmed_red_fucci
log_green_fucci_zeroc_rescale = (log_green_fucci_zeroc - np.min(log_green_fucci_zeroc)) / np.max(log_green_fucci_zeroc)
log_red_fucci_zeroc_rescale = (log_red_fucci_zeroc - np.min(log_red_fucci_zeroc)) / np.max(log_red_fucci_zeroc)
fucci_data = np.column_stack([log_green_fucci_zeroc_rescale,log_red_fucci_zeroc_rescale])
plt.hist2d(log_green_fucci_zeroc_rescale,log_red_fucci_zeroc_rescale,bins=200)
plt.savefig("figures/FucciPlotProteinIFData_unfiltered.png")
plt.show()
plt.close()

#%%
# Idea: Gaussian clustering per plate to identify G1/S/G2 and do kruskal test for variance
# Exec: sklearn.mixture.GaussianMixture & scipy.stats.kruskal
# Output: FDR for cell cycle variation per well per compartment
gaussian = GaussianMixture(n_components=3, random_state=1, max_iter=500)
cluster_labels = gaussian.fit_predict(np.array([log_green_fucci_zeroc_rescale, log_red_fucci_zeroc_rescale]).T)
for cluster in range(3):
    plt.hist2d(log_green_fucci_zeroc_rescale[cluster_labels == cluster],log_red_fucci_zeroc_rescale[cluster_labels == cluster],bins=200)
    plt.title(f"Gaussian clustered data, cluster {cluster}")
    plt.savefig(f"figures/FucciPlotProteinIFData_unfiltered_Gauss{cluster}.png")
    plt.show()
    plt.close()

# G1 is cluster 0; S-ph is cluster 1; G2 is cluster 2
wp_cell_kruskal, wp_nuc_kruskal, wp_cyto_kruskal = [],[],[]
g1 = cluster_labels == 0
sph = cluster_labels == 1
g2 = cluster_labels == 2
for wp in u_well_plates:
    curr_wp_g1 = (well_plate == wp) & g1
    curr_wp_sph = (well_plate == wp) & sph
    curr_wp_g2 = (well_plate == wp) & g2
    wp_cell_kruskal.append(scipy.stats.kruskal(ab_cell[curr_wp_g1], ab_cell[curr_wp_sph], ab_cell[curr_wp_g2])[1])
    wp_nuc_kruskal.append(scipy.stats.kruskal(ab_nuc[curr_wp_g1], ab_nuc[curr_wp_sph], ab_nuc[curr_wp_g2])[1])
    wp_cyto_kruskal.append(scipy.stats.kruskal(ab_cyto[curr_wp_g1], ab_cyto[curr_wp_sph], ab_cyto[curr_wp_g2])[1])

# benjimini-hochberg multiple testing correction
# source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
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

alpha_gauss = 0.05
wp_cell_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_cell = benji_hoch(alpha_gauss, wp_cell_kruskal)
wp_nuc_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_nuc = benji_hoch(alpha_gauss, wp_nuc_kruskal)
wp_cyto_kruskal_gaussccd_adj, wp_pass_gaussccd_bh_cyto = benji_hoch(alpha_gauss, wp_cyto_kruskal)

def values_comp(values_cell, values_nuc, values_cyto, wp_iscell, wp_isnuc, wp_iscyto):    
    values_comp = np.empty_like(values_cell)
    values_comp[wp_iscell] = np.array(values_cell)[wp_iscell]
    values_comp[wp_isnuc] = np.array(values_nuc)[wp_isnuc]
    values_comp[wp_iscyto] = np.array(values_cyto)[wp_iscyto]
    return values_comp

wp_comp_kruskal_gaussccd_p = values_comp(wp_cell_kruskal, wp_nuc_kruskal, wp_cyto_kruskal, wp_iscell, wp_isnuc, wp_iscyto)
wp_comp_kruskal_gaussccd_adj, wp_pass_kruskal_gaussccd_bh_comp = benji_hoch(alpha_gauss, wp_comp_kruskal_gaussccd_p)

# BenjiHoch is actually pretty liberal for this dataset. What about bonferroni?
# But it gets rid of things that seem CCD, so use BH
#wp_comp_kruskal_gaussccd_bonfadj, wp_bonfpass_kruskal_gaussccd_comp = bonf(alpha_gauss, wp_comp_kruskal_gaussccd_p)

print(f"{len(wp_pass_gaussccd_bh_cell)}: number of genes tested")
#print(f"{sum(wp_pass_gaussccd_bh_cell)}: number of passing genes at 5% FDR in cell")
#print(f"{sum(wp_pass_gaussccd_bh_cyto)}: number of passing genes at 5% FDR in cytoplasm")
#print(f"{sum(wp_pass_gaussccd_bh_nuc)}: number of passing genes at 5% FDR in nucleus")
print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp)}: number of passing genes at {alpha_gauss*100}% FDR in compartment")

# address gene redundancy
wp_ensg_counts = np.array([sum([1 for eeee in wp_ensg if eeee == ensg]) for ensg in wp_ensg])
ensg_is_duplicated = wp_ensg_counts > 1
duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
duplicated_ensg_pairs = [u_well_plates[wp_ensg == ensg] for ensg in duplicated_ensg]
print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp[~ensg_is_duplicated])}: number of passing genes at {alpha_gauss*100}% FDR in compartment (no replicate)")
duplicated_ensg_ccd = np.array([sum(wp_pass_kruskal_gaussccd_bh_comp[wp_ensg == ensg]) for ensg in duplicated_ensg])
print(f"{sum(duplicated_ensg_ccd == 2)}: number of CCD genes shown to be CCD in both replicates")
print(f"{sum(duplicated_ensg_ccd == 1)}: number of CCD genes shown to be CCD in just one replicate")
print(f"{sum(duplicated_ensg_ccd == 0)}: number of CCD genes shown to be non-CCD in both replicate")

# any antibody redundancy?
wp_ab_counts = np.array([sum([1 for eeee in wp_ensg if eeee == ensg]) for ensg in wp_ab])
ab_is_duplicated = wp_ab_counts > 1
duplicated_ab = np.unique(wp_ab[ab_is_duplicated])
print(f"{sum(wp_pass_kruskal_gaussccd_bh_comp[~ab_is_duplicated])}: number of passing genes at {alpha_gauss*100}% FDR in compartment (no replicate)")
duplicated_ab_ccd = np.array([sum(wp_pass_kruskal_gaussccd_bh_comp[wp_ab == ensg]) for ab in duplicated_ab])
print(f"{sum(duplicated_ab_ccd == 2)}: number of duplicated antibodies shown to be CCD in both replicates")
print(f"{sum(duplicated_ab_ccd == 1)}: number of duplicated antibodies shown to be CCD in just one replicate")
print(f"{sum(duplicated_ab_ccd == 0)}: number of duplicated antibodies shown to be non-CCD in both replicate")

#%% Pickle the results
def np_save_overwriting(fn, arr):
    with open(fn,"wb") as f:    
        np.save(f, arr, allow_pickle=True)

np_save_overwriting("output/pickles/u_plate.npy", u_plate)
np_save_overwriting("output/pickles/u_well_plates.npy", u_well_plates)
np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)
np_save_overwriting("output/pickles/wp_ab.npy", wp_ab)
np_save_overwriting("output/pickles/area_cell.npy", area_cell)
np_save_overwriting("output/pickles/area_nuc.npy", area_nuc)
np_save_overwriting("output/pickles/area_cyto.npy", area_cyto)
np_save_overwriting("output/pickles/wp_prev_ccd.npy", wp_prev_ccd)
np_save_overwriting("output/pickles/wp_prev_notccd.npy", wp_prev_notccd)
np_save_overwriting("output/pickles/wp_prev_negative.npy", wp_prev_negative)
np_save_overwriting("output/pickles/prev_ccd_ensg.npy", prev_ccd_ensg)
np_save_overwriting("output/pickles/prev_notccd_ensg.npy", prev_notccd_ensg)
np_save_overwriting("output/pickles/prev_negative_ensg.npy", prev_negative_ensg)
np_save_overwriting("output/pickles/well_plate.npy", well_plate)
np_save_overwriting("output/pickles/well_plate_imgnb.npy", well_plate_imgnb)
np_save_overwriting("output/pickles/ab_nuc.npy", ab_nuc)
np_save_overwriting("output/pickles/ab_cyto.npy", ab_cyto)
np_save_overwriting("output/pickles/ab_cell.npy", ab_cell)
np_save_overwriting("output/pickles/mt_cell.npy", mt_cell)
np_save_overwriting("output/pickles/green_fucci.npy", green_fucci)
np_save_overwriting("output/pickles/red_fucci.npy", red_fucci)
np_save_overwriting("output/pickles/log_green_fucci_zeroc.npy", log_green_fucci_zeroc)
np_save_overwriting("output/pickles/log_red_fucci_zeroc.npy", log_red_fucci_zeroc)
np_save_overwriting("output/pickles/log_green_fucci_zeroc_rescale.npy", log_green_fucci_zeroc_rescale)
np_save_overwriting("output/pickles/log_red_fucci_zeroc_rescale.npy", log_red_fucci_zeroc_rescale)
np_save_overwriting("output/pickles/wp_comp_kruskal_gaussccd_adj.npy", wp_comp_kruskal_gaussccd_adj)
np_save_overwriting("output/pickles/wp_pass_kruskal_gaussccd_bh_comp.npy", wp_pass_kruskal_gaussccd_bh_comp)
np_save_overwriting("output/pickles/fucci_data.npy", fucci_data)
np_save_overwriting("output/pickles/wp_iscell.npy", wp_iscell)
np_save_overwriting("output/pickles/wp_isnuc.npy", wp_isnuc)
np_save_overwriting("output/pickles/wp_iscyto.npy", wp_iscyto)

