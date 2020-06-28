# -*- coding: utf-8 -*-
"""
Loading the raw immunofluorescence data for single-cell protein expression.
Filtering the data for:
    - Sample with low cell count
    - Out of focus images removed
    - Antibodies that failed validation are removed
    - Small nuclei ar removed (mitotic cells) prior to this analylsis
    - Large nuclei are removed (segmentation artificats)
    - Only proteins displaying variability are kept for single-cell variability analysis (using manual annotations)

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable


EMPTYWELLS = set(["B11_6745","C11_6745","D11_6745","E11_6745","F11_6745","G11_6745","H11_6745",
    "A12_6745","B12_6745","C12_6745","D12_6745","E12_6745","F12_6745","G12_6745"]) 
# EMPTYWELLS: These wells on the last plate didn't have cells; the segmentation algorithm still annotated some, so remove them
MIN_CELL_COUNT = 60 # Minimum number of cells per sample required for cell cycle analysis with pseudotime

# 0: use mean, (testing intensity that's already normalized for cell size)
# 1: use integrated, (testing integrated because it reflects that small cells are often brighter because they have rounded up and are thicker)
# 2: use integrated / nucleus area (testing a normalization by nucleus size, which may be more consistent in segmentation)
INTENSITY_SWITCH = 0 # cell size does increase into G2 and that has a substantial effect, so use the mean intensity, which is well behaved

def read_raw_data():
    '''Read in the raw protein IF data'''
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
    
    # Pickle the results
    utils.np_save_overwriting("output/pickles/u_plate.npy", u_plate)
    utils.np_save_overwriting("output/pickles/u_well_plates.npy", u_well_plates)
    utils.np_save_overwriting("output/pickles/area_cell.npy", area_cell)
    utils.np_save_overwriting("output/pickles/area_nuc.npy", area_nuc)
    utils.np_save_overwriting("output/pickles/area_cyto.npy", area_cyto)
    utils.np_save_overwriting("output/pickles/well_plate.npy", well_plate)
    utils.np_save_overwriting("output/pickles/well_plate_imgnb.npy", well_plate_imgnb)
    
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
    
    # Pickle the results
    utils.np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)
    utils.np_save_overwriting("output/pickles/wp_ab.npy", wp_ab)
    utils.np_save_overwriting("output/pickles/wp_prev_ccd.npy", wp_prev_ccd)
    utils.np_save_overwriting("output/pickles/wp_prev_notccd.npy", wp_prev_notccd)
    utils.np_save_overwriting("output/pickles/wp_prev_negative.npy", wp_prev_negative)
    utils.np_save_overwriting("output/pickles/prev_ccd_ensg.npy", prev_ccd_ensg)
    utils.np_save_overwriting("output/pickles/prev_notccd_ensg.npy", prev_notccd_ensg)
    utils.np_save_overwriting("output/pickles/prev_negative_ensg.npy", prev_negative_ensg)
    
    return wp_ensg, wp_ab, wp_prev_ccd, wp_prev_notccd, wp_prev_negative, prev_ccd_ensg, prev_notccd_ensg, prev_negative_ensg

def apply_manual_filtering(my_df, result_dict, ab_dict):
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
    
    print("filtering failed antibodies")
    failedab = np.genfromtxt("input/processed/manual/failedab.txt", dtype='str')
    print(f"{len(my_df_filtered)}: number of cells before filtering antibodies failed in HPAv19")
    my_df_filtered = my_df_filtered[~np.isin([ab_dict[wp] for wp in my_df_filtered.well_plate], failedab)]
    print(f"{len(my_df_filtered)}: number of cells after filtering antibodies failed in HPAv19")
    print("finished filtering")
    
    print("filtering mitotic proteins")
    mitoticab = np.genfromtxt("input/processed/manual/mitotic_n_microtubules_to_remove.txt", dtype='str')
    print(f"{len(my_df_filtered)}: number of cells before filtering mitotic/microtubule proteins")
    my_df_filtered = my_df_filtered[~np.isin([ab_dict[wp] for wp in my_df_filtered.well_plate], mitoticab)]
    print(f"{len(my_df_filtered)}: number of cells after filtering mitotic/microtubule proteins")
    print("finished filtering")
    
    return my_df_filtered
    
def plot_areas(areas, title):
    '''Histogram for areas of cell/nuc/cytoplasm'''
    bins = plt.hist(areas, bins=100, alpha=0.5)
    plt.vlines(np.mean(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) - 2 * np.std(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) + 2 * np.std(areas), 0, np.max(bins[0]))
    plt.title(title)
    plt.ylabel("Count")
    plt.xlabel("Area")
    plt.savefig(f"figures/areas{title}.png")
    plt.show()
    plt.close()

def apply_big_nucleus_filter(my_df):
    '''Filter the super big nuclei'''
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

def apply_cell_count_filter(my_df):
    '''Filter low cell counts per sample'''
    my_df_filtered = my_df
    well_plate = np.asarray(my_df_filtered.well_plate)
    u_well_plates = np.unique(my_df_filtered.well_plate)
    cell_count_dict = dict((wp, sum(my_df.well_plate == wp)) for wp in u_well_plates)
    cell_counts = np.array([cell_count_dict[wp] for wp in well_plate])
    print("filtering low cell counts")
    my_df_filtered = my_df_filtered[cell_counts >= MIN_CELL_COUNT]
    print(f"{len(my_df)}: number of cells before filtering out samples with < {MIN_CELL_COUNT} cells")
    print(f"{len(my_df_filtered)}: number of cells after filtering out samples with < {MIN_CELL_COUNT} cells")
    print("finished filtering on cell count")
    return my_df_filtered

def apply_variation_filter(my_df_filtered, result_dict, unfiltered_df):
    '''Separate the varying and nonvarying samples'''
    my_df_filtered_variation, my_df_filtered_novariation = my_df_filtered, my_df_filtered
    variable_firstbatch = np.asarray([wp in result_dict and not result_dict[wp].replace(" ","").startswith("novariation") for wp in my_df_filtered.well_plate])
    
    varann_secondbatch = pd.read_csv("input/processed/excel/SecondBatchVariableLookup.csv")
    variable_ann_secondbatch = np.asarray([str(vv).lower().startswith("yes") for vv in varann_secondbatch["IsVariable"]])
    variable_wp_secondbatch = np.asarray(varann_secondbatch["well_plate"][variable_ann_secondbatch])
    variable_secondbatch = np.isin(my_df_filtered.well_plate, variable_wp_secondbatch)
    
    my_df_filtered_variation = my_df_filtered[variable_firstbatch | variable_secondbatch]
    my_df_filtered_novariation = my_df_filtered[~(variable_firstbatch | variable_secondbatch)]
    print(f"{len(unfiltered_df)}: number of cells before filtering for variation")
    print(f"{len(my_df_filtered_variation)}: number of cells in samples with variation")
    print(f"{len(my_df_filtered_novariation)}: number of cells in samples without variation")
    return my_df_filtered_variation, my_df_filtered_novariation

def metacompartments(u_well_plates, compartment_dict, my_df_filtered_variation):
    '''Get the compartments for the unique wellplates'''
    wp_iscell = np.asarray([compartment_dict[wp].lower().startswith("cell") if wp in compartment_dict else False for wp in u_well_plates])
    wp_isnuc = np.asarray([compartment_dict[wp].lower().startswith("nuc") if wp in compartment_dict else False for wp in u_well_plates])
    wp_iscyto = np.asarray([compartment_dict[wp].lower().startswith("cyto") if wp in compartment_dict else False for wp in u_well_plates])
    
    # Pickle the results
    utils.np_save_overwriting("output/pickles/wp_iscell.npy", wp_iscell)
    utils.np_save_overwriting("output/pickles/wp_isnuc.npy", wp_isnuc)
    utils.np_save_overwriting("output/pickles/wp_iscyto.npy", wp_iscyto)

    wp_nocompartmentinfo = ~wp_iscell & ~wp_isnuc & ~wp_iscyto
    print(f"{sum(wp_nocompartmentinfo)}: samples without compartment information; to be filtered since they're biologically defined as CCD and not included in the analysis")
    print(f"{len(my_df_filtered_variation)}: number of cells before filtering for compartment information")
    my_df_filtered_compartmentvariation = my_df_filtered_variation[~np.isin(my_df_filtered_variation.well_plate, u_well_plates[wp_nocompartmentinfo])]
    print(f"{len(my_df_filtered_compartmentvariation)}: number of cells before filtering for compartment information")
    return wp_iscell, wp_isnuc, wp_iscyto, my_df_filtered_compartmentvariation

def read_sample_data(df):
    '''Read antibody intensity data for each sample and save it to a file for later use.'''
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
    
    # Pickle the results
    utils.np_save_overwriting("output/pickles/ab_nuc.npy", ab_nuc)
    utils.np_save_overwriting("output/pickles/ab_cyto.npy", ab_cyto)
    utils.np_save_overwriting("output/pickles/ab_cell.npy", ab_cell)
    utils.np_save_overwriting("output/pickles/mt_cell.npy", mt_cell)
    utils.np_save_overwriting("output/pickles/green_fucci.npy", green_fucci)
    utils.np_save_overwriting("output/pickles/red_fucci.npy", red_fucci)

    return ab_nuc, ab_cyto, ab_cell, mt_cell, green_fucci, red_fucci