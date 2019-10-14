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
log_green_fucci = np.log10(green_fucci_new)
red_fucci = np.asarray(my_df.Intensity_MeanIntensity_CorrResizedRedFUCCI)
log_red_fucci = np.log10(red_fucci)
fucci_data = np.column_stack([log_green_fucci,log_red_fucci])
# plt.hist2d(np.log10(green_fucci_new),np.log10(red_fucci),bins=200)
# plt.savefig("figures/FucciPlotProteinIFData_unfiltered.png")
# plt.show()
# plt.close()

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

#%% Negative staining (mean intensity max per image)
# Negative staining (mean intensity sum per image)
ab_cell_neg_max = np.asarray([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if wp in result_dict and result_dict[wp].startswith("neg")])
ab_cell_neg_sum = np.asarray([sum(ab_cell[well_plate == wp]) for wp in u_well_plates if wp in result_dict and result_dict[wp].startswith("neg")])
ab_cell_max = np.asarray([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if wp in result_dict and not result_dict[wp].startswith("neg") and not wp.startswith("H12")])
ab_cell_sum = np.asarray([sum(ab_cell[well_plate == wp]) for wp in u_well_plates if wp in result_dict and not result_dict[wp].startswith("neg") and not wp.startswith("H12")])
upper_neg_max_cutoff = np.mean(np.log10(ab_cell_neg_max)) + 0.8 * np.std(np.log10(ab_cell_neg_max))
lower_neg_max_cutoff = np.mean(np.log10(ab_cell_neg_max)) - 0.8 * np.std(np.log10(ab_cell_neg_max))
upper_neg_sum_cutoff = np.mean(np.log10(ab_cell_neg_sum)) + 0.8 * np.std(np.log10(ab_cell_neg_sum))
lower_neg_sum_cutoff = np.mean(np.log10(ab_cell_neg_sum)) - 0.8 * np.std(np.log10(ab_cell_neg_sum))
# bins = plt.hist(np.log10(ab_cell_neg_max), bins=100, alpha=0.5, label="Negative Staining")
# bins = plt.hist(np.log10(ab_cell_max), bins=100, alpha=0.5, label="Positive Staining")
# plt.vlines(upper_neg_max_cutoff, 0, np.max(bins[0]))
# plt.vlines(lower_neg_max_cutoff, 0, np.max(bins[0]))
# plt.title("Staining (max mean intensity per image)")
# plt.show()
# plt.close()
# bins = plt.hist(np.log10(ab_cell_neg_sum), bins=100, alpha=0.5, label="Negative Staining")
# bins = plt.hist(np.log10(ab_cell_sum), bins=100,alpha=0.5,  label="Positive Staining")
# plt.vlines(upper_neg_sum_cutoff, 0, np.max(bins[0]))
# plt.vlines(lower_neg_sum_cutoff, 0, np.max(bins[0]))
# plt.title("Staining (sum mean intensity per image)")
# plt.show()
# plt.close()
plt.scatter(np.log10(ab_cell_neg_max), np.log10(ab_cell_neg_sum), alpha=1)
plt.scatter(np.log10(ab_cell_max), np.log10(ab_cell_sum), alpha=0.1)
plt.vlines(upper_neg_max_cutoff, -0.5, 2)
plt.hlines(upper_neg_sum_cutoff, -2.5, 0)
plt.show()
plt.close()
print(f"{sum((np.log10(ab_cell_neg_max) < upper_neg_max_cutoff) & (np.log10(ab_cell_neg_sum) < upper_neg_sum_cutoff)) / len(ab_cell_neg_sum)}: percent of neg stains removed")
print(f"{sum((np.log10(ab_cell_max) < upper_neg_max_cutoff) & (np.log10(ab_cell_sum) < upper_neg_sum_cutoff)) / len(ab_cell_sum)}: percent of pos stains removed")

#%% Filter negative stains
ab_cell_max = np.asarray([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if not wp.startswith("H12")])
ab_cell_sum = np.asarray([sum(ab_cell[well_plate == wp]) for wp in u_well_plates if not wp.startswith("H12")])
ab_cell_max_negctrl = np.asarray([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if wp.startswith("H12")])
ab_cell_sum_negctrl = np.asarray([sum(ab_cell[well_plate == wp]) for wp in u_well_plates if wp.startswith("H12")])
plt.scatter(np.log10(ab_cell_max_negctrl), np.log10(ab_cell_sum_negctrl), alpha=1)
plt.scatter(np.log10(ab_cell_max), np.log10(ab_cell_sum), alpha=0.1)
plt.vlines(upper_neg_max_cutoff, -0.5, 2)
plt.hlines(upper_neg_sum_cutoff, -2.5, 0)
plt.show()
plt.close()
print(f"{sum((np.log10(ab_cell_max_negctrl) < upper_neg_max_cutoff) & (np.log10(ab_cell_sum_negctrl) < upper_neg_sum_cutoff)) / len(ab_cell_sum_negctrl)}: percent of neg stains removed")
print(f"{sum((np.log10(ab_cell_max) < upper_neg_max_cutoff) & (np.log10(ab_cell_sum) < upper_neg_sum_cutoff)) / len(ab_cell_sum)}: percent of pos stains removed")



#%% Do I need to offset the antibody intensities, too?
for p in u_plate:
    # plt.hist([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if not wp.startswith("H12") and wp.endswith(str(p))], bins=np.logspace(np.log10(0.0001), np.log10(1), 50))
    plt.hist(np.log10([np.max(ab_cell[well_plate == wp]) for wp in u_well_plates if not wp.startswith("H12") and wp.endswith(str(p))]), bins=np.linspace(-2.5, 0, 100))
    plt.title(f"{p} hist")
    plt.show()
    plt.savefig(f"figures/{p}MaxAbIntensities.png")
    plt.close()

#%%
