#%% Imports
from utils import *
import numpy as np
import sklearn.mixture
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

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
pol_sort_fred = np.load("output/pickles/pol_sort_fred.npy", allow_pickle=True)
pol_sort_fgreen = np.load("output/pickles/pol_sort_fgreen.npy", allow_pickle=True)

wp_iscell = np.load("output/pickles/wp_iscell.npy", allow_pickle=True)
wp_isnuc = np.load("output/pickles/wp_isnuc.npy", allow_pickle=True)
wp_iscyto = np.load("output/pickles/wp_iscyto.npy", allow_pickle=True)
print("loaded")

        
#%%
# Idea: process the well data
# Exec: use Devin's code
# Output: the moving average plots for each gene studied




#%% Do the percent variance values match up with what we measured before for the genes?
# Idea: take the perc var values from Devin's analysis and compare them to the ones now
# Execution: run old code and save the data
# Output: plot of percvar vs percvar



    
#%% Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
# Idea: create cutoffs for percent variance and 
# Execution: create cutoffs for perc_var and total variance per compartment and for integrated intensity and mean intensity
# Output: Graphs that illustrate the cutoffs (integrated, mean)
# Output: Overlap of the total variance cutoffs with the original filtering done manually



# Output a list of the genes that show variation


# Accounting for biologically CCD ones



    
pd.DataFrame({"ENSG":np.unique(wp_ensg[wp_prev_ccd & ~ccd_comp])}).to_csv("output/DianaCCDMissingProteins.csv")
pd.DataFrame({"ENSG":np.unique(wp_ensg[~wp_prev_ccd & ccd_comp])}).to_csv("output/NewToDianaCCDProteins.csv") 

#%% Pickle the results needed later
np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp) # removed ones passing in only one replicate
np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp) # removed ones passing in only one replicate
np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)

pd.DataFrame({"gene": wp_ensg[ccd_comp]}).to_csv("output/picklestxt/ccd_compartment_ensg.txt", index=False)
pd.DataFrame({"gene": wp_ensg[nonccd_comp]}).to_csv("output/picklestxt/nonccd_compartment_ensg.txt", index=False)

#%% Figures of merit
