#%% Imports
from utils import *
import numpy as np
from stretch_time import stretch_time
import scipy.optimize
import scipy.stats
import seaborn as sbn
import ProteinFucciPseudotime
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

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

ProteinFucciPseudotime.fucci_polar_coordinate_calculations(fucci_data, 
                           ab_nuc,ab_cyto,ab_cell,mt_cell,area_cell, area_nuc,well_plate,well_plate_imgnb, log_red_fucci_zeroc_rescale,log_green_fucci_zeroc_rescale)


#%% Measures of variance
# Idea: Calculate measures of variance, and show them in plots
# Execution: Now that we already have the data filtered for variability, this is just descriptive.
# Output: scatters of antibody vs microtubule variances by different measures of variaibility

# benjimini-hochberg multiple testing correction
# source: https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html

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
    
    cv_cell.append(scipy.stats.variation(curr_ab_cell))
    cv_nuc.append(scipy.stats.variation(curr_ab_nuc))
    cv_cyto.append(scipy.stats.variation(curr_ab_cyto))
    cv_mt.append(scipy.stats.variation(curr_mt_cell))
    
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
        
        curr_cv_cell_img.append(scipy.stats.variation(curr_ab_cell))
        curr_cv_nuc_img.append(scipy.stats.variation(curr_ab_nuc))
        curr_cv_cyto_img.append(scipy.stats.variation(curr_ab_cyto))
        curr_cv_mt_img.append(scipy.stats.variation(curr_mt_cell))
    
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

mmmm = np.concatenate((gini_comp, gini_mt))
cccc = (["Protein"] * len(gini_comp))
cccc.extend(["Microtubules"] * len(gini_mt))
boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=False, color="grey")
boxplot.set_xlabel("", size=36,fontname='Arial')
boxplot.set_ylabel("Gini", size=18,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=14)
plt.title("")
plt.savefig(f"figures/ProteinMicrotubuleGinis.png")
plt.savefig(f"figures/ProteinMicrotubuleGinis.pdf")
plt.close()
print(f"{scipy.stats.kruskal(gini_comp, gini_mt)[1]}: p-value for difference between protein and microtubule stainings")

mmmm = np.concatenate((var_comp, var_mt))
cccc = (["Protein"] * len(var_comp))
cccc.extend(["Microtubules"] * len(var_mt))
boxplot = sbn.boxplot(x=cccc, y=mmmm, showfliers=False, color="grey")
boxplot.set_xlabel("", size=36,fontname='Arial')
boxplot.set_ylabel("Variance", size=18,fontname='Arial')
boxplot.tick_params(axis="both", which="major", labelsize=14)
plt.title("")
plt.savefig(f"figures/ProteinMicrotubuleVariances.png")
plt.savefig(f"figures/ProteinMicrotubuleVariances.pdf")
plt.close()
print(f"{scipy.stats.kruskal(var_comp, var_mt)[1]}: p-value for difference between protein and microtubule stainings")

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


np_save_overwriting("output/pickles/mean_mean_comp.npy", mean_mean_comp)
np_save_overwriting("output/pickles/cv_comp.npy", cv_comp)
np_save_overwriting("output/pickles/gini_comp.npy", gini_comp)
np_save_overwriting("output/pickles/var_comp.npy", var_comp)
np_save_overwriting("output/pickles/cv_cell.npy", cv_cell)
np_save_overwriting("output/pickles/gini_cell.npy", gini_cell)
np_save_overwriting("output/pickles/var_cell.npy", var_cell)