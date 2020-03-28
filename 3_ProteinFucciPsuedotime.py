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
plt.scatter(mean_mean_comp, mean_diff_from_rng)
plt.xlabel("Mean Mean Intensity")
plt.ylabel("Mean Additional Percent Variance Explained than Random")
plt.savefig(f"figures/IntensityVsMeanDiff.png")
plt.show()
plt.close()
    
#%% Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
# Idea: create cutoffs for percent variance and 
# Execution: create cutoffs for perc_var and total variance per compartment and for integrated intensity and mean intensity
# Output: Graphs that illustrate the cutoffs (integrated, mean)
# Output: Overlap of the total variance cutoffs with the original filtering done manually

plt.figure(figsize=(10,10))
plt.scatter(gini_comp, perc_var_comp, c=-np.log10(wp_comp_kruskal_gaussccd_adj))
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

ccd_comp = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
nonccd_comp = ~ccd_comp

n_tot_variable = len(u_well_plates)

print(f"{n_tot_variable}: # total proteins showing variation")
print(f"{sum(wp_prev_ccd & wp_comp_ccd_use) / len(prev_ccd_ensg)}: fraction of previously annotated CCD genes called CCD")
print(f"{len([ensg for ensg in wp_ensg[wp_comp_ccd_use] if ensg in ensggg1 and ensg in prev_ccd_ensg]) / len(wp_ensg[wp_comp_ccd_use])}: fraction of CCD genes that are previously annotated CCD genes")
print(f"{sum(ccd_comp)}: CCD variable proteins (before addressing redundancy and mitotic structures)")
print(f"{len(prev_ccd_ensg)}: CCD variable proteins (previous)")
print(f"{sum(nonccd_comp)}: non-CCD variable proteins")

### address gene redundancy
wp_ccd_unibimodal = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
wp_ccd_bimodalonecluster = wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2
wp_ccd_bimodaltwocluster = wp_comp_ccd_clust1 & wp_comp_ccd_clust2
wp_ensg_counts = np.array([sum([eeee == ensg for eeee in wp_ensg]) for ensg in wp_ensg])
ensg_is_duplicated = wp_ensg_counts > 1
duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
duplicated_ensg_pairs = [u_well_plates[wp_ensg == ensg] for ensg in duplicated_ensg]
print(f"{sum(wp_ccd_unibimodal[~ensg_is_duplicated])}: number of CCD proteins (no replicate, unimodal and bimodal)")
print(f"{sum(~wp_ccd_unibimodal[~ensg_is_duplicated])}: number of non-CCD proteins (no replicate, unimodal and bimodal)")
print(f"{sum(wp_ccd_bimodalonecluster[~ensg_is_duplicated])}: bimodal samples with one CCD cluster ({sum((wp_ccd_bimodalonecluster & wp_comp_ccd_difffromrng)[~ensg_is_duplicated])}: also CCD unimodally), no replicate")
print(f"{sum(wp_ccd_bimodaltwocluster[~ensg_is_duplicated])}: bimodal samples with two CCD clusters ({sum((wp_ccd_bimodaltwocluster & wp_comp_ccd_difffromrng)[~ensg_is_duplicated])}: also CCD unimodally), no replicate")
duplicated_ensg_ccd = np.array([sum(wp_ccd_unibimodal[wp_ensg == ensg]) for ensg in duplicated_ensg])
print(f"{sum(duplicated_ensg_ccd == 2)}: number of replicated stainings shown to be CCD in both replicates, unimodal and bimodal")
print(f"{sum(duplicated_ensg_ccd == 1)}: number of replicated stainings shown to be CCD in just one replicate, unimodal and bimodal")
print(f"{sum(duplicated_ensg_ccd == 0)}: number of replicated stainings shown to be non-CCD in both replicate, unimodal and bimodal")
duplicated_ensg_ccd_bi1 = np.array([sum(wp_ccd_bimodalonecluster[wp_ensg == ensg]) for ensg in duplicated_ensg])
duplicated_ensg_ccd_plusunimodal = np.array([sum((wp_comp_ccd_difffromrng & wp_ccd_bimodalonecluster)[wp_ensg == ensg]) for ensg in duplicated_ensg])
print(f"{sum(duplicated_ensg_ccd_bi1 == 2)}: number of replicated stainings shown to be CCD in both replicates, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 2)} also unimodally)")
print(f"{sum(duplicated_ensg_ccd_bi1 == 1)}: number of replicated stainings shown to be CCD in just one replicate, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 1)} also unimodally)")
print(f"{sum(duplicated_ensg_ccd_bi1 == 0)}: number of replicated stainings shown to be non-CCD in both replicate, bimodal in one cluster ({sum(duplicated_ensg_ccd_plusunimodal == 0)} also unimodally)")
duplicated_ensg_ccd_bi2 = np.array([sum(wp_ccd_bimodaltwocluster[wp_ensg == ensg]) for ensg in duplicated_ensg])
duplicated_ensg_ccd_plusunimodal = np.array([sum((wp_comp_ccd_difffromrng & wp_ccd_bimodaltwocluster)[wp_ensg == ensg]) for ensg in duplicated_ensg])
print(f"{sum(duplicated_ensg_ccd_bi2 == 2)}: number of replicated stainings shown to be CCD in both replicates, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 2)} also unimodally)")
print(f"{sum(duplicated_ensg_ccd_bi2 == 1)}: number of replicated stainings shown to be CCD in just one replicate, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 1)} also unimodally)")
print(f"{sum(duplicated_ensg_ccd_bi2 == 0)}: number of replicated stainings shown to be non-CCD in both replicate, bimodal in both clusters ({sum(duplicated_ensg_ccd_plusunimodal == 0)} also unimodally)")

bioccd = np.genfromtxt("input/processed/manual/biologically_defined_ccd.txt", dtype='str') # from mitotic structures
protein_ct = len(np.unique(np.concatenate((wp_ensg, bioccd))))
ccd_protein_ct = sum(wp_ccd_unibimodal[~ensg_is_duplicated]) + sum(duplicated_ensg_ccd == 2)
nonccd_protein_ct = sum(~wp_ccd_unibimodal[~ensg_is_duplicated]) + sum(duplicated_ensg_ccd == 0)
duplicated_ensg_bimodal_generally = np.array([sum(wp_isbimodal_generally[wp_ensg == ensg]) for ensg in duplicated_ensg])
unimodal_generally_protein_ct = sum(~wp_isbimodal_generally[~ensg_is_duplicated]) + sum(duplicated_ensg_bimodal_generally < 2)
bimodal_generally_protein_ct = sum(wp_isbimodal_generally[~ensg_is_duplicated]) + sum(duplicated_ensg_bimodal_generally == 2)
print(f"{ccd_protein_ct}: number of ccd proteins; addressed replicates; not including mitotic structures")

# Decision: remove proteins that didn't pass CCD in both replicates
ccd_comp[np.isin(wp_ensg, duplicated_ensg[duplicated_ensg_ccd == 1])] = False
nonccd_comp[np.isin(wp_ensg, duplicated_ensg[duplicated_ensg_ccd == 1])] = False

def analyze_replicates(wp_ccd_with_replicates, wp_ensg, analysis_tag):
    wp_ensg_counts = np.array([sum([eeee == ensg for eeee in wp_ensg]) for ensg in wp_ensg])
    ensg_is_duplicated = wp_ensg_counts > 1
    duplicated_ensg = np.unique(wp_ensg[ensg_is_duplicated])
    duplicated_ensg_ccd = np.array([sum(wp_ccd_with_replicates[wp_ensg == ensg]) for ensg in duplicated_ensg])
    print(f"{sum(wp_ccd_with_replicates[~ensg_is_duplicated])}: number of CCD proteins (no replicate, {analysis_tag})")
    print(f"{sum(~wp_ccd_with_replicates[~ensg_is_duplicated])}: number of non-CCD proteins (no replicate, {analysis_tag})")
    print(f"{sum(duplicated_ensg_ccd == 2)}: number of replicated stainings shown to be CCD in both replicates, {analysis_tag}")
    print(f"{sum(duplicated_ensg_ccd == 1)}: number of replicated stainings shown to be CCD in just one replicate, {analysis_tag}")
    print(f"{sum(duplicated_ensg_ccd == 0)}: number of replicated stainings shown to be non-CCD in both replicate,  {analysis_tag}")

analyze_replicates(wp_comp_ccd_gauss, wp_ensg, "gauss")
analyze_replicates(wp_comp_ccd_gauss & wp_comp_ccd_difffromrng, wp_ensg, "CCD and gaussian")
analyze_replicates(wp_comp_ccd_gauss & ~wp_comp_ccd_difffromrng, wp_ensg, "not CCD and gaussian")
analyze_replicates(~wp_comp_ccd_gauss & wp_comp_ccd_difffromrng, wp_ensg, "CCD and not gaussian")

pd.DataFrame({
        "ENSG":duplicated_ensg,
        "well_plate_pair":[",".join(pair) for pair in duplicated_ensg_pairs],
        "sum(ccd)":duplicated_ensg_ccd,
        "sum(bulk_ccd)":np.array([sum(wp_comp_ccd_gauss[wp_ensg == ensg]) for ensg in duplicated_ensg]),
        "ccd_pair":[",".join([str(wp_comp_ccd_use[u_well_plates == wp][0]) for wp in pair]) for pair in duplicated_ensg_pairs]
    }).to_csv("output/ReplicatedCellCycleDependentProteins.csv")

# Replicates
for ensg in fileprefixes[ensg_is_duplicated]:
    for file in glob.glob(os.path.join(folder, f"{ensg}*_mvavg.png")):
        shutil.copy(file, os.path.join(replicatefolder, os.path.basename(file)))
###

examples=np.loadtxt("input/processed/manual/ensgexamplesfrompaper.txt", dtype=str)
print(f"{sum(~np.isin(examples,wp_ensg[wp_comp_ccd_use]))}: number of examples missing of {len(examples)}, which includes 1 negative control")
print(examples[~np.isin(examples,wp_ensg[wp_comp_ccd_use])])

# Accounting for biolocially CCD ones
knownccd1 = np.genfromtxt("input/processed/manual/knownccd.txt", dtype='str') # from gene ontology, reactome, cyclebase 3.0, NCBI gene from mcm3
knownccd2 = np.genfromtxt("input/processed/manual/known_go_ccd.txt", dtype='str') # from GO cell cycle
knownccd3 = np.genfromtxt("input/processed/manual/known_go_proliferation.txt", dtype='str') # from GO proliferation
print(f"{len(bioccd)}: number of mitotic structure proteins")

ccd_prots_withmitotic = np.unique(np.concatenate((wp_ensg[ccd_comp], bioccd)))
total_proteins_minusmitotic = len(ccd_prots_withmitotic) - sum(np.isin(wp_ensg[nonccd_comp], bioccd))
total_ccd_proteins_withmitotic = len(ccd_prots_withmitotic)
total_nonccd_proteins_minusmitotic = nonccd_protein_ct - sum(np.isin(wp_ensg[nonccd_comp], bioccd))

overlapping_knownccd1 = sum(np.isin(ccd_prots_withmitotic, np.concatenate((knownccd1, knownccd2, knownccd3))))
overlapping_knownccd2 = sum(np.isin(ccd_prots_withmitotic, knownccd2))
overlapping_knownccd3 = sum(np.isin(ccd_prots_withmitotic, knownccd3))

plt.figure(figsize=(10,10))
plt.scatter(gini_comp, perc_var_comp, c=-np.log10(wp_comp_kruskal_gaussccd_adj))
for iii in np.nonzero(np.isin(wp_ensg, examples))[0]:
    print(iii)
    plt.annotate(wp_ensg[iii], (gini_comp[iii], perc_var_comp[iii]))
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

ccdstring = np.array(["No                 "] * len(ccd_comp))
ccdstring[ccd_comp] = "Pseudotime"
ccdstring[np.isin(wp_ensg, bioccd)] = "Mitotic"
ccdstring[ccd_comp & np.isin(wp_ensg, bioccd)] = "Pseudotime&Mitotic"
pd.DataFrame({
    "well_plate" : u_well_plates, 
    "ENSG": wp_ensg,
    "antibody": wp_ab,
    "variance_comp":var_comp,
    "gini_comp":gini_comp,
    "known_by_GoReactomeCyclebaseNcbi":np.isin(wp_ensg, np.concatenate((knownccd1, knownccd2, knownccd3))),
    "mean_percvar_diff_from_random":mean_diff_from_rng,
    "wp_comp_kruskal_gaussccd_adj":wp_comp_kruskal_gaussccd_adj,
    "log_10_pval_eq_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj, wp_comp_eq_percvar_adj+1)),
    "pass_median_diff":wp_comp_ccd_difffromrng,
    "pass_gauss":wp_comp_ccd_gauss,
    "CCD_COMP":ccd_comp,
    "ccd_reason":ccdstring,
    "nonccd_comp":nonccd_comp,
    "wp_prev_ccd":wp_prev_ccd,
    
    # bimodal significance testing
    "ccd_unimodal":wp_comp_ccd_difffromrng,
    "ccd_clust1":wp_comp_ccd_clust1,
    "clust1_difffromrng":mean_diff_from_rng_clust1,
    "clust1_log10pval_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj_clust1, wp_comp_eq_percvar_adj_clust1+1)),
    "ccd_clust2":wp_comp_ccd_clust2,
    "ccd_clust2_difffromrng":mean_diff_from_rng_clust2,
    "ccd_clust2_log10pval_percvar":-np.log10(np.nextafter(wp_comp_eq_percvar_adj_clust2, wp_comp_eq_percvar_adj_clust2+1)),
    }).to_csv("output/CellCycleVariationSummary.csv", index=False)
    
pd.DataFrame({"ENSG":np.unique(wp_ensg[wp_prev_ccd & ~ccd_comp])}).to_csv("output/DianaCCDMissingProteins.csv")
pd.DataFrame({"ENSG":np.unique(wp_ensg[~wp_prev_ccd & ccd_comp])}).to_csv("output/NewToDianaCCDProteins.csv") 

#%% Pickle the results needed later
np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp) # removed ones passing in only one replicate
np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp) # removed ones passing in only one replicate
np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)

pd.DataFrame({"gene": wp_ensg[ccd_comp]}).to_csv("output/picklestxt/ccd_compartment_ensg.txt", index=False)
pd.DataFrame({"gene": wp_ensg[nonccd_comp]}).to_csv("output/picklestxt/nonccd_compartment_ensg.txt", index=False)

#%% Figures of merit
with open("output/figuresofmerit.txt", "w") as file:
    fom = "--- protein pseudotime\n\n"
    fom += f"{len(np.unique(wp_ensg))} proteins that were expressed and exhibited variations in the U-2 OS cell line were selected" + "\n\n"
    fom += f"present the first evidence of cell cycle association for {overlapping_knownccd1} proteins" + "\n\n"
    fom += f"Based on this analysis, we identified {ccd_protein_ct} out of {len(np.unique(wp_ensg))} proteins ({100 * ccd_protein_ct / len(np.unique(wp_ensg))}%) to have variance in expression levels temporally correlated to cell cycle progression, and for which the cell-cycle explained 8% or more variance in expression than random." + "\n\n"
    fom += f"majority of the proteins analyzed ({100 * total_nonccd_proteins_minusmitotic / len(np.unique(wp_ensg))}%) showed cell-to-cell variations that were largely unexplained by cell cycle progression" + "\n\n"
    fom += f"Of the {total_ccd_proteins_withmitotic} proteins ({ccd_protein_ct} in interphase, {len(bioccd)} in mitotic structures, and {-(total_ccd_proteins_withmitotic - ccd_protein_ct - len(bioccd))} in both sets) identified to correlate to cell cycle progression, {overlapping_knownccd1} ({100 * overlapping_knownccd1 / total_ccd_proteins_withmitotic}%) had a known association to the cell cycle as determined either by a GO BP term ... The remaining {total_ccd_proteins_withmitotic - overlapping_knownccd1} proteins ({100* (total_ccd_proteins_withmitotic - overlapping_knownccd1) / total_ccd_proteins_withmitotic}%)," + "\n\n"
    fom += f"The patterns of variability were investigated for these {len(wp_ensg)} proteins for the population of cells measured for each protein. The mean fold change between the highest and lowest expressing cells per protein was {np.mean(wp_bimodal_fcmaxmin)}." + "\n\n"
    fom += f"We determined that {unimodal_generally_protein_ct} proteins ({100 * unimodal_generally_protein_ct / len(np.unique(wp_ensg))}%) had unimodal intensity distributions, and {bimodal_generally_protein_ct} proteins ({100 * bimodal_generally_protein_ct / len(np.unique(wp_ensg))}%) were found to display bimodality" + "\n\n"
    fom += f"Of {sum(wp_isbimodal_fcpadj_pass)} bimodal samples that were analyzed for cell cycle dependence, {sum(wp_ccd_bimodalonecluster)} were CCD in one cluster ({sum(wp_ccd_bimodalonecluster & wp_comp_ccd_difffromrng)} of these were CCD when analyzed unimodally), and {sum(wp_ccd_bimodaltwocluster)} were CCD in both clusters ({sum(wp_ccd_bimodaltwocluster & wp_comp_ccd_difffromrng)} were also CCD when analyzed unimodally), and the remaining {sum(wp_isbimodal_fcpadj_pass & ~wp_ccd_bimodalonecluster & ~wp_ccd_bimodaltwocluster)} were non-CCD in both clusters."
    print(fom)
    file.write(fom)