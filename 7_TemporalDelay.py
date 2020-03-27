#%% Imports
from utils import *
import numpy as np
import collections
import fucci_plotting
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sbn
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'] = 42, 42 #Make PDF text readable

#%% Read in RNA-Seq data again and the CCD gene lists
dd = "All"
count_or_rpkm = "Tpms" # so that the gene-specific results scales match for cross-gene comparisons
print("reading scRNA-Seq data")
biotype_to_use="protein_coding"
adata, phases = read_counts_and_phases(dd, count_or_rpkm, use_spike_ins=False, biotype_to_use=biotype_to_use)
adata, phasesfilt = qc_filtering(adata, do_log_normalize=False, do_remove_blob=True)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
u_plate = np.load("output/pickles/u_plate.npy", allow_pickle=True)
well_plate=np.load("output/pickles/well_plate.npy", allow_pickle=True)
u_well_plates=np.load("output/pickles/u_well_plates.npy", allow_pickle=True)
ab_nuc=np.load("output/pickles/ab_nuc.npy", allow_pickle=True)
ab_cyto=np.load("output/pickles/ab_cyto.npy", allow_pickle=True)
ab_cell=np.load("output/pickles/ab_cell.npy", allow_pickle=True)
mt_cell=np.load("output/pickles/mt_cell.npy", allow_pickle=True)
green_fucci=np.load("output/pickles/green_fucci.npy", allow_pickle=True)
red_fucci=np.load("output/pickles/red_fucci.npy", allow_pickle=True)
log_red_fucci_zeroc=log_green_fucci_zeroc=np.load("output/pickles/log_green_fucci_zeroc.npy", allow_pickle=True)
log_red_fucci_zeroc=np.load("output/pickles/log_red_fucci_zeroc.npy", allow_pickle=True)
log_green_fucci_zeroc_rescale=np.load("output/pickles/log_green_fucci_zeroc_rescale.npy", allow_pickle=True)
log_red_fucci_zeroc_rescale=np.load("output/pickles/log_red_fucci_zeroc_rescale.npy", allow_pickle=True)
fucci_data = np.load("output/pickles/fucci_data.npy", allow_pickle=True)
pol_sort_well_plate = np.load("output/pickles/pol_sort_well_plate.npy", allow_pickle= True)
pol_sort_norm_rev = np.load("output/pickles/pol_sort_norm_rev.npy", allow_pickle=True)
pol_sort_ab_nuc = np.load("output/pickles/pol_sort_ab_nuc.npy", allow_pickle=True)
pol_sort_ab_cyto = np.load("output/pickles/pol_sort_ab_cyto.npy", allow_pickle=True)
pol_sort_ab_cell = np.load("output/pickles/pol_sort_ab_cell.npy", allow_pickle=True)
pol_sort_mt_cell = np.load("output/pickles/pol_sort_mt_cell.npy", allow_pickle=True)
#pol_sort_fred = np.load("output/pickles/pol_sort_fred.npy", allow_pickle=True)
#pol_sort_fgreen = np.load("output/pickles/pol_sort_fgreen.npy", allow_pickle=True)
wp_ensg = np.load("output/pickles/wp_ensg.npy", allow_pickle=True)
wp_iscell = np.load("output/pickles/wp_iscell.npy", allow_pickle=True)
wp_isnuc = np.load("output/pickles/wp_isnuc.npy", allow_pickle=True)
wp_iscyto = np.load("output/pickles/wp_iscyto.npy", allow_pickle=True)
ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
print("loaded")

print("reading RNA data")
ccdtranscript = np.load("output/pickles/ccdtranscript.npy", allow_pickle=True)
ccdprotein_transcript_regulated = np.load("output/pickles/ccdprotein_transcript_regulated.npy", allow_pickle=True)
ccdprotein_nontranscript_regulated = np.load("output/pickles/ccdprotein_nontranscript_regulated.npy", allow_pickle=True)
print("loaded")

#%% Make temporal heatmap and use those peak values to compare to RNA data rank by percent variance explained.
# Idea: make a heatmap showing the relative expression of each protein over the cell cycle, sorted by time of peak expression
# Execution: plt.imshow makes a heatmap if given a 2D array
# Output: heatmap

###PLOTTING
# temporal heatmap
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

def fix_nans(binned_values):
    for i,val in enumerate(binned_values):
        do_fix = np.isnan(val)
        if do_fix:
            is_end = i == len(binned_values)-1
            if is_end:
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

def bin_values(nbins):
    xvals = np.linspace(0,1,num=nbins)
    wp_max_pol = []
    wp_binned_values = []
    for i, well in enumerate(u_well_plates):
        curr_well_inds = pol_sort_well_plate==well
        curr_pol = pol_sort_norm_rev[curr_well_inds]
    #    curr_fred = pol_sort_fred[curr_well_inds]
    #    curr_fgreen = pol_sort_fgreen[curr_well_inds]
        curr_ab_cell, curr_ab_nuc, curr_ab_cyto, curr_mt_cell = pol_sort_ab_cell[curr_well_inds], pol_sort_ab_nuc[curr_well_inds],pol_sort_ab_cyto[curr_well_inds], pol_sort_mt_cell[curr_well_inds]
    
        # Normalize FUCCI colors & mean intensities, normalized for display
    #    curr_fred_norm = curr_fred / np.max(curr_fred)
    #    curr_fgreen_norm = curr_fgreen / np.max(curr_fgreen)
        curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell) 
        curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
        curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto) 
        curr_mt_cell_norm  = curr_mt_cell / np.max(curr_mt_cell)
        
        # Compute binned values
        binned_values = []
        prev_xval = 0
        for xval in xvals:
            if xval==0:
                prev_xval = xval
                continue
            curr_ab_norm = curr_ab_cell_norm if wp_iscell[i] else curr_ab_nuc_norm if wp_isnuc[i] else curr_ab_cyto_norm
            binned_values.append(np.median(curr_ab_norm[(curr_pol < xval) & (curr_pol >= prev_xval)]))
            prev_xval = xval
            
        binned_values = binned_values/np.nanmax(binned_values)
        binned_values = fix_nans(binned_values)
        max_loc = np.nanargmax(binned_values)
        if np.isnan(xvals[max_loc]): print('what')
        
        wp_max_pol.append(xvals[max_loc])
        wp_binned_values.append(binned_values)
    return wp_max_pol, wp_binned_values

# Make an expression array with the CCD proteins
wp_max_pol, wp_binned_values = bin_values(20)
wp_max_pol, wp_binned_values = np.array(wp_max_pol), np.array(wp_binned_values)
wp_max_pol_ccd, wp_binned_values_ccd = wp_max_pol[ccd_comp], wp_binned_values[ccd_comp]
wp_max_sort_inds = np.argsort(wp_max_pol_ccd)
sorted_gene_array = np.take(wp_binned_values_ccd, wp_max_sort_inds, axis=0) # this is the expression values, binned_values, sorted by the binned value at max location (can do in the temporal part)
sorted_maxpol_array = np.take(wp_max_pol_ccd, wp_max_sort_inds)

# Actually making the figure
fig, ax = plt.subplots(figsize=(10, 10))
sc = ax.imshow(sorted_gene_array, interpolation='nearest')

# Do the x ticks
xtick_labels = [str(np.around(x * TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)] #+['G1/S','S/G2']
my_xticks = np.arange(-.5, 20, 2)
num_ticks = 20
xphase_labels = ['G1/S','S/G2']
phase_trans = np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])
ax.set_xticks(my_xticks,minor=True)
ax.set_xticklabels(xtick_labels,minor=True)
ax.set_xticks(phase_trans, minor=False)
ax.set_xticklabels(xphase_labels, minor=False)
ax.tick_params(length=12)

#Do the y ticks
#ytick_locs = [prot_list.index(p) for p in HIGHLIGHTS]
## ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
#ax.set_yticks(ytick_locs,minor=False)
#ax.set_yticklabels(HIGHLIGHTS,minor=False)
## ax.set_yticks(ytick_locs_minor,minor=True)
## ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)

ax.tick_params(direction='out', length=12, width=2, colors='k', axis='x',which='major')
ax.tick_params(direction='out', length=56, width=1, colors='k', axis='y',which='major')
ax.set_aspect('auto')
plt.xlabel('Division Cycle, hrs',size=20,fontname='Arial')
plt.ylabel('Gene',size=20,fontname='Arial')
plt.xticks(size=12,fontname='Arial')
plt.yticks(size=10,fontname='Arial')
divider1 = make_axes_locatable(ax)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(sc, cax = cax1)
cbar.set_label('Relative expression', fontname='Arial', size=20)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(os.path.join("figures",'sorted_heatmap21_sw30_take4.pdf'), transparent=True)
plt.savefig(os.path.join("figures",'sorted_heatmap21_sw30_take4.png'), transparent=True)
plt.show()

np.save("output/pickles/wp_max_pol.npy", wp_max_pol, allow_pickle=True)

#%% Correlations of related genes
plt.rcParams['figure.figsize'] = (5, 5)
plt.rcParams['figure.fontname'] = "Arial"

def scatter_genes(gene1, gene2, r):
    plt.scatter(wp_binned_values[wp_ensg == gene1[0]][0], wp_binned_values[wp_ensg == gene2[0]][0])
    plt.xlabel(f"{gene1[1]} Expression Binned by Pseudotime", fontsize=14, fontname="Arial")
    plt.ylabel(f"{gene2[1]} Expression Binned by Pseudotime", fontsize=14, fontname="Arial")
    plt.text(np.min(wp_binned_values[wp_ensg == gene1[0]][0]), np.max(wp_binned_values[wp_ensg == gene2[0]][0]), f"Pearson's r = {r}", fontsize=14, fontname="Arial")
    plt.savefig(f"figures/Correlations/{gene1[1]}_{gene2[1]}.pdf")
    plt.show()
    plt.close()

for ensg in [("ENSG00000091651", "orc6"),
             ("ENSG00000169740", "znf32"),
             ("ENSG00000105173", "ccne1"),
             ("ENSG00000162999", "dusp19"),
             ("ENSG00000123607", "ttc21b"),
             ("ENSG00000173599", "pc"),
             ("ENSG00000134057", "ccnb1"),#, known
             ("ENSG00000178999", "aurkb"),#, known
             ("ENSG00000156970", "bub1b"),#, known
             ("ENSG00000167065", "dusp18"),#, unknown
             ("ENSG00000138801", "papss1"),#, unknown
             ("ENSG00000156239", "n6amt1"),#, unknown
             ("ENSG00000019144", "phldb1"),#, unknown
             ("ENSG00000151702", "fli1"),#, unknown
             ("ENSG00000132768", "dph2"),#, unknown
             ("ENSG00000102908", "nfat5") # unknown
             ]: 
    print(ensg)
    print(f"number of observations: {sum(pol_sort_well_plate==u_well_plates[wp_ensg == ensg[0]][0])}")
    print(f"time of peak expression: {wp_max_pol[wp_ensg == ensg[0]][0]}")

orc6_znf32 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000091651'][0], wp_binned_values[wp_ensg == 'ENSG00000169740'][0])
scatter_genes(("ENSG00000091651", "ORC6"), ("ENSG00000169740", "ZNF32"), round(orc6_znf32[0], 2))

ccne1_dusp19 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000105173'][0], wp_binned_values[wp_ensg == 'ENSG00000162999'][0])
scatter_genes(("ENSG00000105173", "CCNE1"), ("ENSG00000162999", "DUSP19"), round(ccne1_dusp19[0], 2))

ttc21b_pc = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000123607'][0], wp_binned_values[wp_ensg == 'ENSG00000173599'][0])
scatter_genes(("ENSG00000123607", "TTC21B"), ("ENSG00000173599", "PC"), round(bub1b_dusp18[0], 2))

bub1b_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000156970'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
scatter_genes(("ENSG00000156970", "BUB1B"), ("ENSG00000167065", "DUSP18"), round(bub1b_dusp18[0], 2))

aurkb_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000178999'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
scatter_genes(("ENSG00000178999", "AURKB"), ("ENSG00000167065", "DUSP18"), round(aurkb_dusp18[0], 2))

ccnb1_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000167065", "DUSP18"), round(ccnb1_dusp18[0], 2))

ccnb1_papss1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000138801'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000138801", "PAPSS1"), round(ccnb1_papss1[0], 2))

ccnb1_n6amt1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000156239'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000156239", "N6AMT1"), round(ccnb1_n6amt1[0], 2))

ccnb1_phldb1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000019144'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000019144", "PHLDB1"), round(ccnb1_phldb1[0], 2))

ccnb1_fli1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000151702'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000151702", "FLI1"), round(ccnb1_fli1[0], 2))

ccnb1_dph2 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000132768'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000132768", "DPH2"), round(ccnb1_dph2[0], 2))

ccnb1_nfat5 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000102908'][0])
scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000102908", "NFAT5"), round(ccnb1_nfat5[0], 2))

print(f"correlation of ORC6 and ZNF32: {orc6_znf32[0]}")
print(f"correlation of CCNE1 and DUSP19: {ccne1_dusp19[0]}")
print(f"correlation of TTC21B and PC: {ttc21b_pc[0]}")
print()
print(f"correlation of bub1b_dusp18: {bub1b_dusp18[0]}")
print(f"correlation of aurkb_dusp18: {aurkb_dusp18[0]}")
print(f"correlation of ccnb1_papss1: {ccnb1_papss1[0]}")
print(f"correlation of ccnb1_n6amt1: {ccnb1_n6amt1[0]}")
print(f"correlation of ccnb1_phldb1: {ccnb1_phldb1[0]}")
print(f"correlation of ccnb1_fli1: {ccnb1_fli1[0]}")
print(f"correlation of ccnb1_dph2: {ccnb1_dph2[0]}")

#%% RNA heatmap
# Idea: create a heatmap of peak RNA expression
# Execution: use the same sort of methods as the previous cell
# Output: heatmap
def mvavg(yvals, mv_window):
    return np.convolve(yvals, np.ones((mv_window,))/mv_window, mode='valid')

def binned_median(yvals, nbins):
    binned_medians = []
    for xval in range(nbins):
        startidx = len(yvals) // nbins * xval
        endidx = len(yvals) // nbins * (xval + 1)
        binned_medians.append(np.median(yvals[startidx:endidx]))
    return binned_medians

# Get the peak RNA expression polar locations
expression_data = adata.X # log normalized
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
moving_averages = np.apply_along_axis(mvavg, 0, norm_exp_sort, 100)
max_moving_avg_loc = np.argmax(moving_averages, 0)
max_moving_avg_pol = np.take(fucci_time_sort, max_moving_avg_loc)
max_moving_avg_pol_ccd = max_moving_avg_pol[ccdtranscript]
moving_averages_ccd = moving_averages[:,ccdtranscript]
max_moving_avg_pol_sortinds = np.argsort(max_moving_avg_pol_ccd)
sorted_max_moving_avg_pol_ccd = np.take(max_moving_avg_pol_ccd, max_moving_avg_pol_sortinds)
sorted_rna_array = np.take(moving_averages_ccd, max_moving_avg_pol_sortinds, axis=1).T
sorted_rna_binned = np.apply_along_axis(binned_median, 1, sorted_rna_array, len(xvals))
sorted_rna_binned_norm = sorted_rna_binned / np.max(sorted_rna_binned, axis=1)[:,None]

fig, ax = plt.subplots(figsize=(10, 10))
sc = ax.imshow(sorted_rna_binned_norm, interpolation='nearest')

# Do the x ticks
xtick_labels = [str(np.around(x * TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)] #+['G1/S','S/G2']
my_xticks = np.arange(-.5, 20, 2)
num_ticks = 20
xphase_labels = ['G1/S','S/G2']
phase_trans = np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])
ax.set_xticks(my_xticks,minor=True)
ax.set_xticklabels(xtick_labels,minor=True)
ax.set_xticks(phase_trans, minor=False)
ax.set_xticklabels(xphase_labels, minor=False)
ax.tick_params(length=12)

#Do the y ticks
#ytick_locs = [prot_list.index(p) for p in HIGHLIGHTS]
## ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
#ax.set_yticks(ytick_locs,minor=False)
#ax.set_yticklabels(HIGHLIGHTS,minor=False)
## ax.set_yticks(ytick_locs_minor,minor=True)
## ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)

ax.tick_params(direction='out', length=12, width=2, colors='k', axis='x',which='major')
ax.tick_params(direction='out', length=56, width=1, colors='k', axis='y',which='major')
ax.set_aspect('auto')
plt.xlabel('Division Cycle, hrs',size=20,fontname='Arial')
plt.ylabel('Gene',size=20,fontname='Arial')
plt.xticks(size=12,fontname='Arial')
plt.yticks(size=10,fontname='Arial')
divider1 = make_axes_locatable(ax)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(sc, cax = cax1)
cbar.set_label('Relative expression', fontname='Arial', size=20)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(os.path.join("figures",'sorted_rna_heatmap.pdf'), transparent=True)
plt.savefig(os.path.join("figures",'sorted_rna_heatmap.png'), transparent=True)
plt.show()


#%%
# Idea: calculate the peak RNA expression and compare to the peak protein expression for each gene
# Execution: compare distribution of differences between peaks; grainger test, too
# Output: plot of dist of differences; grainger test results
prot_ccd_ensg = list(wp_ensg[ccd_comp])
rna_ccd_ensg = list(adata.var_names[ccdtranscript])
both_ccd_ensg = np.intersect1d(prot_ccd_ensg, rna_ccd_ensg)
both_prot_ccd_idx = np.array([prot_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
both_rna_ccd_idx = np.array([rna_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
insct_prot_max_pol_ccd = wp_max_pol_ccd[both_prot_ccd_idx]
insct_rna_max_pol_ccd = sorted_max_moving_avg_pol_ccd[both_rna_ccd_idx]
diff_max_pol = insct_prot_max_pol_ccd - insct_rna_max_pol_ccd

print(f"length of prot CCD genes: {len(prot_ccd_ensg)}")
print(f"length of CCD RNA genes: {len(rna_ccd_ensg)}")
print(f"length of intersection betweeen CCD prot and CCD RNA: {len(both_ccd_ensg)}")

# Make link files for Circos using labels g1 [0, 109], g1s [0, 27], g2 [0, 120]
with open("output/differencelinks.txt", "w") as file:
    for iii, diff in enumerate(list(diff_max_pol)):
        name = both_ccd_ensg[iii]
        rnatime = insct_rna_max_pol_ccd[iii]
        prottime = insct_prot_max_pol_ccd[iii]
        startloc = "g1" if rnatime * TOT_LEN < G1_LEN else "g1s" if rnatime * TOT_LEN >= G1_LEN and rnatime * TOT_LEN < G1_LEN + G1_S_TRANS else "g2"
        starttime = rnatime * TOT_LEN * 10 if startloc == "g1" else (rnatime * TOT_LEN - G1_LEN) * 10 if startloc == "g1s" else (rnatime * TOT_LEN - G1_LEN - G1_S_TRANS) * 10
        endloc = "g1" if prottime * TOT_LEN < G1_LEN else "g1s" if prottime * TOT_LEN >= G1_LEN and prottime * TOT_LEN < G1_LEN + G1_S_TRANS else "g2"
        endtime =  prottime * TOT_LEN * 10 if endloc == "g1" else (prottime * TOT_LEN - G1_LEN) * 10 if endloc == "g1s" else (prottime * TOT_LEN - G1_LEN - G1_S_TRANS) * 10
        line = [startloc, str(int(starttime)), str(int(starttime+1)), endloc, str(int(endtime)), str(int(endtime+1))]
        file.write("\t".join(line) + "\n")
        
with open("output/differencelinkends.txt", "w") as file:
      for iii, diff in enumerate(list(diff_max_pol)):
        name = both_ccd_ensg[iii]
        prottime = insct_prot_max_pol_ccd[iii]
        endloc = "g1" if prottime * TOT_LEN < G1_LEN else "g1s" if prottime * TOT_LEN >= G1_LEN and prottime * TOT_LEN < G1_LEN + G1_S_TRANS else "g2"
        endtime =  prottime * TOT_LEN * 10 if endloc == "g1" else (prottime * TOT_LEN - G1_LEN) * 10 if endloc == "g1s" else (prottime * TOT_LEN - G1_LEN - G1_S_TRANS) * 10
        line = [endloc, str(int(endtime)), str(int(endtime+1)), str(0)]
        file.write("\t".join(line) + "\n")  
        
# Make a sankey plot for the transitions
# colordict = {"g1" : 'r', "g1s": 'y', "g2": 'g'}
# transitions = pd.DataFrame({
#         "source": ["g1" if rnatime * TOT_LEN < G1_LEN else "g1s" if rnatime * TOT_LEN >= G1_LEN and rnatime * TOT_LEN < G1_LEN + G1_S_TRANS else "g2" for rnatime in insct_rna_max_pol_ccd],
#         "target": ["g1" if prottime * TOT_LEN < G1_LEN else "g1s" if prottime * TOT_LEN >= G1_LEN and prottime * TOT_LEN < G1_LEN + G1_S_TRANS else "g2" for prottime in insct_prot_max_pol_ccd]})
import alluvial
transitiondict = {}
startphases = ["G1\nRNA\nPeak", "S\nRNA\nPeak", "G2\nRNA\nPeak"]
endphases= ["G1\nProtein\nPeak", "S\nProtein\nPeak", "G2\nProtein\nPeak"]
for iii, diff in enumerate(list(diff_max_pol)):
    rnatime = insct_rna_max_pol_ccd[iii]
    prottime = insct_prot_max_pol_ccd[iii]
    startloc = startphases[0] if rnatime * TOT_LEN < G1_LEN else startphases[1] if rnatime * TOT_LEN >= G1_LEN and rnatime * TOT_LEN < G1_LEN + G1_S_TRANS else startphases[2]
    endloc = endphases[0] if prottime * TOT_LEN < G1_LEN else endphases[1] if prottime * TOT_LEN >= G1_LEN and prottime * TOT_LEN < G1_LEN + G1_S_TRANS else endphases[2]
    if startloc in transitiondict and endloc in transitiondict[startloc]:
        transitiondict[startloc][endloc] += 1
    elif startloc in transitiondict:
        transitiondict[startloc][endloc] = 1
    else:
        transitiondict[startloc] = {endloc:1}
cmap = plt.cm.get_cmap('viridis_r')
ax = alluvial.plot(transitiondict, colors=['r', 'y', 'b'], a_sort=startphases, b_sort=endphases)
fig = ax.get_figure()
fig.set_size_inches(5 ,10)
plt.savefig("figures/transitions.png")
plt.savefig("figures/transitions.pdf")
plt.show()

plt.hist(diff_max_pol * TOT_LEN, alpha=0.5)# kde=False, rug=True, fit=scipy.stats.gamma)
plt.xlabel("Delay in peak protein expression from peak RNA expression, hrs")
plt.ylabel("Count of CCD Proteins")
plt.tight_layout()
plt.savefig("figures/DelayPeakProteinRNA.pdf")
plt.show()
plt.close()

f,ax = plt.subplots(figsize=(6,5))
plt.scatter(x=insct_rna_max_pol_ccd * TOT_LEN, y=insct_prot_max_pol_ccd * TOT_LEN, c=diff_max_pol * TOT_LEN)
# ax.hist(diff_max_pol * TOT_LEN, orientation="vertical", alpha=0.5)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set_xticks([0,5,10,15,20])
cbar = plt.colorbar()
cbar.set_label('Temporal Delay, hrs',fontname='Arial',size=20)
cbar.ax.tick_params(labelsize=18)
plt.xlabel("Peak RNA Expression, hrs",fontname='Arial',size=20)
plt.ylabel("Peak Protein Expression, hrs",fontname='Arial',size=20)
plt.tight_layout()
plt.savefig("figures/TemporalDelayScatter.pdf")
plt.show()
plt.close()

# sbn.jointplot(x=insct_rna_max_pol_ccd * TOT_LEN, y=insct_prot_max_pol_ccd * TOT_LEN, kind="kde", color="k")
# plt.close()

# f,ax = plt.subplots(figsize=(6,6))
# sbn.kdeplot(insct_rna_max_pol_ccd * TOT_LEN, insct_prot_max_pol_ccd * TOT_LEN, ax=ax, color="b")
# sbn.distplot(insct_rna_max_pol_ccd * TOT_LEN, color="orange", kde=False, ax=ax)
# sbn.distplot(insct_prot_max_pol_ccd * TOT_LEN, color="blue", kde=False, vertical=True, ax=ax)
# ax.set_xlim(-4)
# plt.close()

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.hist(insct_prot_max_pol_ccd * TOT_LEN, alpha=0.5, label="Peak Protein Expression Time, hrs")
# ax1.hist(insct_rna_max_pol_ccd * TOT_LEN, alpha=0.5, label="Peak RNA Expression Time, hrs")
# plt.legend(loc="upper left")
# plt.xlabel("Division Cycle, hrs")
# plt.ylabel("Count of Cell Cycle Genes")
# plt.tight_layout()
# plt.savefig(f"figures/DelayPeakProteinRNA_separate.png")
# plt.show()
# plt.close()

mmmm = np.concatenate((insct_prot_max_pol_ccd * TOT_LEN, insct_rna_max_pol_ccd * TOT_LEN))
cccc = (["Protein"] * len(insct_prot_max_pol_ccd))
cccc.extend(["RNA"] * len(insct_rna_max_pol_ccd))
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
print(f"Median RNA expression time for CCD proteins: {TOT_LEN * np.median(insct_rna_max_pol_ccd)}")
print(f"Median protein expression time for CCD proteins: {TOT_LEN * np.median(insct_prot_max_pol_ccd)}")
t, p = scipy.stats.kruskal(insct_rna_max_pol_ccd, insct_prot_max_pol_ccd)
print(f"One-sided kruskal for median protein expression time higher than median RNA expression time: {2*p}")
t, p = scipy.stats.ttest_1samp(diff_max_pol, 0)
print(f"One-sided, one-sample t-test for mean delay in protein expression larger than zero: {2*p}")

#%% Plot the variances against each other
total_variance_rna = np.var(norm_exp_sort, 0)
total_gini_rna = np.apply_along_axis(gini, 0, norm_exp_sort)
total_cv_rna = np.apply_along_axis(scipy.stats.variation, 0, norm_exp_sort)
var_comp_prot = np.load("output/pickles/var_comp.npy", allow_pickle=True)
gini_comp_prot = np.load("output/pickles/gini_comp.npy", allow_pickle=True)
cv_comp_prot = np.load("output/pickles/cv_comp.npy", allow_pickle=True)
var_cell_prot = np.load("output/pickles/var_cell.npy", allow_pickle=True)
gini_cell_prot = np.load("output/pickles/gini_cell.npy", allow_pickle=True)
cv_cell_prot = np.load("output/pickles/cv_cell.npy", allow_pickle=True)

prot_ensg = list(wp_ensg)
rna_ensg = list(adata.var_names)
both_ensg = np.intersect1d(prot_ensg, rna_ensg)
both_prot_idx = np.array([prot_ensg.index(ensg) for ensg in both_ensg])
both_rna_idx = np.array([rna_ensg.index(ensg) for ensg in both_ensg])
insct_prot_variance = var_comp_prot[both_prot_idx]
insct_prot_gini = gini_comp_prot[both_prot_idx]
insct_prot_cv = cv_comp_prot[both_prot_idx]
insct_prot_variance_cell = var_cell_prot[both_prot_idx]
insct_prot_gini_cell = gini_cell_prot[both_prot_idx]
insct_prot_cv_cell = cv_cell_prot[both_prot_idx]
insct_rna_variance = total_variance_rna[both_rna_idx]
insct_rna_gini = total_gini_rna[both_rna_idx]
insct_rna_cv = total_cv_rna[both_rna_idx]

plt.scatter(insct_rna_variance, insct_prot_variance_cell)
plt.xlabel("Total RNA Variance")
plt.ylabel("Total Protein Variance")
plt.savefig("figures/ProteinRNAVariance.png")
plt.show()
plt.close()
plt.scatter(insct_rna_gini, insct_prot_gini_cell)
plt.xlabel("RNA Gini")
plt.ylabel("Protein Gini")
plt.savefig("figures/ProteinRNAGini.png")
plt.show()
plt.close()
plt.scatter(insct_rna_cv, insct_prot_cv_cell)
plt.xlabel("RNA CV")
plt.ylabel("Protein CV")
plt.savefig("figures/ProteinRNACV.png")
plt.show()
plt.close()

pd.DataFrame({
    "gene":both_ensg,
    "variance_rna": insct_rna_variance,
    "gini_rna":insct_rna_gini,
    "cv_rna":insct_rna_cv,
    "variance_comp_prot":insct_prot_variance,
    "gini_comp_prot":insct_prot_gini,
    "cv_comp_prot":insct_prot_cv,
    "variance_cell_prot":insct_prot_variance_cell,
    "gini_cell_prot":insct_prot_gini_cell,
    "cv_cell_prot":insct_prot_cv_cell,
    }).to_csv("output/VarianceRNAProtein.csv",index=False)

#%% Output tables
pd.DataFrame({"gene" : wp_ensg, "max_pol_protein": wp_max_pol, "max_time_protein": wp_max_pol * TOT_LEN}).to_csv("output/max_pol_protein.csv", index=False)
pd.DataFrame({"gene" : adata.var_names, "max_pol_rna": max_moving_avg_pol, "max_time_rna": max_moving_avg_pol * TOT_LEN}).to_csv("output/max_pol_rna.csv", index=False)
# pd.DataFrame({"gene" : })

#%% Sanity checks 
# double check that the names line up
prot_names = np.array(prot_ccd_ensg)[both_prot_ccd_idx]
rna_names =  np.array(rna_ccd_ensg)[both_rna_ccd_idx]
print(f"The name arrays are the same: {all(prot_names == rna_names)}")

# What are the smallest, largest, and median genes and what do they look like?
#smallest = np.argmin(diff_max_pol)
#smallest_gene = name_prot_list[smallest]
#median = np.argsort(diff_max_pol)[len(diff_max_pol)//2]
#median_gene = name_prot_list[median]
#largest = np.argmax(diff_max_pol)
#largest_gene = name_prot_list[largest]
#print(f"smallest delay {diff_max_pol[smallest] * TOT_LEN} hr for {name_prot_list[smallest]}")
#print(f"median delay {diff_max_pol[median] * TOT_LEN} hr for {name_prot_list[median]}")
#print(f"largest delay {diff_max_pol[largest] * TOT_LEN} hr for {name_prot_list[largest]}")
#
#plt.rcParams['figure.figsize'] = (10, 10)
#sorted_gene_array_array = np.array(sorted_gene_array).transpose()
#def plot_avg_rna_and_prot(namelist, outfolder):
#    if not os.path.exists(outfolder): os.mkdir(outfolder)
#    for name in namelist:
#        bin_size = 100
#        mvgavg = moving_averages[:,name_rna_list.index(name)]
#        mvgmax = mvgavg.max()
#        plt.plot(
#            fucci_time_sort[:-(bin_size-1)] * TOT_LEN, 
#            mvgavg / mvgmax,  
#            color="blue", 
#            label=f"RNA Moving Average by {bin_size} Cells")
#        prot_exp = sorted_gene_array_array[:,prot_list.index(name)]
#        plt.plot(
#            xvals[:-1] * TOT_LEN, 
#            prot_exp, 
#            color="red", 
#            label="Protein Expression")
#        plt.xlabel("Cell Division Time, hrs",size=36,fontname='Arial')
#        plt.ylabel("Expression, Normalized by Cell",size=36,fontname='Arial')
#        plt.xticks(size=14)
#        plt.yticks(size=14)
#        plt.title(name,size=36,fontname='Arial')
#        plt.legend(fontsize=14)
#        plt.tight_layout()
#        plt.savefig(f"{outfolder}/{name}.png")
#        plt.close()

#plot_avg_rna_and_prot(name_prot_list, "figures/RNAProteinCCDAvgs")

#%% Figures of merit
peaked_after_g1_prot = sorted_maxpol_array * TOT_LEN > G1_LEN
wp_ensg_counts_ccd = np.array([sum([eeee == ensg for eeee in wp_ensg[ccd_comp]]) for ensg in wp_ensg[ccd_comp]])
duplicated_ensg_ccd = wp_ensg_counts_ccd > 1
duplicated_ensg_peaked_after_g1 = np.array([sum(peaked_after_g1_prot[wp_ensg[ccd_comp] == ensg]) for ensg in duplicated_ensg_ccd])

with open("output/figuresofmerit.txt", "a") as file:
    fom = "--- temporal delay\n\n"
    fom += f"significant delay in peak protein expression compared to transcript expression, {TOT_LEN * np.median(diff_max_pol)} hours on average" + "\n\n"
    fom += f"G1 is the longest period of the cell cycle, in which the majority of RNAs ({100 * sum(sorted_max_moving_avg_pol_ccd * TOT_LEN <= G1_LEN) / len(sorted_max_moving_avg_pol_ccd)}%) peak in expression" + "\n\n"
    fom += f"However, the majority ({100 * (sum(peaked_after_g1_prot[~duplicated_ensg_ccd]) + sum(duplicated_ensg_peaked_after_g1 == 2)) / len(np.unique(wp_ensg[ccd_comp]))}%) of the proteins peaked towards the end of the cell cycle corresponding to the S&G2 phases" + "\n\n"
    fom += f"The delay between peak RNA and protein expression for the 50 CCD proteins that also had CCD transcripts was {TOT_LEN * np.median(diff_max_pol)} hrs on average " + "\n\n"
    fom += f"this delay indicates that it may take a little less than the same amount of time ({12 - TOT_LEN * np.median(diff_max_pol)} hrs) to produce a target metabolite after peak expression of an enzyme." + "\n\n"
    fom += f"" + "\n\n"
    fom += f"" + "\n\n"
    fom += f"" + "\n\n"
    print(fom)
    file.write(fom)