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
adata, phasesfilt = qc_filtering(adata, do_log_normalize=False, do_remove_blob=True)
ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)


#%% Read in the protein data
# Idea: read in the protein data to compare with the RNA seq data
# Exec: pandas
# Output: fucci plot from the immunofluorescence data
print("reading protein IF data")
plate = np.load("output/plate.filterNegStain.npy", allow_pickle=True)
u_plate = np.load("output/u_plate.filterNegStain.npy", allow_pickle=True)
well_plate=np.load("output/well_plate.filterNegStain.npy", allow_pickle=True)
imgnb=np.load("output/imgnb.filterNegStain.npy", allow_pickle=True)
u_well_plates=np.load("output/u_well_plates.filterNegStain.npy", allow_pickle=True)
ab_objnum=np.load("output/ab_objnum.filterNegStain.npy", allow_pickle=True)
ab_nuc=np.load("output/ab_nuc.filterNegStain.npy", allow_pickle=True)
ab_cyto=np.load("output/ab_cyto.filterNegStain.npy", allow_pickle=True)
ab_cell=np.load("output/ab_cell.filterNegStain.npy", allow_pickle=True)
mt_cell=np.load("output/mt_cell.filterNegStain.npy", allow_pickle=True)
green_fucci=np.load("output/green_fucci.filterNegStain.npy", allow_pickle=True)
red_fucci=np.load("output/red_fucci.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc=log_green_fucci_zeroc=np.load("output/log_green_fucci_zeroc.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc=np.load("output/log_red_fucci_zeroc.filterNegStain.npy", allow_pickle=True)
log_green_fucci_zeroc_rescale=np.load("output/log_green_fucci_zeroc_rescale.filterNegStain.npy", allow_pickle=True)
log_red_fucci_zeroc_rescale=np.load("output/log_red_fucci_zeroc_rescale.filterNegStain.npy", allow_pickle=True)
wp_cell_kruskal_adj = np.load("output/wp_cell_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_nuc_kruskal_adj = np.load("output/wp_nuc_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_cyto_kruskal_adj = np.load("output/wp_cyto_kruskal_adj.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_cell = np.load("output/wp_pass_bh_cell.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_nuc = np.load("output/wp_pass_bh_nuc.filterNegStain.npy", allow_pickle=True)
wp_pass_bh_cyto = np.load("output/wp_pass_bh_cyto.filterNegStain.npy", allow_pickle=True)
fucci_data = np.load("output/fucci_data.filterNegStain.npy", allow_pickle=True)
pol_sort_well_plate = np.load("output/pol_sort_well_plate.npy", allow_pickle= True)
pol_sort_norm_rev = np.load("output/pol_sort_norm_rev.npy", allow_pickle=True)
pol_sort_ab_nuc = np.load("output/pol_sort_ab_nuc.npy", allow_pickle=True)
pol_sort_ab_cyto = np.load("output/pol_sort_ab_cyto.npy", allow_pickle=True)
pol_sort_ab_cell = np.load("output/pol_sort_ab_cell.npy", allow_pickle=True)
pol_sort_mt_cell = np.load("output/pol_sort_mt_cell.npy", allow_pickle=True)
pol_sort_fred = np.load("output/pol_sort_fred.npy", allow_pickle=True)
pol_sort_fgreen = np.load("output/pol_sort_fgreen.npy", allow_pickle=True)
wp_ensg = np.load("output/wp_ensg.npy", allow_pickle=True)
wp_iscell = np.load("output/wp_iscell.npy", allow_pickle=True)
wp_isnuc = np.load("output/wp_isnuc.npy", allow_pickle=True)
wp_iscyto = np.load("output/wp_iscyto.npy", allow_pickle=True)
ccd_comp = np.load("output/ccd_comp.npy", allow_pickle=True)
print("loaded")

print("reading RNA data")
ccd_transcript_regulated = np.load("output/ccd_transcript_regulated.npy", allow_pickle=True)
dianaccd_transcript_regulated = np.load("output/dianaccd_transcript_regulated.npy", allow_pickle=True)
dianaccd_nontranscript_regulated = np.load("output/dianaccd_nontranscript_regulated.npy", allow_pickle=True)
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

xvals = np.linspace(0,1,num=21)
wp_max_pol = []
wp_binned_values = []
for i, well in enumerate(u_well_plates):
    curr_well_inds = pol_sort_well_plate==well
    curr_pol = pol_sort_norm_rev[curr_well_inds]
    curr_fred = pol_sort_fred[curr_well_inds]
    curr_fgreen = pol_sort_fgreen[curr_well_inds]
    curr_ab_cell, curr_ab_nuc, curr_ab_cyto, curr_mt_cell = pol_sort_ab_cell[curr_well_inds], pol_sort_ab_nuc[curr_well_inds],pol_sort_ab_cyto[curr_well_inds], pol_sort_mt_cell[curr_well_inds]

    # Normalize FUCCI colors & mean intensities, normalized for display
    curr_fred_norm = curr_fred / np.max(curr_fred)
    curr_fgreen_norm = curr_fgreen / np.max(curr_fgreen)
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

# Make an expression array with the CCD proteins
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
plt.savefig(os.path.join("output_devin",'sorted_heatmap21_sw30_take4.pdf'), transparent=True)
plt.savefig(os.path.join("output_devin",'sorted_heatmap21_sw30_take4.png'), transparent=True)
plt.show()

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
max_moving_avg_pol_ccd = max_moving_avg_pol[ccd_transcript_regulated]
moving_averages_ccd = moving_averages[:,ccd_transcript_regulated]
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
rna_ccd_ensg = list(adata.var_names[ccd_transcript_regulated])
both_ccd_ensg = np.intersect1d(prot_ccd_ensg, rna_ccd_ensg)
both_prot_ccd_idx = np.array([prot_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
both_rna_ccd_idx = np.array([rna_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
insct_prot_max_pol_ccd = wp_max_pol_ccd[both_prot_ccd_idx]
insct_rna_max_pol_ccd = sorted_max_moving_avg_pol_ccd[both_rna_ccd_idx]
diff_max_pol = insct_prot_max_pol_ccd - insct_rna_max_pol_ccd

print(f"length of prot CCD genes: {len(prot_ccd_ensg)}")
print(f"length of CCD RNA genes: {len(rna_ccd_ensg)}")
print(f"length of intersection betweeen CCD prot and CCD RNA: {len(both_ccd_ensg)}")

plt.hist(diff_max_pol * TOT_LEN)
plt.xlabel("Delay in peak protein expression from peak RNA expression, hrs")
plt.ylabel("Count of CCD Proteins")
plt.tight_layout()
plt.savefig("figures/DelayPeakProteinRNA.png")
plt.show()
plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hist(insct_prot_max_pol_ccd * TOT_LEN, alpha=0.5, label="Peak Protein Expression Time, hrs")
ax1.hist(insct_rna_max_pol_ccd * TOT_LEN, alpha=0.5, label="Peak RNA Expression Time, hrs")
plt.legend(loc="upper left")
plt.xlabel("Division Cycle, hrs")
plt.ylabel("Count of Cell Cycle Genes")
plt.tight_layout()
plt.savefig(f"figures/DelayPeakProteinRNA_separate.png")
plt.show()
plt.close()

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

