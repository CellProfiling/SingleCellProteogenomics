# -*- coding: utf-8 -*-
"""
Evaluates the delay between peak protein and RNA expression over the cell cycle:
    - Creates a heatmap of time of peak expression for protein and for RNA
    - Evaluates correlation of temporal expression for known and novel cell cycle dependent (CCD) proteins
    - Illustrates the transitions between phases of peak expression for protein and RNA of the each gene that is CCD for both

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils, FucciCellCycle, MovingAverages, RNADataPreparation, alluvial
from mpl_toolkits.axes_grid1 import make_axes_locatable
import SingleCellProteogenomics.alluvial
import itertools

fucci = FucciCellCycle.FucciCellCycle() # Object representing FUCCI cell cycle phase durations

def protein_heatmap(nbins, highlight_names, highlight_ensg, ccd_comp, u_well_plates, wp_ensg, pol_sort_norm_rev, pol_sort_well_plate, 
                    pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, wp_iscell, wp_isnuc, wp_iscyto):
    '''Make heatmap of the time of peak expression for CCD proteins. Highlight examples on the y-axis if given.'''
    # Make an expression array with the CCD proteins
    wp_max_pol, wp_binned_values, xvals = MovingAverages.bin_values(nbins, u_well_plates, pol_sort_norm_rev, pol_sort_well_plate, pol_sort_ab_cell, pol_sort_ab_nuc, pol_sort_ab_cyto, pol_sort_mt_cell, wp_iscell, wp_isnuc, wp_iscyto)
    wp_max_pol, wp_binned_values = np.array(wp_max_pol), np.array(wp_binned_values)
    wp_max_pol_ccd, wp_binned_values_ccd = wp_max_pol[ccd_comp], wp_binned_values[ccd_comp]
    wp_max_sort_inds = np.argsort(wp_max_pol_ccd)
    sorted_gene_array = np.take(wp_binned_values_ccd, wp_max_sort_inds, axis=0) # these are the expression values (binned_values), sorted by the binned value at max location (can do in the temporal part)
    sorted_maxpol_array = np.take(wp_max_pol_ccd, wp_max_sort_inds)

    # Actually making the figure
    fig, ax = plt.subplots(figsize=(10, 10))
    sc = ax.imshow(sorted_gene_array, interpolation='nearest')

    # Do the x ticks
    xtick_labels = [str(np.around(x * fucci.TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)] #+['G1/S','S/G2']
    my_xticks = np.arange(-.5, nbins, 2)
    num_ticks = nbins
    xphase_labels = ['G1/S','S/G2']
    phase_trans = np.asarray([fucci.G1_PROP*num_ticks-0.5, fucci.G1_S_PROP*num_ticks-0.5])
    ax.set_xticks(my_xticks,minor=True)
    ax.set_xticklabels(xtick_labels,minor=True)
    ax.set_xticks(phase_trans, minor=False)
    ax.set_xticklabels(xphase_labels, minor=False)
    ax.tick_params(length=12)

    # Do the y ticks
    if len(highlight_names) > 0:
        ytick_locs = [wp_ensg.index(ensg) for ensg in highlight_ensg]
        ax.set_yticks(ytick_locs, minor=False)
        ax.set_yticklabels(highlight_names,minor=False)

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

    return sorted_maxpol_array, wp_binned_values, wp_max_pol, wp_max_pol_ccd, xvals

def scatter_genes(gene1, gene2, r, wp_binned_values, wp_ensg):
    '''Make a scatterplot to go along with the correlation of protein expression between two genes'''
    if not os.path.exists("figures/Correlations"): os.mkdir("figures/Correlations")
    plt.scatter(wp_binned_values[wp_ensg == gene1[0]][0], wp_binned_values[wp_ensg == gene2[0]][0])
    plt.xlabel(f"{gene1[1]} Expression Binned by Pseudotime", fontsize=14, fontname="Arial")
    plt.ylabel(f"{gene2[1]} Expression Binned by Pseudotime", fontsize=14, fontname="Arial")
    plt.text(np.min(wp_binned_values[wp_ensg == gene1[0]][0]), np.max(wp_binned_values[wp_ensg == gene2[0]][0]), f"Pearson's r = {r}", fontsize=14, fontname="Arial")
    plt.savefig(f"figures/Correlations/{gene1[1]}_{gene2[1]}.pdf")
    plt.show()
    plt.close()

def peak_expression_correlation_analysis(wp_binned_values, wp_max_pol, wp_ensg, pol_sort_well_plate, u_well_plates):
    '''Perform correlation analysis between known and novel CCD proteins that peak at similar times'''
    prevfigsize = plt.rcParams['figure.figsize']
    plt.rcParams['figure.figsize'] = (5, 5)

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
             ("ENSG00000102908", "nfat5")]: # unknown
        print(ensg)
        print(f"number of observations: {sum(pol_sort_well_plate==u_well_plates[wp_ensg == ensg[0]][0])}")
        print(f"time of peak expression: {wp_max_pol[wp_ensg == ensg[0]][0]}")

    orc6_znf32 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000091651'][0], wp_binned_values[wp_ensg == 'ENSG00000169740'][0])
    scatter_genes(("ENSG00000091651", "ORC6"), ("ENSG00000169740", "ZNF32"), round(orc6_znf32[0], 2), wp_binned_values, wp_ensg)

    ccne1_dusp19 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000105173'][0], wp_binned_values[wp_ensg == 'ENSG00000162999'][0])
    scatter_genes(("ENSG00000105173", "CCNE1"), ("ENSG00000162999", "DUSP19"), round(ccne1_dusp19[0], 2), wp_binned_values, wp_ensg)

    ttc21b_pc = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000123607'][0], wp_binned_values[wp_ensg == 'ENSG00000173599'][0])
    scatter_genes(("ENSG00000123607", "TTC21B"), ("ENSG00000173599", "PC"), round(ttc21b_pc[0], 2), wp_binned_values, wp_ensg)

    bub1b_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000156970'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
    scatter_genes(("ENSG00000156970", "BUB1B"), ("ENSG00000167065", "DUSP18"), round(bub1b_dusp18[0], 2), wp_binned_values, wp_ensg)

    aurkb_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000178999'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
    scatter_genes(("ENSG00000178999", "AURKB"), ("ENSG00000167065", "DUSP18"), round(aurkb_dusp18[0], 2), wp_binned_values, wp_ensg)

    ccnb1_dusp18 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000167065'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000167065", "DUSP18"), round(ccnb1_dusp18[0], 2), wp_binned_values, wp_ensg)

    ccnb1_papss1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000138801'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000138801", "PAPSS1"), round(ccnb1_papss1[0], 2), wp_binned_values, wp_ensg)

    ccnb1_n6amt1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000156239'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000156239", "N6AMT1"), round(ccnb1_n6amt1[0], 2), wp_binned_values, wp_ensg)

    ccnb1_phldb1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000019144'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000019144", "PHLDB1"), round(ccnb1_phldb1[0], 2), wp_binned_values, wp_ensg)

    ccnb1_fli1 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000151702'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000151702", "FLI1"), round(ccnb1_fli1[0], 2), wp_binned_values, wp_ensg)

    ccnb1_dph2 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000132768'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000132768", "DPH2"), round(ccnb1_dph2[0], 2), wp_binned_values, wp_ensg)

    ccnb1_nfat5 = scipy.stats.pearsonr(wp_binned_values[wp_ensg == 'ENSG00000134057'][0], wp_binned_values[wp_ensg == 'ENSG00000102908'][0])
    scatter_genes(("ENSG00000134057", "CCNB1"), ("ENSG00000102908", "NFAT5"), round(ccnb1_nfat5[0], 2), wp_binned_values, wp_ensg)

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

    plt.rcParams['figure.figsize'] = prevfigsize

def binned_median(yvals, nbins):
    '''Compute RNA expression values, binned over pseudotime in `nbins` number of bins'''
    binned_medians = []
    for xval in range(nbins):
        startidx = len(yvals) // nbins * xval
        endidx = len(yvals) // nbins * (xval + 1)
        binned_medians.append(np.median(yvals[startidx:endidx]))
    return binned_medians

def rna_heatmap(adata, highlight_names, highlight_ensg, ccdtranscript, xvals, isIsoformData=False):
    '''Make heatmap of the time of peak expression for CCD transcripts. Highlight examples on the y-axis if given.'''
    # Get the peak RNA expression polar locations
    expression_data = adata.X # log normalized
    normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
    fucci_time_inds = np.argsort(adata.obs["fucci_time"])
    fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
    norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
    moving_averages = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort, 100)
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
    xtick_labels = [str(np.around(x * fucci.TOT_LEN,decimals=2)) for x in np.linspace(0,1,11)] #+['G1/S','S/G2']
    my_xticks = np.arange(-.5, 20, 2)
    num_ticks = 20
    xphase_labels = ['G1/S','S/G2']
    phase_trans = np.asarray([fucci.G1_PROP*num_ticks-0.5, fucci.G1_S_PROP*num_ticks-0.5])
    ax.set_xticks(my_xticks,minor=True)
    ax.set_xticklabels(xtick_labels,minor=True)
    ax.set_xticks(phase_trans, minor=False)
    ax.set_xticklabels(xphase_labels, minor=False)
    ax.tick_params(length=12)

    #Do the y ticks
    if len(highlight_names) > 0:
        ytick_locs = [list(adata.var_names).index(ensg) for ensg in highlight_ensg]
        ax.set_yticks(ytick_locs, minor=False)
        ax.set_yticklabels(highlight_names,minor=False)

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
    plt.savefig(os.path.join("figures",f"sorted_rna_heatmap{'_isoform' if isIsoformData else ''}.pdf"), transparent=True)
    plt.savefig(os.path.join("figures",f"sorted_rna_heatmap{'_isoform' if isIsoformData else ''}.png"), transparent=True)
    plt.show()

    return sorted_max_moving_avg_pol_ccd, norm_exp_sort, max_moving_avg_pol, sorted_rna_binned_norm

def compare_variances_prot_v_rna(adata, norm_exp_sort, wp_ensg, var_comp_prot, gini_comp_prot, cv_comp_prot, var_cell_prot, gini_cell_prot, cv_cell_prot):
    '''Compare the measures of variance for RNA and protein and make scatterplots'''
    total_variance_rna = np.var(norm_exp_sort, 0)
    total_gini_rna = np.apply_along_axis(utils.gini, 0, norm_exp_sort)
    total_cv_rna = np.apply_along_axis(scipy.stats.variation, 0, norm_exp_sort)

    prot_ensg = list(wp_ensg)
    rna_ensg = list(adata.var_names)
    both_ensg = np.intersect1d(prot_ensg, rna_ensg)
    both_prot_idx = np.array([prot_ensg.index(ensg) for ensg in both_ensg])
    both_rna_idx = np.array([rna_ensg.index(ensg) for ensg in both_ensg])
    insct_prot_variance_comp = var_comp_prot[both_prot_idx]
    insct_prot_gini_comp = gini_comp_prot[both_prot_idx]
    insct_prot_cv_comp = cv_comp_prot[both_prot_idx]
    insct_prot_variance_cell = var_cell_prot[both_prot_idx]
    insct_prot_gini_cell = gini_cell_prot[both_prot_idx]
    insct_prot_cv_cell = cv_cell_prot[both_prot_idx]
    insct_rna_variance = total_variance_rna[both_rna_idx]
    insct_rna_gini = total_gini_rna[both_rna_idx]
    insct_rna_cv = total_cv_rna[both_rna_idx]

    utils.general_scatter(insct_rna_variance, insct_prot_variance_cell, "Total RNA Variance", "Total Protein Variance (Cell)", "figures/ProteinRNAVariance.png")
    utils.general_scatter(insct_rna_gini, insct_prot_gini_cell, "RNA Gini", "Protein Gini  (Cell)", "figures/ProteinRNAGini.png")
    utils.general_scatter(insct_rna_cv, insct_prot_cv_cell, "RNA CV", "Protein CV  (Cell)", "figures/ProteinRNACV.png")
    
    utils.general_scatter(insct_rna_variance, insct_prot_variance_comp, "Total RNA Variance", "Total Protein Variance (Compartment)", "figures/ProteinRNAVariance.png")
    utils.general_scatter(insct_rna_gini, insct_prot_gini_comp, "RNA Gini", "Protein Gini  (Compartment)", "figures/ProteinRNAGini.png")
    utils.general_scatter(insct_rna_cv, insct_prot_cv_comp, "RNA CV", "Protein CV  (Compartment)", "figures/ProteinRNACV.png")

    pd.DataFrame({
        "gene":both_ensg,
        "variance_rna": insct_rna_variance,
        "gini_rna":insct_rna_gini,
        "cv_rna":insct_rna_cv,
        "variance_comp_prot":insct_prot_variance_comp,
        "gini_comp_prot":insct_prot_gini_comp,
        "cv_comp_prot":insct_prot_cv_comp,
        "variance_cell_prot":insct_prot_variance_cell,
        "gini_cell_prot":insct_prot_gini_cell,
        "cv_cell_prot":insct_prot_cv_cell,
        }).to_csv("output/VarianceRNAProtein.csv",index=False)
    
def peak_expression_alluvial(diff_max_pol, insct_rna_max_pol_ccd, insct_prot_max_pol_ccd):
    '''Generate alluvial plot comparing the phase of peak expression for each gene'''
    transitiondict = {}
    startphases = ["G1\nRNA\nPeak", "S\nRNA\nPeak", "G2\nRNA\nPeak"]
    endphases= ["G1\nProtein\nPeak", "S\nProtein\nPeak", "G2\nProtein\nPeak"]
    for iii, diff in enumerate(list(diff_max_pol)):
        rnatime = insct_rna_max_pol_ccd[iii]
        prottime = insct_prot_max_pol_ccd[iii]
        startloc = startphases[0] if rnatime * fucci.TOT_LEN < fucci.G1_LEN else startphases[1] if rnatime * fucci.TOT_LEN >= fucci.G1_LEN and rnatime * fucci.TOT_LEN < fucci.G1_LEN + fucci.G1_S_TRANS else startphases[2]
        endloc = endphases[0] if prottime * fucci.TOT_LEN < fucci.G1_LEN else endphases[1] if prottime * fucci.TOT_LEN >= fucci.G1_LEN and prottime * fucci.TOT_LEN < fucci.G1_LEN + fucci.G1_S_TRANS else endphases[2]
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

def peak_expression_delay_scatter(insct_rna_max_pol_ccd, insct_prot_max_pol_ccd, diff_max_pol):
    '''Make a scatterplot with time of peak RNA and protein expression colored by the temporal delay'''
    f,ax = plt.subplots(figsize=(6,5))
    plt.scatter(x=insct_rna_max_pol_ccd * fucci.TOT_LEN, y=insct_prot_max_pol_ccd * fucci.TOT_LEN, c=diff_max_pol * fucci.TOT_LEN)
    # ax.hist(diff_max_pol * TOT_LEN, orientation="vertical", alpha=0.5)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_xticks([0,5,10,15,20])
    cbar = plt.colorbar(pad=0.15)
    cbar.set_label('Temporal Delay, hrs',fontname='Arial',size=20)
    cbar.ax.tick_params(labelsize=18)
    plt.xlabel("Peak RNA Expression, hrs",fontname='Arial',size=20)
    plt.ylabel("Peak Protein Expression, hrs",fontname='Arial',size=20)
    plt.tight_layout()
    plt.savefig("figures/TemporalDelayScatter.pdf")
    plt.show()
    plt.close()

def compare_peak_expression_prot_v_rna(adata, wp_ensg, ccd_comp, ccdtranscript, wp_max_pol, wp_max_pol_ccd, sorted_maxpol_array, max_moving_avg_pol, sorted_max_moving_avg_pol_ccd):
    '''Compare the time of peak expression of protein and RNA'''
    prot_ccd_ensg = list(wp_ensg[ccd_comp])
    rna_ccd_ensg = list(adata.var_names[ccdtranscript])
    both_ccd_ensg = np.intersect1d(prot_ccd_ensg, rna_ccd_ensg)
    both_prot_ccd_idx = np.array([prot_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
    both_rna_ccd_idx = np.array([rna_ccd_ensg.index(ensg) for ensg in both_ccd_ensg])
    insct_prot_max_pol_ccd = wp_max_pol_ccd[both_prot_ccd_idx]
    insct_rna_max_pol_ccd = sorted_max_moving_avg_pol_ccd[both_rna_ccd_idx]
    diff_max_pol = insct_prot_max_pol_ccd - insct_rna_max_pol_ccd

    #% Sanity check: double check that the names line up
    prot_names = np.array(prot_ccd_ensg)[both_prot_ccd_idx]
    rna_names =  np.array(rna_ccd_ensg)[both_rna_ccd_idx]
    print(f"The name arrays are the same: {all(prot_names == rna_names)}")

    # Alluvial plot showing RNA-protein phase of peak expression for each genes
    peak_expression_alluvial(diff_max_pol, insct_rna_max_pol_ccd, insct_prot_max_pol_ccd)
    # Histogram for peak expression
    utils.general_histogram(diff_max_pol * fucci.TOT_LEN, "Delay in peak protein expression from peak RNA expression, hrs", "Count of CCD Proteins", 0.5, "figures/DelayPeakProteinRNA.pdf")
    # Scatter for peak expression with colorbar for the delay
    peak_expression_delay_scatter(insct_rna_max_pol_ccd, insct_prot_max_pol_ccd, diff_max_pol)
    # Boxplot for delay of peak expression
    utils.general_boxplot((insct_prot_max_pol_ccd * fucci.TOT_LEN, insct_rna_max_pol_ccd * fucci.TOT_LEN), ("Protein", "RNA"), 
        "", "Peak Expression, hrs", "", True, "figures/DelayPeakProteinRNA_boxplot.png")

    print(f"Count of prot CCD genes: {len(prot_ccd_ensg)}")
    print(f"Count of CCD RNA genes: {len(rna_ccd_ensg)}")
    print(f"Count of intersection betweeen CCD prot and CCD RNA: {len(both_ccd_ensg)}")
    print(f"Median delay of RNA and protein expression time for CCD proteins: {fucci.TOT_LEN * np.median(diff_max_pol)}")
    print(f"Median RNA expression time for CCD proteins: {fucci.TOT_LEN * np.median(insct_rna_max_pol_ccd)}")
    print(f"Median protein expression time for CCD proteins: {fucci.TOT_LEN * np.median(insct_prot_max_pol_ccd)}")
    t, p = scipy.stats.kruskal(insct_rna_max_pol_ccd, insct_prot_max_pol_ccd)
    print(f"One-sided kruskal for median protein expression time higher than median RNA expression time: {2*p}")
    t, p = scipy.stats.ttest_1samp(diff_max_pol, 0)
    print(f"One-sided, one-sample t-test for mean delay in protein expression larger than zero: {2*p}")

    #% Output tables
    pd.DataFrame({"gene" : wp_ensg, "max_pol_protein": wp_max_pol, "max_time_protein": wp_max_pol * fucci.TOT_LEN}).to_csv("output/max_pol_protein.csv", index=False)
    pd.DataFrame({"gene" : adata.var_names, "max_pol_rna": max_moving_avg_pol, "max_time_rna": max_moving_avg_pol * fucci.TOT_LEN}).to_csv("output/max_pol_rna.csv", index=False)

    #% Figures of merit
    peaked_after_g1_prot = sorted_maxpol_array * fucci.TOT_LEN > fucci.G1_LEN
    wp_ensg_counts_ccd = np.array([sum([eeee == ensg for eeee in wp_ensg[ccd_comp]]) for ensg in wp_ensg[ccd_comp]])
    duplicated_ensg_ccd = wp_ensg_counts_ccd > 1
    duplicated_ensg_peaked_after_g1 = np.array([sum(peaked_after_g1_prot[wp_ensg[ccd_comp] == ensg]) for ensg in duplicated_ensg_ccd])
    with open("output/figuresofmerit.txt", "a") as file:
        fom = "--- temporal delay\n\n"
        fom += f"significant delay in peak protein expression compared to transcript expression, {fucci.TOT_LEN * np.median(diff_max_pol)} hours on average" + "\n\n"
        fom += f"G1 is the longest period of the cell cycle, in which the majority of RNAs ({100 * sum(sorted_max_moving_avg_pol_ccd * fucci.TOT_LEN <=fucci. G1_LEN) / len(sorted_max_moving_avg_pol_ccd)}%) peak in expression" + "\n\n"
        fom += f"However, the majority ({100 * (sum(peaked_after_g1_prot[~duplicated_ensg_ccd]) + sum(duplicated_ensg_peaked_after_g1 == 2)) / len(np.unique(wp_ensg[ccd_comp]))}%) of the proteins peaked towards the end of the cell cycle corresponding to the S&G2 phases" + "\n\n"
        fom += f"The delay between peak RNA and protein expression for the 50 CCD proteins that also had CCD transcripts was {fucci.TOT_LEN * np.median(diff_max_pol)} hrs on average " + "\n\n"
        fom += f"this delay indicates that it may take a little less than the same amount of time ({12 - fucci.TOT_LEN * np.median(diff_max_pol)} hrs) to produce a target metabolite after peak expression of an enzyme." + "\n\n"
        fom += f"" + "\n\n"
        fom += f"" + "\n\n"
        fom += f"" + "\n\n"
        print(fom)
        file.write(fom)
        
def analyze_ccd_isoform_correlations(adata, adata_isoform, ccdtranscript, ccdtranscript_isoform, xvals):
    '''Evaluate the pearson correlations for CCD isoforms from genes with multiple CCD isoforms'''
    gene_varnames, isoform_varnames = list(adata.var_names), list(adata_isoform.var_names)
    isoformToGene = pd.read_csv("input/processed/python/IsoformToGene.csv", index_col=False, header=None, names=["transcript_id", "gene_id"])
    isoformIdList = list(isoformToGene["transcript_id"])
    isoform_varnames_geneids = np.array([isoformToGene["gene_id"][isoformIdList.index(t)] for t in isoform_varnames])
    ccdIsoformWithCcdGene = ccdtranscript_isoform[np.isin(isoform_varnames_geneids, gene_varnames)] & np.array([ccdtranscript[gene_varnames.index(gene_id)] for gene_id in isoform_varnames_geneids if gene_id in gene_varnames])
    numIsoformsPerGene = isoformToGene.groupby("gene_id")["gene_id"].value_counts()
    perGene_geneIds = np.array([x[0] for x in numIsoformsPerGene.index])
    useGene = np.isin(perGene_geneIds, gene_varnames)
    numIsoformsPerGene = np.array(numIsoformsPerGene[useGene])
    ccdIsoformsPerGene = np.array([sum(ccdtranscript_isoform[isoform_varnames_geneids == gene_id]) for gene_id in perGene_geneIds[useGene]])
    ccdAndNonCcdIsoformsPerGene = np.array([numIsoformsPerGene[ii] != ccdIsoformsPerGene[ii] for ii, gene_id in enumerate(numIsoformsPerGene)])
    
    isoformsFromGenesWithMultipleCCD = [adata_isoform.var_names[(isoform_varnames_geneids == gene_id) & ccdtranscript_isoform] for gene_id in perGene_geneIds[useGene][ccdIsoformsPerGene > 1]]
    # isoformsFromGenesWithMultipleCCD_maxpol = [max_moving_avg_pol_isoform[(isoform_varnames_geneids == gene_id) & ccdtranscript_isoform] for gene_id in perGene_geneIds[useGene][ccdIsoformsPerGene > 1]]
   
    expression_data = adata_isoform.X # log normalized
    normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
    fucci_time_inds = np.argsort(adata_isoform.obs["fucci_time"])
    fucci_time_sort = np.take(np.array(adata_isoform.obs["fucci_time"]), fucci_time_inds)
    norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)
    moving_averages = np.apply_along_axis(MovingAverages.mvavg, 0, norm_exp_sort, 100)
    binned=np.apply_along_axis(binned_median, 1, moving_averages.T, len(xvals))
    binned_norm = binned / np.max(binned, axis=1)[:,None]
    
    isoformsFromGenesWithMultipleCCD_binned = [binned_norm[(isoform_varnames_geneids == gene_id) & ccdtranscript_isoform] for gene_id in perGene_geneIds[useGene][ccdIsoformsPerGene > 1]]
    pearsonCorrelations = [[scipy.stats.pearsonr(combo[0], combo[1])[0] for combo in itertools.combinations(a, 2)] for a in isoformsFromGenesWithMultipleCCD_binned]
    adata_isoform_raw, phases_isoform_raw = RNADataPreparation.read_counts_and_phases("Tpms", False, "protein_coding", use_isoforms=True)
    print("ISOFORM PEARSON CORRELATIONS FOR CCD GENES WITH MULTIPLE CCD ISOFORMS")
    for ii, arr in enumerate(pearsonCorrelations):
        if sum(np.array(arr) < 0.5) == 0: continue
        print(f"{','.join(isoformsFromGenesWithMultipleCCD[ii])}\tPearson R's:{arr}\tTPM values:{' '.join([str(np.mean(adata_isoform_raw.X[:, list(adata_isoform_raw.var_names).index(t)])) for t in isoformsFromGenesWithMultipleCCD[ii]])}")
    return pearsonCorrelations