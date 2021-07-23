#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 16:51:51 2021

Methods for analyzing the protein cell cycle depence results

@author: anthony.cesnik
"""

from SingleCellProteogenomics import (
    utils,
    MovingAverages,
)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
import scipy
import warnings, pickle, shutil, os
from sklearn.linear_model import MultiTaskLassoCV

# Number of bins for creating UMAPs/LASSO model. Chosen for the best stability.
BINS_FOR_UMAP_AND_LASSO = 400

# Used for getting median FUCCI marker intensity for LASSO analysis
WINDOW_FUCCI_MARKERS = 100

# Cutoffs chosen for umaps
chosen_nn = 5
chosen_md = 0.2


def global_plots_protein(protein_ccd, alpha_ccd):
    """Illustrate the CCD variances of all proteins"""
    if not protein_ccd.do_plotting:
        return

    utils.general_scatter(
        protein_ccd.perc_var_comp,
        protein_ccd.protein_variability.mean_mean_comp,
        "percent variance",
        "mean mean intensity",
        "figures/PercVarVsMeanMeanIntensity_comp.png",
    )
    utils.general_scatter(
        protein_ccd.protein_variability.mean_mean_comp,
        protein_ccd.mean_diff_from_rng,
        "Mean Mean Intensity",
        "Mean Additional Percent Variance Explained than Random",
        "figures/IntensityVsMeanDiff.png",
    )
    utils.general_scatter_color(
        protein_ccd.protein_variability.gini_comp,
        protein_ccd.perc_var_comp,
        "Gini of Protein Expression",
        "Fraction of Variance Due to Cell Cycle",
        -np.log10(protein_ccd.protein_bimodality.wp_comp_kruskal_gaussccd_adj),
        "FDR for Cell Cycle Dependence",
        True,
        "Compartment - Fraction of Variance Due to Cell Cycle",
        "figures/CompartmentProteinFractionVariance.png",
    )
    utils.general_scatter_color(
        protein_ccd.protein_variability.gini_comp,
        protein_ccd.perc_var_comp,
        "Gini of Protein Expression",
        "Fraction of Variance Due to Cell Cycle",
        protein_ccd.wp_ccd_unibimodal,
        "",
        False,
        "Compartment - Fraction of Variance Due to Cell Cycle",
        "figures/CompartmentProteinFractionVarianceTF.png",
        "bwr_r",
        0.5,
    )
    utils.general_scatter_color(
        protein_ccd.protein_variability.cv_comp,
        protein_ccd.perc_var_comp,
        "CV of Protein Expression",
        "Fraction of Variance Due to Cell Cycle",
        -np.log10(protein_ccd.protein_bimodality.wp_comp_kruskal_gaussccd_adj),
        "FDR for Cell Cycle Dependence",
        True,
        "Compartment - Fraction of Variance Due to Cell Cycle",
        "figures/CompartmentCVProteinFractionVariance.png",
    )

    pervar_eq_percvar_adj = np.nextafter(
        protein_ccd.wp_comp_eq_percvar_adj, protein_ccd.wp_comp_eq_percvar_adj + 1
    )
    plt.figure(figsize=(10, 10))
    plt.scatter(protein_ccd.perc_var_comp, -np.log10(protein_ccd.pervar_eq_percvar_adj))
    plt.xlabel("percent variance new")
    plt.ylabel("-log10 FDR for CCD")
    plt.hlines(
        -np.log10(alpha_ccd),
        np.min(protein_ccd.perc_var_comp),
        np.max(protein_ccd.perc_var_comp),
    )
    plt.savefig("figures/PercVarVsLog10FdrCCD_comp.png")
    plt.close()


def compare_to_lasso_analysis(protein_ccd):
    """Comparison of pseudotime alignment to LASSO for finding CCD proteins"""
    proteins = pd.read_csv("output/ProteinPseudotimePlotting.csv.gz", sep="\t")
    numCells = np.array([len(x.split(",")) for x in proteins["cell_fred"]])
    do_normalize = False
    wp_max_pol, wp_binned_values, xvals = MovingAverages.bin_values(
        BINS_FOR_UMAP_AND_LASSO,
        protein_ccd.u_well_plates,
        protein_ccd.fucci_polar_coords.pol_sort_norm_rev,
        protein_ccd.fucci_polar_coords.pol_sort_well_plate,
        protein_ccd.fucci_polar_coords.pol_sort_ab_cell,
        protein_ccd.fucci_polar_coords.pol_sort_ab_nuc,
        protein_ccd.fucci_polar_coords.pol_sort_ab_cyto,
        protein_ccd.fucci_polar_coords.pol_sort_mt_cell,
        protein_ccd.wp_iscell,
        protein_ccd.wp_isnuc,
        protein_ccd.wp_iscyto,
        do_normalize=do_normalize,
    )
    wp_binned_values = np.array(wp_binned_values)

    # Bin the FUCCI coordinates
    mvavg_red = MovingAverages.mvavg(
        protein_ccd.fucci_polar_coords.pol_sort_fred, WINDOW_FUCCI_MARKERS
    )
    mvavg_green = MovingAverages.mvavg(
        protein_ccd.fucci_polar_coords.pol_sort_fgreen, WINDOW_FUCCI_MARKERS
    )
    mvavg_xvals = MovingAverages.mvavg(
        protein_ccd.fucci_polar_coords.pol_sort_norm_rev, WINDOW_FUCCI_MARKERS
    )

    # plot fucci coordinate averages
    plt.scatter(mvavg_red, mvavg_green, c=mvavg_xvals)
    plt.xlabel("Centered CDT1 log10 Intensity")
    plt.ylabel("Centered GMNN log10 Intensity")
    cb = plt.colorbar()
    cb.set_label("Pseudotime")
    plt.savefig("figures/FucciCoordinateAverages.png")
    plt.close()

    xvals = np.linspace(0, 1, num=BINS_FOR_UMAP_AND_LASSO)
    wp_max_pol = []
    binned_values_fred, binned_values_fgreen = [], []
    for xval in xvals:
        if xval == 0:
            prev_xval = xval
            continue
        binned_values_fred.append(
            np.median(mvavg_red[(mvavg_xvals < xval) & (mvavg_xvals >= prev_xval)])
        )
        binned_values_fgreen.append(
            np.median(mvavg_green[(mvavg_xvals < xval) & (mvavg_xvals >= prev_xval)])
        )
        prev_xval = xval
    binned_values_fred, binned_values_fgreen = np.array(binned_values_fred), np.array(
        binned_values_fgreen
    )

    # plot binned fucci coordinate averages
    plt.scatter(binned_values_fred, binned_values_fgreen, c=xvals[1:])
    plt.xlabel("Binned Centered CDT1 log10 Intensity")
    plt.ylabel("Binned Centered GMNN log10 Intensity")
    cb = plt.colorbar()
    cb.set_label("Pseudotime")
    plt.savefig("figures/FucciCoordinateBinnedAverages.png")
    plt.close()

    protein_fucci = np.vstack((binned_values_fred, binned_values_fgreen))
    fucci_protein_path = f"output/pickles/fucci_protein_lasso_binned{BINS_FOR_UMAP_AND_LASSO}{'Norm' if do_normalize else 'NoNorm'}.pkl"
    if os.path.exists(fucci_protein_path):
        fucci_protein = np.load(open(fucci_protein_path, "rb"), allow_pickle=True)
    else:
        fucci_protein = MultiTaskLassoCV()
        fucci_protein.fit(np.array(wp_binned_values).T, protein_fucci.T)
        pickle.dump(fucci_protein, open(fucci_protein_path, "wb"))
    plt.scatter(fucci_protein.alphas_, np.mean(fucci_protein.mse_path_, axis=1))
    plt.xlim((np.min(fucci_protein.alphas_), np.max(fucci_protein.alphas_)))
    lasso_statements = [
        f"{sum(np.sum(fucci_protein.coef_, axis=0) != 0)}: number of nonzero lasso coefficients",
        f"{protein_ccd.wp_ensg[np.sum(fucci_protein.coef_, axis=0) != 0]}: genes with nonzero lasso coeff",
        f"{sum(protein_ccd.ccd_comp[np.sum(fucci_protein.coef_, axis=0) != 0])}: CCD protein with nonzero lasso coeff",
        f"{np.sum(fucci_protein.coef_, axis=0)[np.sum(fucci_protein.coef_, axis=0) != 0]}",
    ]
    print("\n".join(lasso_statements))

    # Make a UMAPs for the LASSO analysis to demonstrate higher false negative rate
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nz_coef_protein = np.sum(fucci_protein.coef_, axis=0) != 0
        reducer = umap.UMAP(n_neighbors=chosen_nn, min_dist=chosen_md, random_state=0)
        embeddingCcd = reducer.fit_transform(wp_binned_values[nz_coef_protein, :].T)
        plt.scatter(embeddingCcd[:, 0], embeddingCcd[:, 1], c=xvals[1:])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        cb = plt.colorbar()
        cb.set_label("Pseudotime")
        plt.savefig("figures/umapProteinLassoCCD.pdf")
        plt.close()

        embeddingNonCcd = reducer.fit_transform(wp_binned_values[~nz_coef_protein, :].T)
        plt.scatter(embeddingNonCcd[:, 0], embeddingNonCcd[:, 1], c=xvals[1:])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        cb = plt.colorbar()
        cb.set_label("Pseudotime")
        plt.savefig("figures/umapProteinLassoNonCCD.pdf")
        plt.close()


def generate_protein_umaps(protein_ccd):
    if not os.path.exists("figures/ProteinUmaps"):
        os.mkdir("figures/ProteinUmaps")
    if not os.path.exists("figures/ProteinUmapStability"):
        os.mkdir("figures/ProteinUmapStability")

    warnings.filterwarnings("ignore")
    nneighbors = [5, 10, 15, 20, 50]  # used nn=10 in the paper
    mindists = [
        0,
        0.01,
        0.02,
        0.05,
        0.1,
        0.2,
        0.5,
        1,
    ]  # used 0.5 (default) in the paper
    nbinses = [50, 100, 200, 300, 400, 500]
    do_normalize = False
    plt.close()
    for nbins in nbinses:
        wp_max_pol, wp_binned_values, xvals = MovingAverages.bin_values(
            nbins,
            protein_ccd.u_well_plates,
            protein_ccd.fucci_polar_coords.pol_sort_norm_rev,
            protein_ccd.fucci_polar_coords.pol_sort_well_plate,
            protein_ccd.fucci_polar_coords.pol_sort_ab_cell,
            protein_ccd.fucci_polar_coords.pol_sort_ab_nuc,
            protein_ccd.fucci_polar_coords.pol_sort_ab_cyto,
            protein_ccd.fucci_polar_coords.pol_sort_mt_cell,
            protein_ccd.wp_iscell,
            protein_ccd.wp_isnuc,
            protein_ccd.wp_iscyto,
            do_normalize=do_normalize,
        )
        wp_binned_values = np.array(wp_binned_values)
        for nn in nneighbors:
            for md in mindists:
                cutoff = protein_ccd.chosen_cutoff
                reducer = umap.UMAP(n_neighbors=nn, min_dist=md, random_state=0)
                embeddingCcd = reducer.fit_transform(
                    wp_binned_values[protein_ccd.mean_diff_from_rng > cutoff, :].T
                )
                plt.scatter(embeddingCcd[:, 0], embeddingCcd[:, 1], c=xvals[1:])
                plt.xlabel("UMAP1")
                plt.ylabel("UMAP2")
                cb = plt.colorbar()
                cb.set_label("Pseudotime")
                plt.savefig(
                    f"figures/ProteinUmapStability/proteinUmap_nbins{nbins}_nn{nn}_md{md}_{cutoff}CCD.pdf"
                )
                plt.close()
                embeddingNonCcd = reducer.fit_transform(
                    wp_binned_values[protein_ccd.mean_diff_from_rng <= cutoff, :].T
                )
                plt.scatter(embeddingNonCcd[:, 0], embeddingNonCcd[:, 1], c=xvals[1:])
                plt.xlabel("UMAP1")
                plt.ylabel("UMAP2")
                cb = plt.colorbar()
                cb.set_label("Pseudotime")
                plt.savefig(
                    f"figures/ProteinUmapStability/proteinUmap_nbins{nbins}_nn{nn}_md{md}_{cutoff}NonCCD.pdf"
                )
                plt.close()

    chosen_nb = BINS_FOR_UMAP_AND_LASSO
    if not os.path.exists("figures/ProteinUmapNumBins"):
        os.mkdir("figures/ProteinUmapNumBins")
    for nbins in nbinses:
        orig = f"figures/ProteinUmapStability/proteinUmap_nbins{nbins}_nn{chosen_nn}_md{chosen_md}_{protein_ccd.chosen_cutoff}CCD.pdf"
        new = f"figures/ProteinUmapNumBins/proteinUmap_nbins{nbins}_nn{chosen_nn}_md{chosen_md}_{protein_ccd.chosen_cutoff}CCD.pdf"
        shutil.copy(orig, new)
    if not os.path.exists("figures/ProteinUmapStabilityChoice"):
        os.mkdir("figures/ProteinUmapStabilityChoice")
    for nn in nneighbors:
        for md in mindists:
            orig = f"figures/ProteinUmapStability/proteinUmap_nbins{chosen_nb}_nn{nn}_md{md}_{protein_ccd.chosen_cutoff}CCD.pdf"
            new = f"figures/ProteinUmapStabilityChoice/proteinUmap_nbins{chosen_nb}_nn{nn}_md{md}_{protein_ccd.chosen_cutoff}CCD.pdf"
            shutil.copy(orig, new)

    nbins = BINS_FOR_UMAP_AND_LASSO  # Seems to be the most stable given all genes. When the subsets get small, though, it'll fall apart.
    wp_max_pol, wp_binned_values, xvals = MovingAverages.bin_values(
        nbins,
        protein_ccd.u_well_plates,
        protein_ccd.fucci_polar_coords.pol_sort_norm_rev,
        protein_ccd.fucci_polar_coords.pol_sort_well_plate,
        protein_ccd.fucci_polar_coords.pol_sort_ab_cell,
        protein_ccd.fucci_polar_coords.pol_sort_ab_nuc,
        protein_ccd.fucci_polar_coords.pol_sort_ab_cyto,
        protein_ccd.fucci_polar_coords.pol_sort_mt_cell,
        protein_ccd.wp_iscell,
        protein_ccd.wp_isnuc,
        protein_ccd.wp_iscyto,
        do_normalize=do_normalize,
    )
    wp_binned_values = np.array(wp_binned_values)
    reducer = umap.UMAP(n_neighbors=chosen_nn, min_dist=chosen_md, random_state=0)
    for cutoff in (np.arange(20) + 1) / 100:
        embeddingCcd = reducer.fit_transform(
            wp_binned_values[protein_ccd.mean_diff_from_rng > cutoff, :].T
        )

        plt.scatter(embeddingCcd[:, 0], embeddingCcd[:, 1], c=xvals[1:])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        cb = plt.colorbar()
        cb.set_label("Pseudotime")
        plt.savefig(f"figures/ProteinUmaps/proteinUmap{round(cutoff, 2)}CCD.pdf")
        plt.close()

        embeddingNonCcd = reducer.fit_transform(
            wp_binned_values[protein_ccd.mean_diff_from_rng <= cutoff, :].T
        )
        plt.scatter(embeddingNonCcd[:, 0], embeddingNonCcd[:, 1], c=xvals[1:])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        cb = plt.colorbar()
        cb.set_label("Pseudotime")
        plt.savefig(f"figures/ProteinUmaps/proteinUmap{round(cutoff, 2)}NonCCD.pdf")
        plt.close()
    warnings.filterwarnings("default")


def make_plotting_dataframe(protein_ccd):
    """Make a single table for HPA website figures on protein pseudotime, boxplots, and fucci plots"""
    mvperc_10p = [x[0] for x in protein_ccd.mvperc_comps]
    mvperc_90p = [x[-1] for x in protein_ccd.mvperc_comps]
    mvperc_25p = [x[1] for x in protein_ccd.mvperc_comps]
    mvperc_75p = [x[-2] for x in protein_ccd.mvperc_comps]
    ccdStrings = utils.get_ccd_strings(
        protein_ccd.ccd_comp, protein_ccd.wp_ensg, protein_ccd.bioccd
    )

    # Save plotting dataframe
    compartment = utils.get_compartment_strings(
        protein_ccd.wp_iscell, protein_ccd.wp_iscyto, protein_ccd.wp_isnuc
    )
    cell_pseudotime = [
        ",".join([str(ppp) for ppp in pp]) for pp in protein_ccd.curr_pols
    ]
    cell_intensity = [
        ",".join([str(yyy) for yyy in yy]) for yy in protein_ccd.curr_ab_norms
    ]
    cell_fred = [",".join([str(rrr) for rrr in rr]) for rr in protein_ccd.curr_freds]
    cell_fgreen = [
        ",".join([str(ggg) for ggg in gg]) for gg in protein_ccd.curr_fgreens
    ]
    mvavg_x = [",".join([str(xxx) for xxx in xx]) for xx in protein_ccd.mvavgs_x]
    mvavg_y = [",".join([str(yyy) for yyy in yy]) for yy in protein_ccd.mvavgs_comp]
    df = pd.DataFrame(
        {
            "ENSG": protein_ccd.wp_ensg,
            "Antibody": protein_ccd.wp_ab,
            "Compartment": compartment,
            "CCD": ccdStrings,
            "cell_pseudotime": cell_pseudotime,
            "cell_intensity": cell_intensity,
            "cell_fred": cell_fred,
            "cell_fgreen": cell_fgreen,
            "mvavg_x": mvavg_x,
            "mvavg_y": mvavg_y,
            "mvavgs_10p": [",".join([str(yyy) for yyy in yy]) for yy in mvperc_10p],
            "mvavgs_90p": [",".join([str(yyy) for yyy in yy]) for yy in mvperc_90p],
            "mvavgs_25p": [",".join([str(yyy) for yyy in yy]) for yy in mvperc_25p],
            "mvavgs_75p": [",".join([str(yyy) for yyy in yy]) for yy in mvperc_75p],
            "phase": [",".join(pp) for pp in protein_ccd.curr_mockbulk_phases],
            "gini": protein_ccd.protein_variability.gini_comp,
            "percent_variance": protein_ccd.percvar_comp,
            "WellPlate": protein_ccd.u_well_plates,
        }
    )
    only_these_no_mt = ~np.isin(
        protein_ccd.u_well_plates, protein_ccd.removeThese
    ) & np.array([not xx.startswith("Mitotic") for xx in ccdStrings])
    df_filtered = df[only_these_no_mt]
    df_filtered.to_csv("output/ProteinPseudotimePlotting.csv.gz", index=False, sep="\t")

    # Save information on mitotic proteins
    df_mt = pd.DataFrame(
        {
            "ENSG": protein_ccd.wp_ensg,
            "Antibody": protein_ccd.wp_ab,
            "Compartment": utils.get_compartment_strings(
                protein_ccd.wp_iscell, protein_ccd.wp_iscyto, protein_ccd.wp_isnuc
            ),
            "CCD": ccdStrings,
            "WellPlate": protein_ccd.u_well_plates,
        }
    )
    only_these_with_mt = ~np.isin(
        protein_ccd.u_well_plates, protein_ccd.removeThese
    ) & np.array([xx.startswith("Mitotic") for xx in ccdStrings])
    df_mt_filtered = df_mt[only_these_with_mt]
    df_mt_filtered.to_csv("output/ProteinMitoticOnly.csv.gz", index=False, sep="\t")


def additional_ccd_analysis(protein_ccd):
    """Additional analysis for after the CCD calculations"""
    n_tot_variable = len(protein_ccd.u_well_plates)
    general_statements = [
        f"{n_tot_variable}: # total samples",
        f"{sum(protein_ccd.ccd_comp)}: CCD variable proteins (before addressing redundancy and mitotic structures)",
        f"{sum(protein_ccd.nonccd_comp)}: non-CCD variable proteins",
    ]
    print("\n".join(general_statements))

    ### address gene redundancy
    wp_ccd_bimodalonecluster = (
        protein_ccd.wp_comp_ccd_clust1 ^ protein_ccd.wp_comp_ccd_clust2
    )
    wp_ccd_bimodaltwocluster = (
        protein_ccd.wp_comp_ccd_clust1 & protein_ccd.wp_comp_ccd_clust2
    )
    wp_removeReplicate = np.isin(protein_ccd.u_well_plates, protein_ccd.removeThese)
    protein_ct = len(
        np.unique(np.concatenate((protein_ccd.wp_ensg, protein_ccd.bioccd)))
    )
    ccd_protein_ct = sum(protein_ccd.wp_ccd_unibimodal[~wp_removeReplicate])
    nonccd_protein_ct = sum(~protein_ccd.wp_ccd_unibimodal[~wp_removeReplicate])
    unimodal_generally_protein_ct = sum(
        ~protein_ccd.wp_isbimodal_generally[~wp_removeReplicate]
    )
    bimodal_generally_protein_ct = sum(
        protein_ccd.wp_isbimodal_generally[~wp_removeReplicate]
    )
    print(
        f"{ccd_protein_ct}: number of ccd proteins; addressed replicates; not including mitotic structures"
    )

    # Decision: remove replicate antibodies manually
    protein_ccd.ccd_comp[wp_removeReplicate] = False
    protein_ccd.nonccd_comp[wp_removeReplicate] = False

    # Accounting for biologically CCD ones
    knownccd1 = np.genfromtxt(
        "input/ProteinData/knownccd.txt", dtype="str"
    )  # from gene ontology, reactome, cyclebase 3.0, NCBI gene from mcm3
    knownccd2 = np.genfromtxt(
        "input/ProteinData/known_go_ccd.txt", dtype="str"
    )  # from GO cell cycle
    knownccd3 = np.genfromtxt(
        "input/ProteinData/known_go_proliferation.txt", dtype="str"
    )  # from GO proliferation
    print(f"{len(protein_ccd.bioccd)}: number of mitotic structure proteins")

    ccd_prots_withmitotic = np.unique(
        np.concatenate((protein_ccd.wp_ensg[protein_ccd.ccd_comp], protein_ccd.bioccd))
    )
    total_proteins_minusmitotic = len(ccd_prots_withmitotic) - sum(
        np.isin(protein_ccd.wp_ensg[protein_ccd.nonccd_comp], protein_ccd.bioccd)
    )
    total_ccd_proteins_withmitotic = len(ccd_prots_withmitotic)
    total_nonccd_proteins_minusmitotic = nonccd_protein_ct - sum(
        np.isin(protein_ccd.wp_ensg[protein_ccd.nonccd_comp], protein_ccd.bioccd)
    )
    overlapping_knownccd1 = sum(
        np.isin(
            ccd_prots_withmitotic, np.concatenate((knownccd1, knownccd2, knownccd3))
        )
    )

    with open("output/figuresofmerit.txt", "w") as file:
        fom = "--- protein pseudotime\n\n"
        protCount = len(np.unique(protein_ccd.wp_ensg))
        fom += f"{protCount} proteins that were expressed and exhibited variations in the U-2 OS cell line were selected"
        fom += "\n\n"

        novelCount = len(ccd_prots_withmitotic) - overlapping_knownccd1
        fom += f"present the first evidence of cell cycle association for {novelCount} proteins"
        fom += "\n\n"

        ccdPercent = 100 * ccd_protein_ct / len(np.unique(protein_ccd.wp_ensg))
        fom += f"Based on this analysis, we identified {ccd_protein_ct} out of {protCount} proteins ({ccdPercent}%) to have variance in expression levels temporally correlated to cell cycle progression, "
        fom += f"and for which the cell-cycle explained {protein_ccd.chosen_cutoff}% or more variance in expression than random."
        fom += "\n\n"

        nonCcdPercent = (
            100
            * total_nonccd_proteins_minusmitotic
            / len(np.unique(protein_ccd.wp_ensg))
        )
        fom += f"majority of the proteins analyzed ({total_nonccd_proteins_minusmitotic}, {nonCcdPercent}%) showed cell-to-cell variations that were largely unexplained by cell cycle progression"
        fom += "\n\n"
        bothPseudotimeMitoticCount = (
            ccd_protein_ct + len(protein_ccd.bioccd) - total_ccd_proteins_withmitotic
        )
        knownCcdPercent = 100 * overlapping_knownccd1 / total_ccd_proteins_withmitotic
        novelCcdPercent = (
            100
            * (total_ccd_proteins_withmitotic - overlapping_knownccd1)
            / total_ccd_proteins_withmitotic
        )
        fom += f"Of the {total_ccd_proteins_withmitotic} proteins ({ccd_protein_ct} in interphase, {len(protein_ccd.bioccd)} in mitotic structures"
        fom += f"and {bothPseudotimeMitoticCount} in both sets) identified to correlate to cell cycle progression, {overlapping_knownccd1} ({knownCcdPercent}%) "
        fom += "had a known association to the cell cycle as determined either by a GO BP term ... "
        fom += f"The remaining {total_ccd_proteins_withmitotic - overlapping_knownccd1} proteins ({novelCcdPercent}%),"
        fom += "\n\n"

        fom += f"The patterns of variability were investigated for these {sum(~wp_removeReplicate)} proteins for the population of cells measured for each protein. "
        fom += f"The mean fold change between the highest and lowest expressing cells per protein was {np.mean(np.array(protein_ccd.protein_bimodality.wp_bimodal_fcmaxmin)[~wp_removeReplicate])}."
        fom += "\n\n"

        fom += f"We determined that {unimodal_generally_protein_ct} proteins ({100 * unimodal_generally_protein_ct / len(np.unique(protein_ccd.wp_ensg))}%) had unimodal intensity distributions, "
        fom += f"and {bimodal_generally_protein_ct} proteins ({100 * bimodal_generally_protein_ct / len(np.unique(protein_ccd.wp_ensg))}%) were found to display bimodality"
        fom += "\n\n"

        fom += f"Of {sum(protein_ccd.protein_bimodality.wp_isbimodal_fcpadj_pass)} bimodal samples that were analyzed for cell cycle dependence, {sum(protein_ccd.wp_ccd_bimodalonecluster)} were CCD in one cluster "
        fom += f"({sum(protein_ccd.wp_ccd_bimodalonecluster & protein_ccd.wp_comp_ccd_difffromrng)} of these were CCD when analyzed unimodally), and {sum(protein_ccd.wp_ccd_bimodaltwocluster)} were CCD in both clusters "
        fom += f"({sum(protein_ccd.wp_ccd_bimodaltwocluster & protein_ccd.wp_comp_ccd_difffromrng)} were also CCD when analyzed unimodally), and the remaining "
        fom += f"{sum(protein_ccd.protein_bimodality.wp_isbimodal_fcpadj_pass & ~protein_ccd.wp_ccd_bimodalonecluster & ~protein_ccd.wp_ccd_bimodaltwocluster)} were non-CCD in both clusters."
        fom += "\n\n"

        ccdPseudotimeNotGauss = sum(
            protein_ccd.ccd_comp & ~protein_ccd.protein_clustering.wp_comp_ccd_gauss
        )
        fom += f"We aggregated single-cell measurements by cell cycle phase to simulate a bulk experiment, and {ccdPseudotimeNotGauss} CCD proteins detected at the single-cell level were not detectable in bulk phases, such as TRNT1"
        print(fom)
        file.write(fom)

    # read in reliability scores
    wp_ab_list = list(protein_ccd.wp_ab)
    ab_scores = list(np.zeros(protein_ccd.wp_ab.shape, dtype=str))
    with open("input/ProteinData/ReliabilityScores.txt") as file:
        for line in file:
            if line.startswith("Antibody RRID"):
                continue
            score = line.split("\t")[1].strip()
            ablist = line.split("\t")[0].replace(":", "").replace(",", "").split()
            for ab in ablist:
                if ab in protein_ccd.wp_ab:
                    ab_scores[wp_ab_list.index(ab)] = score

    pd.DataFrame(
        {
            "well_plate": protein_ccd.u_well_plates,
            "ENSG": protein_ccd.wp_ensg,
            "antibody": protein_ccd.wp_ab,
            "antibody_hpa_scores": ab_scores,
            "compartment": utils.get_compartment_strings(
                protein_ccd.wp_iscell, protein_ccd.wp_iscyto, protein_ccd.wp_isnuc
            ),
            "variance_comp": protein_ccd.protein_variability.var_comp,
            "gini_comp": protein_ccd.protein_variability.gini_comp,
            "percent_variance_explained": protein_ccd.percvar_comp,
            "known_by_GoReactomeCyclebaseNcbi": np.isin(
                protein_ccd.wp_ensg, np.concatenate((knownccd1, knownccd2, knownccd3))
            ),
            "mean_percvar_diff_from_random": protein_ccd.mean_diff_from_rng,
            "wp_comp_kruskal_gaussccd_adj": protein_ccd.protein_clustering.wp_comp_kruskal_gaussccd_adj,
            "log_10_pval_eq_percvar": -np.log10(
                np.nextafter(
                    protein_ccd.wp_comp_eq_percvar_adj,
                    protein_ccd.wp_comp_eq_percvar_adj + 1,
                )
            ),
            "pass_median_diff": protein_ccd.wp_comp_ccd_difffromrng,
            "pass_gauss": protein_ccd.protein_clustering.wp_comp_ccd_gauss,
            "CCD_COMP": protein_ccd.ccd_comp,
            "ccd_reason": utils.get_ccd_strings(
                protein_ccd.ccd_comp, protein_ccd.wp_ensg, protein_ccd.bioccd
            ),
            "nonccd_comp": protein_ccd.nonccd_comp,
            # bimodal significance testing
            "ccd_unimodal": protein_ccd.wp_comp_ccd_difffromrng,
            "ccd_clust1": protein_ccd.wp_comp_ccd_clust1,
            "clust1_difffromrng": protein_ccd.mean_diff_from_rng_clust1,
            "clust1_log10pval_percvar": -np.log10(
                np.nextafter(
                    protein_ccd.wp_comp_eq_percvar_adj_clust1,
                    protein_ccd.wp_comp_eq_percvar_adj_clust1 + 1,
                )
            ),
            "ccd_clust2": protein_ccd.wp_comp_ccd_clust2,
            "ccd_clust2_difffromrng": protein_ccd.mean_diff_from_rng_clust2,
            "ccd_clust2_log10pval_percvar": -np.log10(
                np.nextafter(
                    protein_ccd.wp_comp_eq_percvar_adj_clust2,
                    protein_ccd.wp_comp_eq_percvar_adj_clust2 + 1,
                )
            ),
        }
    ).to_csv("output/CellCycleVariationSummary.csv", index=False)
