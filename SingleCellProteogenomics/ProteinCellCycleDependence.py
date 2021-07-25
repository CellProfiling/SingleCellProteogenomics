# -*- coding: utf-8 -*-
"""
Methods for assessing cell cycle dependence of protein abundance in single cells.
-  Percent variance attributed to the cell cycle was calculated using the (variance of moving average / total variance)
-  Randomization analysis was used to determine statistical significance of high percent variances due to the cell cycle

@author: Anthony J. Cesnik, cesnik@stanford.edu
@author: devinsullivan
"""

from SingleCellProteogenomics import (
    utils,
    FucciPseudotime,
    ProteinVariability,
    ProteinBimodality,
    ProteinCellCycleMovingAverageResult,
)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import shutil, os

# Cutoff used for percent additional variance explained by the cell cycle than random
MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM = 0.08


def clust_to_wp(clust, clust_idx):
    """
    clust: boolean-typed numpy array
    Gather results for either high- and low-expressing cell populations from combined results.
    (The high- and low-expressing populations were combined with all single-population results for multiple testing correction.)
    """
    wp_clust = np.array([False] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust


def clust_to_wp_doub(clust, clust_idx):
    """
    clust: double-typed numpy array
    Gather results for either high- and low-expressing cell populations from combined results.
    (The high- and low-expressing populations were combined with all single-population results for multiple testing correction.)
    """
    wp_clust = np.array([0.0] * len(clust_idx))
    wp_clust[clust_idx] = clust
    return wp_clust


def get_fileprefixes(wp_ensg):
    """Generate the file prefixes for given genes"""
    return np.array(
        [f"{ensg}_{sum(wp_ensg[:ei] == ensg)}" for ei, ensg in enumerate(wp_ensg)]
    )


class ProteinCellCycleDependence:
    def __init__(self, protein_data, protein_clustering, do_plotting):
        self.fucci_polar_coords = FucciPseudotime.FucciPolarCoords()
        self.protein_variability = ProteinVariability.ProteinVariability()
        self.protein_bimodality = ProteinBimodality.ProteinBimodality()

        # load duplicate antibodies marked for removal
        self.removeThese = pd.read_csv(
            "input/ProteinData/ReplicatesToRemove.txt", header=None
        )[0]

        self.u_well_plates = protein_data.u_well_plates
        self.wp_ensg = protein_data.wp_ensg
        self.wp_ab = protein_data.wp_ab
        self.ab_nuc = protein_data.ab_nuc
        self.ab_cyto = protein_data.ab_cyto
        self.ab_cell = protein_data.ab_cell
        self.mt_cell = protein_data.mt_cell
        self.area_cell = protein_data.area_cell
        self.area_nuc = protein_data.area_nuc
        self.well_plate = protein_data.well_plate
        self.well_plate_imgnb = protein_data.well_plate_imgnb
        self.well_plate_imgnb_objnb = protein_data.well_plate_imgnb_objnb
        self.log_red_fucci_zeroc_rescale = (
            protein_clustering.log_red_fucci_zeroc_rescale
        )
        self.log_green_fucci_zeroc_rescale = (
            protein_clustering.log_green_fucci_zeroc_rescale
        )
        self.mockbulk_phases = protein_clustering.mockbulk_phases
        self.wp_iscell = protein_data.wp_iscell
        self.wp_isnuc = protein_data.wp_isnuc
        self.wp_iscyto = protein_data.wp_iscyto
        self.fucci_data = protein_clustering.fucci_data

        self.chosen_cutoff = MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM

        self.do_plotting = do_plotting

        self.folder = "figures/TemporalMovingAverages"
        self.folder_mt = "figures/TemporalMovingAveragesMicrotubules"
        self.folder_rng = "figures/TemporalMovingAverageRandomizationExamples"

        self.fileprefixes = get_fileprefixes(self.wp_ensg)
        self.prepare_for_plotting()

    def calculate_cell_cycle_dependence_protein(
        self,
        threads,
        alpha_ccd=0.01,
        use_log_ccd=False,
        keep_outliers=False,
    ):
        """
        Use a moving average model of protein expression over the cell cycle to determine cell cycle dependence.
        Generates plots for each gene
        """
        self.use_log_ccd = use_log_ccd
        self.do_remove_outliers = not keep_outliers

        # Moving average calculations
        self.curr_well_inds = [
            self.fucci_polar_coords.pol_sort_well_plate == well
            for well in self.u_well_plates
        ]
        moving_avg_results = [
            ProteinCellCycleMovingAverageResult.ProteinCellCycleMovingAverageResult(self, idx)
            for idx, well in enumerate(self.u_well_plates)
        ]
        utils.parmap(
            ProteinCellCycleMovingAverageResult.ProteinCellCycleMovingAverageResult.calculate_moving_averages,
            moving_avg_results,
            threads,
        )
        self.parse_moving_average_results(moving_avg_results)
        self.analyze_ccd_variation_protein(alpha_ccd)

    def parse_moving_average_results(self, results):
        """Move the moving average results to separate lists"""
        # For plotting dataframe
        self.curr_freds = np.array([m.curr_fred for m in results])
        self.curr_fgreens = np.array([m.curr_fgreen for m in results])
        self.curr_mockbulk_phases = np.array([m.curr_mockbulk_phase for m in results])

        # percent variances for bimodal
        self.perc_var_comp_clust1 = np.array(
            [m.perc_var_comp_clust1_val for m in results]
        )
        self.perc_var_comp_clust2 = np.array(
            [m.perc_var_comp_clust2_val for m in results]
        )
        self.mvavgs_comp_clust1 = np.array([m.mvavg_clust1 for m in results])
        self.mvavgs_comp_clust2 = np.array([m.mvavg_clust2 for m in results])
        self.perc_var_comp_clust1_rng = np.array(
            [m.curr_clust1_percvar_rng_comp for m in results]
        )
        self.perc_var_comp_clust2_rng = np.array(
            [m.curr_clust2_percvar_rng_comp for m in results]
        )
        self.mvavgs_x_clust1 = np.array([m.mvavg_x_clust1 for m in results])
        self.mvavgs_x_clust2 = np.array([m.mvavg_x_clust2 for m in results])

        # moving average y values & x value
        self.mvavgs_comp = np.array([m.mvavg_comp for m in results])
        self.mvavgs_mt = np.array([m.mvavg_mt for m in results])
        self.mvavgs_x = np.array([m.mvavg_xvals for m in results])

        # percent variance attributed to cell cycle (mean POI intensities)
        self.perc_var_cell = np.array([m.curr_perc_var_cell for m in results])
        self.perc_var_nuc = np.array([m.curr_perc_var_nuc for m in results])
        self.perc_var_cyto = np.array([m.curr_perc_var_cyto for m in results])
        self.perc_var_mt = np.array([m.curr_perc_var_mt for m in results])
        self.perc_var_comp = utils.values_comp(
            self.perc_var_cell,
            self.perc_var_nuc,
            self.perc_var_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )

        # randomized in pseudotime; percent variances
        self.perc_var_comp_rng = np.array([m.curr_percvar_rng_comp for m in results])
        self.perc_var_mt_rng = np.array([m.curr_percvar_rng_mt for m in results])

        self.curr_ab_norms = np.array([m.curr_comp_norm for m in results])
        self.mvperc_comps = np.array([m.mvperc_comp for m in results])

        # Get other values with the current inds
        self.cell_counts = np.array([len(m.curr_pol) for m in results])
        self.curr_area_cell = np.array(
            [
                self.fucci_polar_coords.pol_sort_area_cell[self.curr_well_inds[idx]]
                for idx, m in enumerate(results)
            ]
        )
        self.curr_area_nuc = np.array(
            [
                self.fucci_polar_coords.pol_sort_area_nuc[self.curr_well_inds[idx]]
                for idx, m in enumerate(results)
            ]
        )
        self.curr_well_plate_imgnb = np.array(
            [
                self.fucci_polar_coords.pol_sort_well_plate_imgnb[
                    self.curr_well_inds[idx]
                ]
                for idx, m in enumerate(results)
            ]
        )

    def prepare_for_plotting(self):
        """Prepares folders for plotting results"""
        folders = [self.folder, self.folder_mt, self.folder_rng]
        for folder in folders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def analyze_ccd_variation_protein(self, alpha_ccd):
        """Analyze the cell cycle dependence for all the proteins"""
        perc_var_comp_withbimodal = np.concatenate(
            (self.perc_var_comp, self.perc_var_comp_clust1, self.perc_var_comp_clust2)
        )
        perc_var_comp_rng_withbimodal = np.concatenate(
            (
                self.perc_var_comp_rng,
                self.perc_var_comp_clust1_rng,
                self.perc_var_comp_clust2_rng,
            )
        )
        perc_diff_from_rng_withbimodal = (perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T
        perc_diff_from_rng_withbimodal_mt = (self.perc_var_mt - self.perc_var_mt_rng.T).T
        ccd_var_comp_rng_wilcoxp_withbimodal = (
            np.apply_along_axis(
                scipy.stats.wilcoxon,
                1,
                perc_diff_from_rng_withbimodal,
                None,
                "wilcox",
                False,
                "greater",
            )
            .T[1]
            .T
        )
        ccd_var_mt_rng_wilcoxp = (
            np.apply_along_axis(
                scipy.stats.wilcoxon,
                1,
                perc_diff_from_rng_withbimodal_mt,
                None,
                "wilcox",
                False,
                "greater",
            )
            .T[1]
            .T
        )
        mean_diff_from_rng_mt = np.mean(
            perc_diff_from_rng_withbimodal_mt, 1
        )

        # randomization tests, try being a bit more stringent, try drawing the cutoff based on microtubules per sample
        (
            wp_comp_eq_percvar_adj_withbimodal,
            wp_comp_pass_eq_percvar_adj_withbimodal,
        ) = utils.bonf(alpha_ccd, ccd_var_comp_rng_wilcoxp_withbimodal)
        wp_comp_gtpass_eq_percvar_adj_withbimodal = (
            wp_comp_pass_eq_percvar_adj_withbimodal
            & (
                perc_var_comp_withbimodal
                > np.median(perc_var_comp_rng_withbimodal, axis=1)
            )
        )

        # median differences from random
        req_blurb = f"Requiring {MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM*100}% additional percent variance explained than random."
        print(req_blurb)
        mean_diff_from_rng_withbimodal = np.mean(
            perc_diff_from_rng_withbimodal, 1
        )
        wp_comp_ccd_difffromrng_withbimodal = (
            mean_diff_from_rng_withbimodal >= MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM
        )

        ###### Calculate the cutoffs for total intensity and percent variance attributed to the cell cycle
        # separate unimodal from bimodal again
        end = len(self.perc_var_comp)
        end2 = end + len(self.perc_var_comp_clust1)
        wp_comp_eq_percvar_adj = wp_comp_eq_percvar_adj_withbimodal[:end]
        wp_comp_pass_eq_percvar_adj = wp_comp_pass_eq_percvar_adj_withbimodal[:end]
        wp_comp_gtpass_eq_percvar_adj = wp_comp_gtpass_eq_percvar_adj_withbimodal[:end]
        mean_diff_from_rng = mean_diff_from_rng_withbimodal[:end]
        wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng_withbimodal[:end]

        wp_comp_pass_eq_percvar_adj_clust1 = clust_to_wp(
            wp_comp_pass_eq_percvar_adj_withbimodal[end:end2],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        wp_comp_eq_percvar_adj_clust1 = clust_to_wp(
            wp_comp_eq_percvar_adj_withbimodal[end:end2],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        mean_diff_from_rng_clust1 = clust_to_wp_doub(
            mean_diff_from_rng_withbimodal[end:end2],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        wp_comp_ccd_difffromrng_clust1 = clust_to_wp(
            wp_comp_ccd_difffromrng_withbimodal[end:end2],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        wp_comp_pass_eq_percvar_adj_clust2 = clust_to_wp(
            wp_comp_pass_eq_percvar_adj_withbimodal[end2:],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        wp_comp_eq_percvar_adj_clust2 = clust_to_wp(
            wp_comp_eq_percvar_adj_withbimodal[end2:],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        mean_diff_from_rng_clust2 = clust_to_wp_doub(
            mean_diff_from_rng_withbimodal[end2:],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )
        wp_comp_ccd_difffromrng_clust2 = clust_to_wp(
            wp_comp_ccd_difffromrng_withbimodal[end2:],
            self.protein_bimodality.wp_isbimodal_fcpadj_pass,
        )

        wp_normal_randompercvar_p = (
            np.apply_along_axis(
                scipy.stats.normaltest,
                1,
                (self.perc_var_comp - self.perc_var_comp_rng.T).T,
            )
            .T[1]
            .T
        )
        wp_randompercvarnorm_adj, wp_randompercvarnorm_pass = utils.benji_hoch(
            0.05, wp_normal_randompercvar_p
        )

        self.wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng
        wp_comp_ccd_gauss = (
            self.protein_bimodality.wp_comp_kruskal_gaussccd_adj <= self.alpha_ccd
        )
        percvar_statements = [
            f"{sum(wp_randompercvarnorm_pass)}: number of genes with randomized percvars that form normal distributions",
            f"{sum(wp_comp_ccd_difffromrng)}: # proteins showing CCD variation unimodally, comp, percvar rng median diff",
            f"{sum(wp_comp_ccd_difffromrng) / len(wp_comp_ccd_difffromrng)}: fraction of variable proteins showing CCD variation, comp, percvar rng median diff",
            f"{sum(wp_comp_ccd_gauss)}: # proteins showing CCD variation, comp, gaussian analysis",
            f"{sum(wp_comp_ccd_gauss) / len(wp_comp_ccd_gauss)}: fraction of variable proteins showing CCD variation, comp, gaussian analysis",
            f"{sum(wp_comp_ccd_gauss & wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, CCD and gaussian",
            f"{sum(wp_comp_ccd_gauss & ~wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, not CCD and gaussian",
            f"{sum(~wp_comp_ccd_gauss & wp_comp_ccd_difffromrng)}: # proteins showing CCD variation, CCD and not gaussian",
        ]
        print("\n".join(percvar_statements))

        # Address bimodal ones
        # 1) The number that are bimodal in one cluster (number of those that are also CCD as unimodal)
        # 2) The number that are bimodal in both clusters (number of those that are also CCD as unimodal)
        wp_comp_ccd_clust1 = wp_comp_ccd_difffromrng_clust1
        wp_comp_ccd_clust2 = wp_comp_ccd_difffromrng_clust2
        bimodal_statements = [
            f"{sum(self.protein_bimodality.wp_isbimodal_fcpadj_pass)}: samples with bimodal antibody intensities",
            f"{sum(wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2)}: bimodal samples with one CCD cluster ({sum((wp_comp_ccd_clust1 ^ wp_comp_ccd_clust2) & wp_comp_ccd_difffromrng)}: also CCD unimodally)",
            f"{sum(wp_comp_ccd_clust1 & wp_comp_ccd_clust2)}: bimodal samples with two CCD clusters ({sum((wp_comp_ccd_clust1 & wp_comp_ccd_clust2) & wp_comp_ccd_difffromrng)}: also CCD unimodally)",
        ]
        print("\n".join(bimodal_statements))

        self.wp_ccd_unibimodal = (
            wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
        )

        ccd_comp = wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
        nonccd_comp = ~ccd_comp
        bioccd = np.genfromtxt(
            "input/ProteinData/BiologicallyDefinedCCD.txt", dtype="str"
        )  # from mitotic structures

        self.wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng
        self.mean_diff_from_rng_mt = mean_diff_from_rng_mt
        self.wp_comp_ccd_clust1 = wp_comp_ccd_clust1
        self.wp_comp_ccd_clust2 = wp_comp_ccd_clust2
        self.wp_comp_ccd_gauss = wp_comp_ccd_gauss
        self.mean_diff_from_rng = mean_diff_from_rng
        self.wp_comp_eq_percvar_adj = wp_comp_eq_percvar_adj
        self.mean_diff_from_rng_clust1 = mean_diff_from_rng_clust1
        self.wp_comp_eq_percvar_adj_clust1 = wp_comp_eq_percvar_adj_clust1
        self.mean_diff_from_rng_clust2 = mean_diff_from_rng_clust2
        self.wp_comp_eq_percvar_adj_clust2 = wp_comp_eq_percvar_adj_clust2

        self.ccd_comp = ccd_comp
        self.nonccd_comp = nonccd_comp
        self.bioccd = bioccd

        # pickle the results (removed ones passing in only one replicate)
        utils.np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp)
        utils.np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp)

        # Do some plotting if desired
        if self.do_plotting:
            self.protein_ccd.copy_mvavg_plots_protein()

            plt.figure(figsize=(10, 10))
            plt.scatter(
                perc_var_comp_withbimodal,
                mean_diff_from_rng_withbimodal,
                c=wp_comp_ccd_difffromrng_withbimodal,
                cmap="bwr_r",
            )
            plt.hlines(
                MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM,
                np.min(self.perc_var_comp),
                np.max(self.perc_var_comp),
                color="gray",
            )
            plt.xlabel("Percent Variance Explained by Cell Cycle")
            plt.ylabel("Mean Difference from Random")
            plt.savefig("figures/MedianDiffFromRandom.png")
            plt.savefig("figures/MedianDiffFromRandom.pdf")
            plt.close()

            pervar_adj_withbimodal_nextafter = np.nextafter(
                wp_comp_eq_percvar_adj_withbimodal,
                wp_comp_eq_percvar_adj_withbimodal + 1,
            )
            plt.figure(figsize=(10, 10))
            plt.scatter(
                mean_diff_from_rng_withbimodal,
                -np.log10(pervar_adj_withbimodal_nextafter),
                c=wp_comp_ccd_difffromrng_withbimodal,
                cmap="bwr_r",
            )
            plt.vlines(
                MIN_MEAN_PERCVAR_DIFF_FROM_RANDOM,
                np.min(-np.log10(pervar_adj_withbimodal_nextafter)),
                np.max(-np.log10(pervar_adj_withbimodal_nextafter)),
                color="gray",
            )
            plt.xlabel("Mean Difference from Random")
            plt.ylabel("-log10 adj p-value from randomization")
            plt.savefig("figures/MedianDiffFromRandomVolcano.png")
            plt.savefig("figures/MedianDiffFromRandomVolcano.pdf")
            plt.close()

    def copy_mvavg_plots_protein(self):
        """Copy the plots generated for each gene to a more informative location"""
        folder = self.folder
        wp_comp_ccd_difffromrng = self.wp_comp_ccd_difffromrng
        wp_isbimodal_fcpadj_pass = self.protein_bimodality.wp_isbimodal_fcpadj_pass
        wp_comp_ccd_clust1 = self.wp_comp_ccd_clust1
        wp_comp_ccd_clust2 = self.wp_comp_ccd_clust2
        wp_ccd_unibimodal = self.wp_ccd_unibimodal
        wp_comp_ccd_gauss = self.wp_comp_ccd_gauss

        # Copy profiles to the right place:
        # 1) CCD Unimodal
        # 2) CCD Bimodal
        # 3) Non-CCD
        # 3) Gaussian analysis and CCD (unimodal or bimodal)
        # 4) Gaussian analysis and non-CCD (unimodal or bimodal)
        ccdunifolder = "figures/CCDUnimodal"
        ccdunibifolder = "figures/CCDUnimodalAndBimodal"
        ccdpbifolder = "figures/CCDBimodal"
        ccdgaussccdfolder = "figures/GaussAndCCD"
        ccdgaussnonccdfolder = "figures/GaussAndNonCCD"
        nongaussccdfolder = "figures/CCDAndNonGauss"
        nonccdfolder = "figures/NonCCD"
        bimodalnonccdfolder = "figures/NonCCDBimodal"
        examplesfolder = "figures/Examples"
        folderlist = [
            ccdunifolder,
            ccdunibifolder,
            ccdpbifolder,
            ccdgaussccdfolder,
            ccdgaussnonccdfolder,
            nongaussccdfolder,
            nonccdfolder,
            bimodalnonccdfolder,
            examplesfolder,
        ]
        for f in folderlist:
            if not os.path.exists(f):
                os.mkdir(f)

        # CCD Unimodal
        for ensg in self.fileprefixes[wp_comp_ccd_difffromrng]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdunifolder, ensg + "_mvavg.pdf"),
            )

        # CCD Unimodal and Bimodal
        for ensg in self.fileprefixes[
            wp_comp_ccd_difffromrng & (wp_comp_ccd_clust1 | wp_comp_ccd_clust2)
        ]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_mvavg.pdf"),
            )
        for ensg in self.fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_clust1_mvavg.pdf"),
            )
        for ensg in self.fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust2_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_clust2_mvavg.pdf"),
            )

        # CCD Bimodal
        for ensg in self.fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust1_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1_mvavg.pdf"),
            )
        for ensg in self.fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust2_mvavg.pdf"),
            )

        # Non-CCD
        for ensg in self.fileprefixes[
            ~wp_comp_ccd_difffromrng & ~wp_comp_ccd_clust1 & ~wp_comp_ccd_clust2
        ]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(nonccdfolder, ensg + "_mvavg.pdf"),
            )

        # Non-CCD Bimodal
        for ensg in self.fileprefixes[
            wp_isbimodal_fcpadj_pass & ~wp_comp_ccd_clust1 & ~wp_comp_ccd_clust2
        ]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(bimodalnonccdfolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust1_mvavg.pdf"),
                os.path.join(bimodalnonccdfolder, ensg + "_clust1_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust2_mvavg.pdf"),
                os.path.join(bimodalnonccdfolder, ensg + "_clust2_mvavg.pdf"),
            )

        # Gauss and CCD
        for ensg in self.fileprefixes[wp_comp_ccd_gauss & wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdgaussccdfolder, ensg + "_mvavg.pdf"),
            )
        # Gauss and Non-CCD
        for ensg in self.fileprefixes[wp_comp_ccd_gauss & ~wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdgaussnonccdfolder, ensg + "_mvavg.pdf"),
            )
        # Non-Gauss and CCD
        for ensg in self.fileprefixes[~wp_comp_ccd_gauss & wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(nongaussccdfolder, ensg + "_mvavg.pdf"),
            )
