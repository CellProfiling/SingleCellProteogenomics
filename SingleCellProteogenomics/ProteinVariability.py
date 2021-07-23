# -*- coding: utf-8 -*-
"""
The variance of protein expression in the cell populations was evaluated using:
    - variance
    - coefficient of variation (CV)
    - Gini index

@author: Anthony J. Cesnik, cesnik@stanford.edu
@author: devinsullivan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from SingleCellProteogenomics import utils

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300


class ProteinVariability:
    """Class for investigating protein variability"""

    def __init__(self):
        pass

    def plot_average_intensities_by_batch(self):
        """Plot average intensities by batch to look at agreement"""
        firstbatch = np.asarray(
            [not str(p).split("_")[1].startswith("67") for p in self.u_well_plates]
        )
        plt.hist(
            np.array(self.mean_mean_comp)[firstbatch], bins=200, label="firstbatch"
        )
        plt.hist(
            np.array(self.mean_mean_comp)[~firstbatch], bins=200, label="secondbatch"
        )
        plt.legend()
        plt.savefig("figures/MeanMeancComp.png")
        plt.close()

        firstbatch = np.asarray(
            [not str(p).split("_")[1].startswith("67") for p in self.u_well_plates]
        )
        plt.hist(np.array(self.mean_mean_mt)[firstbatch], bins=200, label="firstbatch")
        plt.hist(
            np.array(self.mean_mean_mt)[~firstbatch], bins=200, label="secondbatch"
        )
        plt.legend()
        plt.savefig("figures/MeanMeanMt.png")
        plt.close()

    def variability_plots(self):
        """Make boxplots and scatterplots comparing the different types of variance"""
        # Boxplots
        utils.general_boxplot(
            (self.var_cell, self.var_cyto, self.var_nuc, self.var_mt),
            ("var_cell", "var_cyto", "var_nuc", "var_mt"),
            "Metacompartment",
            f"Variance using {'log' if self.use_log else 'natural'} intensity values",
            "",
            True,
            "figures/VarianceBoxplot.png",
        )
        utils.general_boxplot(
            (self.cv_cell, self.cv_cyto, self.cv_nuc, self.cv_mt),
            ("cv_cell", "cv_cyto", "cv_nuc", "cv_mt"),
            "Metacompartment",
            f"Coeff. of Var. using {'log' if self.use_log else 'natural'} intensity values",
            "",
            True,
            "figures/CVBoxplot.png",
        )
        utils.general_boxplot(
            (self.gini_cell, self.gini_cyto, self.gini_nuc, self.gini_mt),
            ("gini_cell", "gini_cyto", "gini_nuc", "gini_mt"),
            "Metacompartment",
            f"Gini using {'log' if self.use_log else 'natural'} intensity values",
            "",
            True,
            "figures/GiniBoxplot.png",
        )

        # Scatters
        utils.general_scatter(
            self.var_comp, self.var_mt, "var_comp", "var_mt", "figures/var_comp_mt.png"
        )
        utils.general_scatter(
            self.cv_comp, self.cv_mt, "cv_comp", "cv_mt", "figures/cv_comp_mt.png"
        )
        utils.general_scatter(
            self.gini_comp,
            self.gini_mt,
            "gini_comp",
            "gini_mt",
            "figures/gini_comp_mt.png",
        )
        utils.general_scatter(
            self.var_comp,
            self.mean_mean_comp,
            "var_comp",
            f"{'log10' if self.use_log else 'natural'} intensity",
            "figures/VarianceVsIntensityComp.png",
        )

    def compare_image_to_sample_level_variances(self):
        '''Compare the variances of image-level and sample-level variability'''
        s_maxvariance = (
            np.concatenate(self.wpi_img)[np.argmax(np.concatenate(self.var_comp_img))]
            + ": the image with the max variance"
        )
        print(s_maxvariance)
        high_var_img = np.concatenate(self.wpi_img)[
            np.concatenate(
                [vvv > 4 * self.var_comp[i] for i, vvv in enumerate(self.var_comp_img)]
            )
        ]
        s_highvar = f"{high_var_img}: the images with greater than 4x the variance of the whole sample"
        print(s_highvar)
        
        norm_cv_img = np.concatenate(
            [vvv / self.cv_comp[i] for i, vvv in enumerate(self.cv_comp_img)]
        )
        cutoff = np.mean(norm_cv_img) + 3 * np.std(norm_cv_img)
        high_cv_img = np.concatenate(self.wpi_img)[norm_cv_img > cutoff]
        s_highvar = f"{high_cv_img}: the images with greater than 4x the variance of the whole sample"
        print(s_highvar)
        # np.intersect1d(high_var_img, high_cv_img)

        if self.do_plotting:
            utils.general_scatter(
                np.concatenate(
                    [
                        [self.var_comp[i]] * len(vvv)
                        for i, vvv in enumerate(self.var_comp_img)
                    ]
                ),
                np.concatenate(self.var_comp_img),
                "variance within compartment",
                "variance for each image",
                "figures/VarianceByImage.png",
            )
            utils.general_scatter(
                np.concatenate(
                    [
                        [self.gini_comp[i]] * len(vvv)
                        for i, vvv in enumerate(self.gini_comp_img)
                    ]
                ),
                np.concatenate(self.gini_comp_img),
                "gini within compartment",
                "gini for each image",
                "figures/GiniByImage.png",
            )
            utils.general_scatter(
                np.concatenate(
                    [
                        [self.cv_comp[i]] * len(vvv)
                        for i, vvv in enumerate(self.cv_comp_img)
                    ]
                ),
                np.concatenate(self.cv_comp_img),
                "cv within compartment",
                "cv for each image",
                "figures/CVByImage.png",
            )
            # plt.hist(
            #     np.concatenate(
            #         [vvv / self.var_comp[i] for i, vvv in enumerate(self.var_comp_img)]
            #     )
            # )
            # plt.close()
            # plt.hist(norm_cv_img)
            # plt.close()

    def compare_variances(self):
        """Compare variances for protein and microtubules, the internal control for each image"""
        # p-value (one-sided) calculations for comparing variances to microtubule intensity variances
        p_varProt_varMt = 2 * scipy.stats.kruskal(
            self.var_comp[~np.isin(self.u_well_plates, self.removeThese)],
            self.var_mt[~np.isin(self.u_well_plates, self.removeThese)],
        )[1]
        p_cvProt_cvMt = 2 * scipy.stats.kruskal(
            self.cv_comp[~np.isin(self.u_well_plates, self.removeThese)],
            self.cv_mt[~np.isin(self.u_well_plates, self.removeThese)],
        )[1]
        p_giniProt_giniMt = 2 * scipy.stats.kruskal(
            self.gini_comp[~np.isin(self.u_well_plates, self.removeThese)],
            self.gini_mt[~np.isin(self.u_well_plates, self.removeThese)],
        )[1]
        p_statements = [
            f"{p_varProt_varMt}: one-sided p-value for difference between protein and microtubule variances",
            f"{p_cvProt_cvMt}: one-sided p-value for difference between protein and microtubule CVs",
            f"{p_giniProt_giniMt}: one-sided p-value for difference between protein and microtubule Gini indices",
        ]
        print("\n".join(p_statements))

        if self.do_plotting:
            # Boxplots comparing variances to microtubule variances
            # Remove duplicated antibodies to make these independent samples for one-sided Kruskal-Wallis tests
            utils.general_boxplot(
                (
                    self.var_comp[~np.isin(self.u_well_plates, self.removeThese)],
                    self.var_mt[~np.isin(self.u_well_plates, self.removeThese)],
                ),
                ("Protein", "Microtubules"),
                "",
                "Variance",
                "",
                False,
                "figures/ProteinMicrotubuleVariances.pdf",
            )
            utils.general_boxplot(
                (
                    self.cv_comp[~np.isin(self.u_well_plates, self.removeThese)],
                    self.gini_mt[~np.isin(self.u_well_plates, self.removeThese)],
                ),
                ("Protein", "Microtubules"),
                "",
                "CV",
                "",
                False,
                "figures/ProteinMicrotubuleCVs.pdf",
            )
            utils.general_boxplot(
                (
                    self.gini_comp[~np.isin(self.u_well_plates, self.removeThese)],
                    self.gini_mt[~np.isin(self.u_well_plates, self.removeThese)],
                ),
                ("Protein", "Microtubules"),
                "",
                "Gini",
                "",
                False,
                "figures/ProteinMicrotubuleGinis.pdf",
            )

    def calculate_variation(self, protein_ccd, use_log=False):
        """Calculate overall variation of protein staining intensity in single cells"""
        self.do_plotting = protein_ccd.do_plotting
        self.use_log = use_log
        self.removeThese = protein_ccd.removeThese
        self.u_well_plates = protein_ccd.u_well_plates
        self.wp_iscell = protein_ccd.wp_iscell
        self.wp_isnuc = protein_ccd.wp_isnuc
        self.wp_iscyto = protein_ccd.wp_iscyto
        self.pol_sort_well_plate = protein_ccd.fucci_polar_coords.pol_sort_well_plate
        self.pol_sort_ab_cell = protein_ccd.fucci_polar_coords.pol_sort_ab_cell
        self.pol_sort_ab_nuc = protein_ccd.fucci_polar_coords.pol_sort_ab_nuc
        self.pol_sort_ab_cyto = protein_ccd.fucci_polar_coords.pol_sort_ab_cyto
        self.pol_sort_mt_cell = protein_ccd.fucci_polar_coords.pol_sort_mt_cell
        self.pol_sort_well_plate_imgnb = (
            protein_ccd.fucci_polar_coords.pol_sort_well_plate_imgnb
        )

        # mean intensity variances per antibody
        self.var_cell = []
        self.var_nuc = []
        self.var_cyto = []
        self.var_mt = []

        # mean intensity CVs per antibody
        self.cv_cell = []
        self.cv_nuc = []
        self.cv_cyto = []
        self.cv_mt = []

        # mean intensity ginis per antibody
        self.gini_cell = []
        self.gini_nuc = []
        self.gini_cyto = []
        self.gini_mt = []

        # mean mean-intensity
        self.mean_mean_cell = []
        self.mean_mean_nuc = []
        self.mean_mean_cyto = []
        self.mean_mean_mt = []

        self.cell_counts = []

        self.wpi_img = []

        # mean intensity gini per field of view
        self.gini_cell_img = []
        self.gini_nuc_img = []
        self.gini_cyto_img = []
        self.gini_mt_img = []

        # mean intensity variances per field of view
        self.var_cell_img = []
        self.var_nuc_img = []
        self.var_cyto_img = []
        self.var_mt_img = []

        # mean intensity CVs per field of view
        self.cv_cell_img = []
        self.cv_nuc_img = []
        self.cv_cyto_img = []
        self.cv_mt_img = []

        # The variance needs to be calculated separately for each well because they all have different numbers of cells
        for well in self.u_well_plates:
            curr_well_inds = self.pol_sort_well_plate == well
            curr_ab_cell = (
                self.pol_sort_ab_cell[curr_well_inds]
                if not use_log
                else np.log10(self.pol_sort_ab_cell[curr_well_inds])
            )
            curr_ab_nuc = (
                self.pol_sort_ab_nuc[curr_well_inds]
                if not use_log
                else np.log10(self.pol_sort_ab_nuc[curr_well_inds])
            )
            curr_ab_cyto = (
                self.pol_sort_ab_cyto[curr_well_inds]
                if not use_log
                else np.log10(self.pol_sort_ab_cyto[curr_well_inds])
            )
            curr_mt_cell = (
                self.pol_sort_mt_cell[curr_well_inds]
                if not use_log
                else np.log10(self.pol_sort_mt_cell[curr_well_inds])
            )

            self.cell_counts.append(len(curr_ab_cell))

            self.var_cell.append(np.var(curr_ab_cell))
            self.var_nuc.append(np.var(curr_ab_nuc))
            self.var_cyto.append(np.var(curr_ab_cyto))
            self.var_mt.append(np.var(curr_mt_cell))

            self.cv_cell.append(scipy.stats.variation(curr_ab_cell))
            self.cv_nuc.append(scipy.stats.variation(curr_ab_nuc))
            self.cv_cyto.append(scipy.stats.variation(curr_ab_cyto))
            self.cv_mt.append(scipy.stats.variation(curr_mt_cell))

            self.gini_cell.append(utils.gini(curr_ab_cell))
            self.gini_nuc.append(utils.gini(curr_ab_nuc))
            self.gini_cyto.append(utils.gini(curr_ab_cyto))
            self.gini_mt.append(utils.gini(curr_mt_cell))

            # Save the mean mean intensities
            self.mean_mean_cell.append(np.mean(curr_ab_cell))
            self.mean_mean_nuc.append(np.mean(curr_ab_nuc))
            self.mean_mean_cyto.append(np.mean(curr_ab_cyto))
            self.mean_mean_mt.append(np.mean(curr_mt_cell))

            curr_well_plate_imgnbs = self.pol_sort_well_plate_imgnb[curr_well_inds]
            curr_wpi_img = []

            # mean intensity ginis per field of view
            curr_gini_cell_img = []
            curr_gini_nuc_img = []
            curr_gini_cyto_img = []
            curr_gini_mt_img = []

            # mean intensity variances per field of view
            curr_var_cell_img = []
            curr_var_nuc_img = []
            curr_var_cyto_img = []
            curr_var_mt_img = []

            # mean intensity CVs per field of view
            curr_cv_cell_img = []
            curr_cv_nuc_img = []
            curr_cv_cyto_img = []
            curr_cv_mt_img = []

            for wpi in np.unique(curr_well_plate_imgnbs):
                curr_wpis = self.pol_sort_well_plate_imgnb == wpi
                curr_ab_cell = (
                    self.pol_sort_ab_cell[curr_wpis]
                    if not use_log
                    else np.log10(self.pol_sort_ab_cell[curr_wpis])
                )
                curr_ab_nuc = (
                    self.pol_sort_ab_nuc[curr_wpis]
                    if not use_log
                    else np.log10(self.pol_sort_ab_nuc[curr_wpis])
                )
                curr_ab_cyto = (
                    self.pol_sort_ab_cyto[curr_wpis]
                    if not use_log
                    else np.log10(self.pol_sort_ab_cyto[curr_wpis])
                )
                curr_mt_cell = (
                    self.pol_sort_mt_cell[curr_wpis]
                    if not use_log
                    else np.log10(self.pol_sort_mt_cell[curr_wpis])
                )

                curr_wpi_img.append(wpi)

                curr_var_cell_img.append(np.var(curr_ab_cell))
                curr_var_nuc_img.append(np.var(curr_ab_nuc))
                curr_var_cyto_img.append(np.var(curr_ab_cyto))
                curr_var_mt_img.append(np.var(curr_mt_cell))

                curr_gini_cell_img.append(utils.gini(curr_ab_cell))
                curr_gini_nuc_img.append(utils.gini(curr_ab_nuc))
                curr_gini_cyto_img.append(utils.gini(curr_ab_cyto))
                curr_gini_mt_img.append(utils.gini(curr_mt_cell))

                curr_cv_cell_img.append(scipy.stats.variation(curr_ab_cell))
                curr_cv_nuc_img.append(scipy.stats.variation(curr_ab_nuc))
                curr_cv_cyto_img.append(scipy.stats.variation(curr_ab_cyto))
                curr_cv_mt_img.append(scipy.stats.variation(curr_mt_cell))

            self.wpi_img.append(curr_wpi_img)
            self.var_cell_img.append(curr_var_cell_img)
            self.var_nuc_img.append(curr_var_nuc_img)
            self.var_cyto_img.append(curr_var_cyto_img)
            self.var_mt_img.append(curr_var_mt_img)

            self.gini_cell_img.append(curr_gini_cell_img)
            self.gini_nuc_img.append(curr_gini_nuc_img)
            self.gini_cyto_img.append(curr_gini_cyto_img)
            self.gini_mt_img.append(curr_gini_mt_img)

            self.cv_cell_img.append(curr_cv_cell_img)
            self.cv_nuc_img.append(curr_cv_nuc_img)
            self.cv_cyto_img.append(curr_cv_cyto_img)
            self.cv_mt_img.append(curr_cv_mt_img)

        # Make numpy arrays and store
        self.var_cell, self.var_nuc, self.var_cyto, self.var_mt = (
            np.array(self.var_cell),
            np.array(self.var_nuc),
            np.array(self.var_cyto),
            np.array(self.var_mt),
        )
        self.gini_cell, self.gini_nuc, self.gini_cyto, self.gini_mt = (
            np.array(self.gini_cell),
            np.array(self.gini_nuc),
            np.array(self.gini_cyto),
            np.array(self.gini_mt),
        )
        self.cv_cell, self.cv_nuc, self.cv_cyto, self.cv_mt = (
            np.array(self.cv_cell),
            np.array(self.cv_nuc),
            np.array(self.cv_cyto),
            np.array(self.cv_mt),
        )
        self.var_comp_img = utils.values_comp(
            self.var_cell_img,
            self.var_nuc_img,
            self.var_cyto_img,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.gini_comp_img = utils.values_comp(
            self.gini_cell_img,
            self.gini_nuc_img,
            self.gini_cyto_img,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.cv_comp_img = utils.values_comp(
            self.cv_cell_img,
            self.cv_nuc_img,
            self.cv_cyto_img,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.mean_mean_comp = utils.values_comp(
            self.mean_mean_cell,
            self.mean_mean_nuc,
            self.mean_mean_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.cv_comp = utils.values_comp(
            self.cv_cell,
            self.cv_nuc,
            self.cv_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.gini_comp = utils.values_comp(
            self.gini_cell,
            self.gini_nuc,
            self.gini_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
        self.var_comp = utils.values_comp(
            self.var_cell,
            self.var_nuc,
            self.var_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )

        if self.do_plotting:
            print("Plotting average intensities of proteins and microtubules by batch.")
            self.plot_average_intensities_by_batch()

            print("Making general plots for variance, CV, and gini by compartment")
            self.variability_plots()

        print("Comparing image to sample variance")
        self.compare_image_to_sample_level_variances()
        self.compare_variances()

        # Pickle and return main results
        utils.np_save_overwriting(
            "output/pickles/mean_mean_comp.npy", self.mean_mean_comp
        )
        utils.np_save_overwriting("output/pickles/cv_comp.npy", self.cv_comp)
        utils.np_save_overwriting("output/pickles/gini_comp.npy", self.gini_comp)
        utils.np_save_overwriting("output/pickles/var_comp.npy", self.var_comp)
        utils.np_save_overwriting("output/pickles/cv_cell.npy", self.cv_cell)
        utils.np_save_overwriting("output/pickles/gini_cell.npy", self.gini_cell)
        utils.np_save_overwriting("output/pickles/var_cell.npy", self.var_cell)
