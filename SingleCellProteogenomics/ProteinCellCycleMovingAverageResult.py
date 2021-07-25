#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 12:22:17 2021

@author: anthony.cesnik
"""

from SingleCellProteogenomics import MovingAverages
import numpy as np
import pandas as pd

# Get the same results each time
np.random.seed(0)

# Number of points for moving average window for protein analysis
WINDOW = 10

# Number of permutations used for randomization analysis
PERMUTATIONS = 10000


class ProteinCellCycleMovingAverageResult:
    def __init__(self, protein_ccd, idx):
        self.idx = idx
        self.protein_ccd = protein_ccd
        self.use_log_ccd = protein_ccd.use_log_ccd
        self.curr_well_inds = protein_ccd.curr_well_inds
        self.fucci_polar_coords = protein_ccd.fucci_polar_coords
        self.do_remove_outliers = protein_ccd.do_remove_outliers
        self.wp_iscell = protein_ccd.wp_iscell
        self.wp_isnuc = protein_ccd.wp_isnuc
        self.wp_iscyto = protein_ccd.wp_iscyto

    def vals_at_inds(self, source, curr_well_inds):
        if self.use_log_ccd:
            return np.log10(source[curr_well_inds])
        else:
            return source[curr_well_inds]

    def get_curr_comp_vals(self, cell, nuc, cyto):
        return (
            cell
            if self.wp_iscell[self.idx]
            else nuc
            if self.wp_isnuc[self.idx]
            else cyto
        )

    def calculate_moving_averages(self):
        curr_well_inds = self.curr_well_inds[self.idx]
        self.curr_pol = self.fucci_polar_coords.pol_sort_norm_rev[curr_well_inds]
        curr_ab_cell = self.vals_at_inds(
            self.fucci_polar_coords.pol_sort_ab_cell, curr_well_inds
        )
        curr_ab_nuc = self.vals_at_inds(
            self.fucci_polar_coords.pol_sort_ab_nuc, curr_well_inds
        )
        curr_ab_cyto = self.vals_at_inds(
            self.fucci_polar_coords.pol_sort_ab_cyto, curr_well_inds
        )
        curr_mt_cell = self.vals_at_inds(
            self.fucci_polar_coords.pol_sort_mt_cell, curr_well_inds
        )
        self.curr_fred = self.fucci_polar_coords.pol_sort_fred[curr_well_inds]
        self.curr_fgreen = self.fucci_polar_coords.pol_sort_fgreen[curr_well_inds]
        self.curr_mockbulk_phase = self.fucci_polar_coords.pol_sort_mockbulk_phases[
            curr_well_inds
        ]
        if self.do_remove_outliers:
            curr_comp = self.get_curr_comp_vals(curr_ab_cell, curr_ab_nuc, curr_ab_cyto)
            self.curr_pol = MovingAverages.remove_outliers(curr_comp, self.curr_pol)
            curr_ab_cell = MovingAverages.remove_outliers(curr_comp, curr_ab_cell)
            curr_ab_nuc = MovingAverages.remove_outliers(curr_comp, curr_ab_nuc)
            curr_ab_cyto = MovingAverages.remove_outliers(curr_comp, curr_ab_cyto)
            curr_mt_cell = MovingAverages.remove_outliers(curr_comp, curr_mt_cell)
            self.curr_fred = MovingAverages.remove_outliers(curr_comp, self.curr_fred)
            self.curr_fgreen = MovingAverages.remove_outliers(
                curr_comp, self.curr_fgreen
            )
            self.curr_mockbulk_phase = MovingAverages.remove_outliers(
                curr_comp, self.curr_mockbulk_phase
            )

        # Normalize mean intensities, normalized for display
        curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell)
        curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
        curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto)
        self.curr_mt_cell_norm = curr_mt_cell / np.max(curr_mt_cell)

        # Original method from Devin's work
        self.perc_var_cell_val, self.mvavg_cell = MovingAverages.mvavg_perc_var(
            curr_ab_cell_norm, WINDOW
        )
        self.perc_var_nuc_val, self.mvavg_nuc = MovingAverages.mvavg_perc_var(
            curr_ab_nuc_norm, WINDOW
        )
        self.perc_var_cyto_val, self.mvavg_cyto = MovingAverages.mvavg_perc_var(
            curr_ab_cyto_norm, WINDOW
        )
        self.perc_var_mt_val, self.mvavg_mt = MovingAverages.mvavg_perc_var(
            self.curr_mt_cell_norm, WINDOW
        )

        self.mvavg_xvals = MovingAverages.mvavg(self.curr_pol, WINDOW)

        self.curr_comp_norm = self.get_curr_comp_vals(
            curr_ab_cell_norm, curr_ab_nuc_norm, curr_ab_cyto_norm
        )
        self.curr_mvavg_comp = self.get_curr_comp_vals(
            self.mvavg_cell, self.mvavg_nuc, self.mvavg_cyto
        )
        self.curr_comp_percvar = self.get_curr_comp_vals(
            self.perc_var_cell_val, self.perc_var_nuc_val, self.perc_var_cyto_val
        )

        # Perform permutation analysis
        self.permutation_analysis_protein()

        # Assess CCD in high- and low-expressing cell populations of proteins determined to have bimodal population intensity
        self.clusters = self.assess_bimodal_ccd(curr_comp)
        self.plotting_bimodals_and_microtubules()

    def permutation_analysis_protein(self):
        """Randomization analysis of cell cycle dependence: permute cell order and calculate percent variance due to the cell cycle"""
        perms = np.asarray(
            [
                np.random.permutation(len(self.curr_pol))
                for nnn in np.arange(PERMUTATIONS)
            ]
        )
        curr_comp_perm = np.asarray([self.curr_comp_norm[perm] for perm in perms])
        curr_mt_perm = np.asarray([self.curr_mt_cell_norm[perm] for perm in perms])
        curr_mvavg_rng_comp = np.apply_along_axis(
            MovingAverages.mvavg, 1, curr_comp_perm, WINDOW
        )
        curr_mvavg_rng_mt = np.apply_along_axis(
            MovingAverages.mvavg, 1, curr_mt_perm, WINDOW
        )
        self.curr_percvar_rng_comp = np.var(curr_mvavg_rng_comp, axis=1) / np.var(
            curr_comp_perm, axis=1
        )
        self.curr_percvar_rng_mt = np.var(curr_mvavg_rng_mt, axis=1) / np.var(
            curr_mt_perm, axis=1
        )
        self.permutation_example_plot(curr_comp_perm, curr_mvavg_rng_comp)

    def bimodal_permutation_analysis(self, clust1_idx, clust2_idx):
        """Randomization analysis for bimodal samples"""
        perms1 = np.asarray(
            [np.random.permutation(sum(clust1_idx)) for nnn in np.arange(PERMUTATIONS)]
        )
        perms2 = np.asarray(
            [np.random.permutation(sum(clust2_idx)) for nnn in np.arange(PERMUTATIONS)]
        )
        curr_comp_perm1 = np.asarray(
            [self.curr_comp_norm[clust1_idx][perm] for perm in perms1]
        )
        curr_comp_perm2 = np.asarray(
            [self.curr_comp_norm[clust2_idx][perm] for perm in perms2]
        )
        curr_clust1_mvavg_rng_comp = np.apply_along_axis(
            MovingAverages.mvavg, 1, curr_comp_perm1, WINDOW
        )
        curr_clust2_mvavg_rng_comp = np.apply_along_axis(
            MovingAverages.mvavg, 1, curr_comp_perm2, WINDOW
        )
        self.curr_clust1_percvar_rng_comp = np.var(
            curr_clust1_mvavg_rng_comp, axis=1
        ) / np.var(curr_comp_perm1, axis=1)
        self.curr_clust2_percvar_rng_comp = np.var(
            curr_clust2_mvavg_rng_comp, axis=1
        ) / np.var(curr_comp_perm2, axis=1)

    def assess_bimodal_ccd(self, curr_comp):
        """Assess CCD in high- and low-expressing cell populations of proteins determined to have bimodal population intensity"""
        self.clusters = None
        bimodality = self.protein_ccd.protein_bimodality
        if bimodality.wp_isbimodal_fcpadj_pass[self.idx]:
            # Calculate moving averages
            clustidxs = bimodality.wp_bimodal_cluster_idxs[self.idx]
            clust1_idx, clust2_idx = clustidxs
            if self.protein_ccd.do_remove_outliers:
                clust1_idx = MovingAverages.remove_outliers(curr_comp, clust1_idx)
                clust2_idx = MovingAverages.remove_outliers(curr_comp, clust2_idx)

            mvmv1 = MovingAverages.mvavg_perc_var(
                self.curr_comp_norm[clust1_idx], WINDOW
            )
            mvmv2 = MovingAverages.mvavg_perc_var(
                self.curr_comp_norm[clust2_idx], WINDOW
            )
            self.perc_var_comp_clust1_val, self.mvavg_clust1 = mvmv1
            self.perc_var_comp_clust2_val, self.mvavg_clust2 = mvmv2

            self.mvavg_x_clust1 = MovingAverages.mvavg(
                self.curr_pol[clust1_idx], WINDOW
            )
            self.mvavg_x_clust2 = MovingAverages.mvavg(
                self.curr_pol[clust2_idx], WINDOW
            )

            clust1gt = np.mean(self.mvavg_clust1) > np.mean(self.mvavg_clust2)
            self.clusters = np.array([clust1gt] * len(clust1_idx))
            self.clusters[clust2_idx] = not clust1gt

            # Permutation analysis
            self.bimodal_permutation_analysis(clust1_idx, clust2_idx)

            # Get percentiles and plot if desired
            clust1_range = np.arange(sum(clust1_idx) - WINDOW + 1)
            clust2_range = np.arange(sum(clust2_idx) - WINDOW + 1)
            windows1 = np.asarray(
                [np.arange(start, start + WINDOW) for start in clust1_range]
            )
            windows2 = np.asarray(
                [np.arange(start, start + WINDOW) for start in clust2_range]
            )
            mvperc1 = MovingAverages.mvpercentiles(
                self.curr_comp_norm[clust1_idx][windows1]
            )
            mvperc2 = MovingAverages.mvpercentiles(
                self.curr_comp_norm[clust2_idx][windows2]
            )

            if self.protein_ccd.do_plotting:
                MovingAverages.temporal_mov_avg_protein(
                    self.curr_pol[clust1_idx],
                    self.curr_comp_norm[clust1_idx],
                    self.mvavg_x_clust1,
                    self.mvavg_clust1,
                    mvperc1,
                    None,
                    self.protein_ccd.folder,
                    self.protein_ccd.fileprefixes[self.idx] + "_clust1",
                )
                MovingAverages.temporal_mov_avg_protein(
                    self.curr_pol[clust2_idx],
                    self.curr_comp_norm[clust2_idx],
                    self.mvavg_x_clust2,
                    self.mvavg_clust2,
                    mvperc2,
                    None,
                    self.protein_ccd.folder,
                    self.protein_ccd.fileprefixes[self.idx] + "_clust2",
                )

    def permutation_example_plot(self, curr_comp_perm, curr_mvavg_rng_comp):
        """Make example plots for the randomization trials for the NFAT5 example in manuscript"""
        if (
            self.protein_ccd.wp_ensg[self.idx] == "ENSG00000102908"
            and self.protein_ccd.do_plotting
        ):
            median_rng_idx = np.argsort(self.curr_percvar_rng_comp)
            intervals = [
                0,
                len(self.curr_percvar_rng_comp) // 4,
                len(self.curr_percvar_rng_comp) // 2,
                3 * len(self.curr_percvar_rng_comp) // 4,
                len(self.curr_percvar_rng_comp) - 1,
            ]
            for iii, idx in enumerate(intervals):
                MovingAverages.temporal_mov_avg_randomization_example_protein(
                    self.curr_pol,
                    self.curr_comp_norm,
                    curr_comp_perm[median_rng_idx[self.idx]],
                    self.mvavg_xvals,
                    self.mvavg_comp,
                    curr_mvavg_rng_comp[median_rng_idx[self.idx]],
                    self.protein_ccd.folder_rng,
                    f"{self.protein_ccd.fileprefixes[self.idx]}_withrandomization_{iii}",
                )

    def plotting_bimodals_and_microtubules(self):
        """Plotting bimodals and microtubules over pseudotime"""
        window_range = np.arange(len(self.curr_pol) - WINDOW + 1)
        windows = np.asarray(
            [np.arange(start, start + WINDOW) for start in window_range]
        )
        self.mvperc_comp = MovingAverages.mvpercentiles(self.curr_comp_norm[windows])

        if self.protein_ccd.do_plotting:
            # Make the plots for each protein (takes 10 mins for all of them)
            MovingAverages.temporal_mov_avg_protein(
                self.curr_pol,
                self.curr_comp_norm,
                self.mvavg_xvals,
                self.mvavg_comp,
                self.mvperc_comp,
                None,
                self.protein_ccd.folder,
                self.protein_ccd.fileprefixes[self.idx],
            )

            # Make the plots for microtubules (takes 10 mins for all, so just do the arbitrary one in the Fig 1)
            wellplate = self.protein_ccd.u_well_plates[self.idx]
            if wellplate == "C07_55405991":
                mvperc_mt = MovingAverages.mvpercentiles(
                    self.curr_mt_cell_norm[windows]
                )
                MovingAverages.temporal_mov_avg_protein(
                    self.curr_pol,
                    self.curr_mt_cell_norm,
                    self.mvavg_xvals,
                    self.mvavg_mt,
                    mvperc_mt,
                    None,
                    self.protein_ccd.folder_mt,
                    f"{self.protein_ccd.fileprefixes[self.idx]}_mt",
                )

                # keep here in case making all pseudotime plots
                if wellplate == "C07_55405991":
                    pd.DataFrame(
                        {
                            "ENSG": self.protein_ccd.wp_ensg[self.idx],
                            "Antibody": self.protein_ccd.wp_ab[self.idx],
                            "Compartment": "Cell",
                            "CCD": "No",
                            "cell_pseudotime": [
                                ",".join([str(ppp) for ppp in pp])
                                for pp in [self.curr_pol]
                            ],
                            "cell_intensity": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [self.curr_mt_cell_norm]
                            ],
                            "mvavg_x": [
                                ",".join([str(xxx) for xxx in xx])
                                for xx in [self.mvavg_xvals]
                            ],
                            "mvavg_y": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [self.mvavg_mt]
                            ],
                            "mvavgs_10p": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [mvperc_mt[0]]
                            ],
                            "mvavgs_90p": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [mvperc_mt[-1]]
                            ],
                            "mvavgs_25p": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [mvperc_mt[1]]
                            ],
                            "mvavgs_75p": [
                                ",".join([str(yyy) for yyy in yy])
                                for yy in [mvperc_mt[-2]]
                            ],
                            "phase": [
                                ",".join(pp) for pp in [self.curr_mockbulk_phase]
                            ],
                            "WellPlate": wellplate,
                        }
                    ).to_csv("output/mtplottingline.tsv", index=False, sep="\t")

            # Make the plots for each bimodal protein
            if self.clusters is not None:
                MovingAverages.temporal_mov_avg_protein(
                    self.curr_pol,
                    self.curr_comp_norm,
                    self.mvavg_xvals,
                    self.mvavg_comp,
                    self.mvperc_comp,
                    self.clusters,
                    self.folder,
                    self.protein_ccd.fileprefixes[self.idx] + "_clust1&2",
                )
