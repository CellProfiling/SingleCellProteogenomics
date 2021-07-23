# -*- coding: utf-8 -*-
"""
Methods for assessing bimodality of protein intensity distributions.

Distinct high- and low-expressing cell populations that had no correlation to cell division time
were evaluated separately for cell cycle dependence in subsequent analysis.

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

from SingleCellProteogenomics import utils, ProteinBimodality, MovingAverages
import sklearn.mixture
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy


class ProteinBimodality:
    def __init__(self):
        pass

    def identify_bimodal_intensity_distributions(self, protein_ccd):
        """
        Some proteins display bimodal intensity distributions.
        This method seeks to identify distributions with high- and low-expressing cells,
            so that they may be assessed for CCD independently in `ProteinCellCycleDependence.py`.
        """
        u_well_plates = protein_ccd.u_well_plates
        wp_iscell = protein_ccd.wp_iscell
        wp_isnuc = protein_ccd.wp_isnuc
        pol_sort_well_plate = protein_ccd.fucci_polar_coords.pol_sort_well_plate
        pol_sort_norm_rev = protein_ccd.fucci_polar_coords.pol_sort_norm_rev
        pol_sort_ab_cell = protein_ccd.fucci_polar_coords.pol_sort_ab_cell
        pol_sort_ab_nuc = protein_ccd.fucci_polar_coords.pol_sort_ab_nuc
        pol_sort_ab_cyto = protein_ccd.fucci_polar_coords.pol_sort_ab_cyto

        wp_bimodal_cluster_idxs = []
        wp_bimodal_diffmeans = []
        wp_bimodal_fcmeans = []
        wp_bimodal_fcmaxmin = []
        wp_bimodal_clusterlabels = []
        wp_isbimodal_p = []
        wp_timebimodal_p = []
        wp_intensities = []

        # Use Gaussian clustering to investigate if there is bimodality
        gaussian = sklearn.mixture.GaussianMixture(
            n_components=2, random_state=1, max_iter=500
        )
        for i, well in enumerate(u_well_plates):
            curr_well_inds = pol_sort_well_plate == well
            # the reversal isn't really helpful here
            curr_pol = pol_sort_norm_rev[curr_well_inds]
            curr_ab_cell = pol_sort_ab_cell[curr_well_inds]
            curr_ab_nuc = pol_sort_ab_nuc[curr_well_inds]
            curr_ab_cyto = pol_sort_ab_cyto[curr_well_inds]

            # Normalize mean intensities, normalized for display
            curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell)
            curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
            curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto)
            curr_comp_norm = np.asarray(
                curr_ab_cell_norm
                if wp_iscell[i]
                else curr_ab_nuc_norm
                if wp_isnuc[i]
                else curr_ab_cyto_norm
            )
            wp_intensities.append(curr_comp_norm)

            cluster_labels = gaussian.fit_predict(curr_comp_norm.reshape(1, -1).T)
            wp_bimodal_clusterlabels.append(cluster_labels)
            c1 = cluster_labels == 0
            c2 = cluster_labels == 1
            wp_bimodal_cluster_idxs.append([c1, c2])
            wp_bimodal_diffmeans.append(
                np.mean(curr_comp_norm[c2]) - np.mean(curr_comp_norm[c1])
            )
            wp_bimodal_fcmeans.append(
                np.mean(curr_comp_norm[c2]) / np.mean(curr_comp_norm[c1])
            )
            wp_bimodal_fcmaxmin.append(np.max(curr_comp_norm) / np.min(curr_comp_norm))

            # Use a kruskal-wallis test to assess whether there's a significant difference of intensities between clusters
            k, p = scipy.stats.kruskal(curr_comp_norm[c1], curr_comp_norm[c2])
            wp_isbimodal_p.append(p)

            # Use a kruskal-wallis test to assess whether there's (not) a sigificant difference in pseudotime between clusters,
            # since strongly CCD proteins will produce bimodal intensity distributions that should still be assessed as one population
            k, p = scipy.stats.kruskal(curr_pol[c1], curr_pol[c2])
            wp_timebimodal_p.append(p)

        # Multiple testing corrections
        wp_isbimodal_padj, wp_isbimodal_pass = utils.benji_hoch(0.01, wp_isbimodal_p)
        wp_timebimodal_padj, wp_timebimodal_pass = utils.benji_hoch(
            0.01, wp_timebimodal_p
        )

        # Criteria for calling bimodality
        wp_enoughcellsinbothclusters = np.array(
            [sum(c1[0]) > 50 and sum(c1[1]) > 50 for c1 in wp_bimodal_cluster_idxs]
        )
        wp_isbimodal_generally = (
            np.abs(np.log(wp_bimodal_fcmeans) / np.log(2)) > 1
        ) & wp_isbimodal_pass
        wp_isbimodal_fcpadj_pass = (
            (np.abs(np.log(wp_bimodal_fcmeans) / np.log(2)) > 1)
            & wp_isbimodal_pass
            & ~wp_timebimodal_pass
            & wp_enoughcellsinbothclusters
        )

        wp_removeReplicate = np.isin(u_well_plates, protein_ccd.removeThese)
        s_unimodal = f"{sum(~wp_isbimodal_generally[~wp_removeReplicate])}: number of proteins displaying unimodal distributions ({sum(~wp_isbimodal_generally)/len(wp_isbimodal_generally)}%)"
        s_bimodal = f"{sum(wp_isbimodal_generally[~wp_removeReplicate])}: number of proteins displaying bimodal distributions ({sum(wp_isbimodal_generally)/len(wp_isbimodal_generally)}%)"
        print(s_unimodal)
        print(s_bimodal)

        self.wp_isbimodal_fcpadj_pass = wp_isbimodal_fcpadj_pass
        self.wp_bimodal_cluster_idxs = wp_bimodal_cluster_idxs
        self.wp_isbimodal_generally = wp_isbimodal_generally
        self.wp_bimodal_fcmaxmin = wp_bimodal_fcmaxmin

        if protein_ccd.do_plotting:
            # Show that the intensity measurements are reasonable for these bimodal samples
            plt.hist(
                np.concatenate(
                    np.array(wp_intensities, dtype=object)[wp_isbimodal_generally]
                )
            )
            plt.xlabel("Mean intensity")
            plt.ylabel("Count")
            plt.title(
                "Intensities of Cells within Bimodal Distributions\nAre Similar to those Overall"
            )
            plt.close()

            # Illustrate the significantly distinct high- and low-expressing cell populations
            plt.scatter(
                np.log(wp_bimodal_fcmeans) / np.log(2),
                -np.log10(wp_isbimodal_padj),
                c=wp_isbimodal_generally,
                alpha=0.5,
                cmap="bwr_r",
            )
            plt.xlabel("Log2 Fold Change Between Gaussian Clusters")
            plt.ylabel("-Log10 Adj. p-Value for Difference Between Clusters")
            plt.savefig("figures/BimodalSignificance_GeneralBimodality.png")
            plt.close()

            # Illustrate the significantly distinct high- and low-expressing cell populations
            # with no difference in pseudotime. These are evaluated separately for CCD.
            plt.scatter(
                np.log(wp_bimodal_fcmeans) / np.log(2),
                -np.log10(wp_isbimodal_padj),
                c=wp_isbimodal_fcpadj_pass,
                alpha=0.5,
                cmap="bwr_r",
            )
            plt.xlabel("Log2 Fold Change Between Gaussian Clusters")
            plt.ylabel("-Log10 Adj. p-Value for Difference Between Clusters")
            plt.savefig("figures/BimodalSignificance.png")
            plt.savefig("figures/BimodalSignificance.pdf")
            plt.close()

            # Illustrate the samples with sufficient cell count for CCD evaluation of high- and low-expressing cell populations.
            plt.scatter(
                [sum(c1[0]) for c1 in wp_bimodal_cluster_idxs],
                [sum(c1[1]) for c1 in wp_bimodal_cluster_idxs],
                c=wp_enoughcellsinbothclusters,
                alpha=0.5,
                cmap="bwr_r",
            )
            plt.xlabel("Cell Count, Cluster 1")
            plt.ylabel("Cell Count, Cluster 2")
            plt.savefig("figures/BimodalCellCount.png")
            plt.close()
            
        def assess_bimodal_ccd(self, protein_ccd, i, curr_comp, fileprefixes, mvavg_window):
            '''Assess CCD in high- and low-expressing cell populations of proteins determined to have bimodal population intensity'''
            clusters = None
            if self.wp_isbimodal_fcpadj_pass[i]:
                (
                    clust1_idx,
                    clust2_idx,
                ) = self.wp_bimodal_cluster_idxs[i]
                if protein_ccd.do_remove_outliers:
                    clust1_idx, clust2_idx = MovingAverages.remove_outliers(
                        curr_comp, clust1_idx
                    ), MovingAverages.remove_outliers(curr_comp, clust2_idx)
                perc_var_comp_clust1_val, mvavg_clust1 = MovingAverages.mvavg_perc_var(
                    curr_comp_norm[clust1_idx], mvavg_window
                )
                perc_var_comp_clust2_val, mvavg_clust2 = MovingAverages.mvavg_perc_var(
                    curr_comp_norm[clust2_idx], mvavg_window
                )
                self.mvavgs_x_clust1.append(
                    MovingAverages.mvavg(curr_pol[clust1_idx], mvavg_window)
                )
                self.mvavgs_x_clust2.append(
                    MovingAverages.mvavg(curr_pol[clust2_idx], mvavg_window)
                )
                self.perc_var_comp_clust1.append(perc_var_comp_clust1_val)
                self.perc_var_comp_clust2.append(perc_var_comp_clust2_val)
                self.mvavgs_comp_clust1.append(mvavg_clust1)
                self.mvavgs_comp_clust2.append(mvavg_clust2)

                clust1gt = np.mean(mvavg_clust1) > np.mean(mvavg_clust2)
                clusters = np.array([clust1gt] * len(clust1_idx))
                clusters[clust2_idx] = not clust1gt

                perms1 = np.asarray(
                    [
                        np.random.permutation(sum(clust1_idx))
                        for nnn in np.arange(PERMUTATIONS)
                    ]
                )
                perms2 = np.asarray(
                    [
                        np.random.permutation(sum(clust2_idx))
                        for nnn in np.arange(PERMUTATIONS)
                    ]
                )
                curr_comp_perm1 = np.asarray(
                    [curr_comp_norm[clust1_idx][perm] for perm in perms1]
                )
                curr_comp_perm2 = np.asarray(
                    [curr_comp_norm[clust2_idx][perm] for perm in perms2]
                )
                curr_clust1_mvavg_rng_comp = np.apply_along_axis(
                    MovingAverages.mvavg, 1, curr_comp_perm1, mvavg_window
                )
                curr_clust2_mvavg_rng_comp = np.apply_along_axis(
                    MovingAverages.mvavg, 1, curr_comp_perm2, mvavg_window
                )
                curr_clust1_percvar_rng_comp = np.var(
                    curr_clust1_mvavg_rng_comp, axis=1
                ) / np.var(curr_comp_perm1, axis=1)
                curr_clust2_percvar_rng_comp = np.var(
                    curr_clust2_mvavg_rng_comp, axis=1
                ) / np.var(curr_comp_perm2, axis=1)
                self.perc_var_comp_clust1_rng.append(curr_clust1_percvar_rng_comp)
                self.perc_var_comp_clust2_rng.append(curr_clust2_percvar_rng_comp)

                windows1 = np.asarray(
                    [
                        np.arange(start, start + mvavg_window)
                        for start in np.arange(sum(clust1_idx) - mvavg_window + 1)
                    ]
                )
                windows2 = np.asarray(
                    [
                        np.arange(start, start + mvavg_window)
                        for start in np.arange(sum(clust2_idx) - mvavg_window + 1)
                    ]
                )
                mvperc1 = MovingAverages.mvpercentiles(
                    curr_comp_norm[clust1_idx][windows1]
                )
                mvperc2 = MovingAverages.mvpercentiles(
                    curr_comp_norm[clust2_idx][windows2]
                )
                if self.do_plotting:
                    MovingAverages.temporal_mov_avg_protein(
                        curr_pol[clust1_idx],
                        curr_comp_norm[clust1_idx],
                        self.mvavgs_x_clust1[-1],
                        self.mvavgs_comp_clust1[-1],
                        mvperc1,
                        None,
                        self.folder,
                        fileprefixes[i] + "_clust1",
                    )
                    MovingAverages.temporal_mov_avg_protein(
                        curr_pol[clust2_idx],
                        curr_comp_norm[clust2_idx],
                        self.mvavgs_x_clust2[-1],
                        self.mvavgs_comp_clust2[-1],
                        mvperc2,
                        None,
                        self.folder,
                        fileprefixes[i] + "_clust2",
                    )
                return clusters
