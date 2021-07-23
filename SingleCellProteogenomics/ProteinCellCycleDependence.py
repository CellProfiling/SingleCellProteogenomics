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
    MovingAverages,
    FucciPseudotime,
    ProteinVariability,
    ProteinBimodality,
)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import shutil, os

# Get the same results each time
np.random.seed(0)

# Number of points for moving average window for protein analysis
WINDOW = 10

# Number of permutations used for randomization analysis
PERMUTATIONS = 10000

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


def permutation_analysis_protein(
    idx,
    curr_pol,
    curr_ab_cell_norm,
    curr_ab_nuc_norm,
    curr_ab_cyto_norm,
    curr_mt_cell_norm,
    perc_var_cell_val,
    perc_var_nuc_val,
    perc_var_cyto_val,
    wp_iscell,
    wp_isnuc,
    wp_iscyto,
    mvavg_cell,
    mvavg_nuc,
    mvavg_cyto,
):
    """Randomization analysis of cell cycle dependence: permute cell order and calculate percent variance due to the cell cycle"""
    perms = np.asarray(
        [np.random.permutation(len(curr_pol)) for nnn in np.arange(PERMUTATIONS)]
    )
    curr_comp_norm = np.asarray(
        curr_ab_cell_norm
        if wp_iscell[idx]
        else curr_ab_nuc_norm
        if wp_isnuc[idx]
        else curr_ab_cyto_norm
    )
    curr_comp_percvar = np.asarray(
        perc_var_cell_val
        if wp_iscell[idx]
        else perc_var_nuc_val
        if wp_isnuc[idx]
        else perc_var_cyto_val
    )
    curr_comp_mvavg = np.asarray(
        mvavg_cell if wp_iscell[idx] else mvavg_nuc if wp_isnuc[idx] else mvavg_cyto
    )
    curr_comp_perm = np.asarray([curr_comp_norm[perm] for perm in perms])
    curr_mt_perm = np.asarray([curr_mt_cell_norm[perm] for perm in perms])
    curr_mvavg_rng_comp = np.apply_along_axis(
        MovingAverages.mvavg, 1, curr_comp_perm, WINDOW
    )
    curr_mvavg_rng_mt = np.apply_along_axis(
        MovingAverages.mvavg, 1, curr_mt_perm, WINDOW
    )
    curr_percvar_rng_comp = np.var(curr_mvavg_rng_comp, axis=1) / np.var(
        curr_comp_perm, axis=1
    )
    curr_percvar_rng_mt = np.var(curr_mvavg_rng_mt, axis=1) / np.var(
        curr_mt_perm, axis=1
    )
    return (
        curr_comp_norm,
        curr_comp_percvar,
        curr_comp_mvavg,
        curr_comp_perm,
        curr_mt_perm,
        curr_mvavg_rng_comp,
        curr_mvavg_rng_mt,
        curr_percvar_rng_comp,
        curr_percvar_rng_mt,
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

    def copy_mvavg_plots_protein(self):
        """Copy the plots generated for each gene to a more informative location"""
        folder = self.folder
        wp_ensg = self.wp_ensg
        wp_comp_ccd_difffromrng = self.wp_comp_ccd_difffromrng
        wp_isbimodal_fcpadj_pass = self.protein_bimodality.wp_isbimodal_fcpadj_pass
        wp_comp_ccd_clust1 = self.wp_comp_ccd_clust1
        wp_comp_ccd_clust2 = self.wp_comp_ccd_clust2
        wp_ccd_unibimodal = self.wp_ccd_unibimodal
        wp_comp_ccd_gauss = self.wp_comp_ccd_gauss

        fileprefixes = get_fileprefixes(wp_ensg)

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
        for ensg in fileprefixes[wp_comp_ccd_difffromrng]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdunifolder, ensg + "_mvavg.pdf"),
            )

        # CCD Unimodal and Bimodal
        for ensg in fileprefixes[
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
        for ensg in fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_clust1_mvavg.pdf"),
            )
        for ensg in fileprefixes[wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust2_mvavg.pdf"),
                os.path.join(ccdunibifolder, ensg + "_clust2_mvavg.pdf"),
            )

        # CCD Bimodal
        for ensg in fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust1]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust1_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1_mvavg.pdf"),
            )
        for ensg in fileprefixes[~wp_comp_ccd_difffromrng & wp_comp_ccd_clust2]:
            shutil.copy(
                os.path.join(folder, ensg + "_clust1&2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust1&2_mvavg.pdf"),
            )
            shutil.copy(
                os.path.join(folder, ensg + "_clust2_mvavg.pdf"),
                os.path.join(ccdpbifolder, ensg + "_clust2_mvavg.pdf"),
            )

        # Non-CCD
        for ensg in fileprefixes[
            ~wp_comp_ccd_difffromrng & ~wp_comp_ccd_clust1 & ~wp_comp_ccd_clust2
        ]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(nonccdfolder, ensg + "_mvavg.pdf"),
            )

        # Non-CCD Bimodal
        for ensg in fileprefixes[
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
        for ensg in fileprefixes[wp_comp_ccd_gauss & wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdgaussccdfolder, ensg + "_mvavg.pdf"),
            )
        # Gauss and Non-CCD
        for ensg in fileprefixes[wp_comp_ccd_gauss & ~wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(ccdgaussnonccdfolder, ensg + "_mvavg.pdf"),
            )
        # Non-Gauss and CCD
        for ensg in fileprefixes[~wp_comp_ccd_gauss & wp_ccd_unibimodal]:
            shutil.copy(
                os.path.join(folder, ensg + "_mvavg.pdf"),
                os.path.join(nongaussccdfolder, ensg + "_mvavg.pdf"),
            )

    def calculate_cell_cycle_dependence_protein(
        self,
        alpha_ccd=0.01,
        use_log_ccd=False,
        keep_outliers=False,
    ):
        """
        Use a moving average model of protein expression over the cell cycle to determine cell cycle dependence.
        Generates plots for each gene
        """
        self.do_remove_outliers = not keep_outliers

        # percent variance attributed to cell cycle (mean POI intensities)
        self.perc_var_cell = []
        self.perc_var_nuc = []
        self.perc_var_cyto = []
        self.perc_var_mt = []

        # randomized in pseudotime; percent variances
        self.perc_var_comp_rng = []
        self.perc_var_mt_rng = []

        # percent variances for bimodal
        self.perc_var_comp_clust1 = []
        self.perc_var_comp_clust2 = []
        self.mvavgs_comp_clust1 = []
        self.mvavgs_comp_clust2 = []
        self.perc_var_comp_clust1_rng = []
        self.perc_var_comp_clust2_rng = []
        self.mvavgs_x_clust1 = []
        self.mvavgs_x_clust2 = []

        # moving average y values & x value
        self.mvavgs_comp = []
        self.mvavgs_mt = []
        self.mvavgs_x = []

        # for plotting dataframe
        self.curr_pols = []
        self.curr_ab_norms = []
        self.mvperc_comps = []
        self.curr_freds = []
        self.curr_fgreens = []
        self.curr_mockbulk_phases = []

        self.curr_area_cell = []
        self.curr_area_nuc = []
        self.curr_well_plate_imgnb = []
        cell_counts = []

        self.folder = "figures/TemporalMovingAverages"
        folder_mt = "figures/TemporalMovingAveragesMicrotubules"
        folder_rng = "figures/TemporalMovingAverageRandomizationExamples"
        if not os.path.exists(self.folder):
            os.mkdir(self.folder)
        if not os.path.exists(folder_mt):
            os.mkdir(folder_mt)
        if not os.path.exists(folder_rng):
            os.mkdir(folder_rng)
        fileprefixes = get_fileprefixes(self.wp_ensg)

        for i, well in enumerate(self.u_well_plates):
            #    print(well)
            plt.close("all")
            if i % 100 == 0:
                print(f"well {i} of {len(self.u_well_plates)}")
            #    well = 'H05_55405991'#GMNN well, used for testing
            curr_well_inds = (
                self.fucci_polar_coords.pol_sort_well_plate == well
            )  # the reversal isn't really helpful here
            curr_pol = self.fucci_polar_coords.pol_sort_norm_rev[curr_well_inds]
            curr_ab_cell = (
                self.fucci_polar_coords.pol_sort_ab_cell[curr_well_inds]
                if not use_log_ccd
                else np.log10(self.fucci_polar_coords.pol_sort_ab_cell[curr_well_inds])
            )
            curr_ab_nuc = (
                self.fucci_polar_coords.pol_sort_ab_nuc[curr_well_inds]
                if not use_log_ccd
                else np.log10(self.fucci_polar_coords.pol_sort_ab_nuc[curr_well_inds])
            )
            curr_ab_cyto = (
                self.fucci_polar_coords.pol_sort_ab_cyto[curr_well_inds]
                if not use_log_ccd
                else np.log10(self.fucci_polar_coords.pol_sort_ab_cyto[curr_well_inds])
            )
            curr_mt_cell = (
                self.fucci_polar_coords.pol_sort_mt_cell[curr_well_inds]
                if not use_log_ccd
                else np.log10(self.fucci_polar_coords.pol_sort_mt_cell[curr_well_inds])
            )
            curr_fred = self.fucci_polar_coords.pol_sort_fred[curr_well_inds]
            curr_fgreen = self.fucci_polar_coords.pol_sort_fgreen[curr_well_inds]
            curr_mockbulk_phase = self.fucci_polar_coords.pol_sort_mockbulk_phases[
                curr_well_inds
            ]
            if self.do_remove_outliers:
                curr_comp = (
                    curr_ab_cell
                    if self.wp_iscell[i]
                    else curr_ab_nuc
                    if self.wp_isnuc[i]
                    else curr_ab_cyto
                )
                curr_pol = MovingAverages.remove_outliers(curr_comp, curr_pol)
                curr_ab_cell = MovingAverages.remove_outliers(curr_comp, curr_ab_cell)
                curr_ab_nuc = MovingAverages.remove_outliers(curr_comp, curr_ab_nuc)
                curr_ab_cyto = MovingAverages.remove_outliers(curr_comp, curr_ab_cyto)
                curr_mt_cell = MovingAverages.remove_outliers(curr_comp, curr_mt_cell)
                curr_fred = MovingAverages.remove_outliers(curr_comp, curr_fred)
                curr_fgreen = MovingAverages.remove_outliers(curr_comp, curr_fgreen)
                curr_mockbulk_phase = MovingAverages.remove_outliers(
                    curr_comp, curr_mockbulk_phase
                )

            # Normalize mean intensities, normalized for display
            curr_ab_cell_norm = curr_ab_cell / np.max(curr_ab_cell)
            curr_ab_nuc_norm = curr_ab_nuc / np.max(curr_ab_nuc)
            curr_ab_cyto_norm = curr_ab_cyto / np.max(curr_ab_cyto)
            curr_mt_cell_norm = curr_mt_cell / np.max(curr_mt_cell)

            # Original method from Devin's work
            perc_var_cell_val, mvavg_cell = MovingAverages.mvavg_perc_var(
                curr_ab_cell_norm, WINDOW
            )
            perc_var_nuc_val, mvavg_nuc = MovingAverages.mvavg_perc_var(
                curr_ab_nuc_norm, WINDOW
            )
            perc_var_cyto_val, mvavg_cyto = MovingAverages.mvavg_perc_var(
                curr_ab_cyto_norm, WINDOW
            )
            perc_var_mt_val, mvavg_mt = MovingAverages.mvavg_perc_var(
                curr_mt_cell_norm, WINDOW
            )
            mvavg_xvals = MovingAverages.mvavg(curr_pol, WINDOW)

            # Permutation analysis
            (
                curr_comp_norm,
                curr_comp_percvar,
                curr_comp_mvavg,
                curr_comp_perm,
                curr_mt_perm,
                curr_mvavg_rng_comp,
                curr_mvavg_rng_mt,
                curr_percvar_rng_comp,
                curr_percvar_rng_mt,
            ) = permutation_analysis_protein(
                i,
                curr_pol,
                curr_ab_cell_norm,
                curr_ab_nuc_norm,
                curr_ab_cyto_norm,
                curr_mt_cell_norm,
                perc_var_cell_val,
                perc_var_nuc_val,
                perc_var_cyto_val,
                self.wp_iscell,
                self.wp_isnuc,
                self.wp_iscyto,
                mvavg_cell,
                mvavg_nuc,
                mvavg_cyto,
            )

            self.perc_var_comp_rng.append(curr_percvar_rng_comp)
            self.perc_var_mt_rng.append(curr_percvar_rng_mt)

            # Assess CCD in high- and low-expressing cell populations of proteins determined to have bimodal population intensity
            clusters = self.protein_bimodality.assess_bimodal_ccd(
                self, i, curr_comp, fileprefixes, WINDOW
            )

            # Make example plots for the randomization trials for the NFAT5 example in manuscript
            if self.wp_ensg[i] == "ENSG00000102908" and self.do_plotting:
                median_rng_idx = np.argsort(curr_percvar_rng_comp)
                for iii, idx in enumerate(
                    [
                        0,
                        len(curr_percvar_rng_comp) // 4,
                        len(curr_percvar_rng_comp) // 2,
                        3 * len(curr_percvar_rng_comp) // 4,
                        len(curr_percvar_rng_comp) - 1,
                    ]
                ):
                    MovingAverages.temporal_mov_avg_randomization_example_protein(
                        curr_pol,
                        curr_comp_norm,
                        curr_comp_perm[median_rng_idx[idx]],
                        mvavg_xvals,
                        curr_comp_mvavg,
                        curr_mvavg_rng_comp[median_rng_idx[idx]],
                        folder_rng,
                        f"{fileprefixes[i]}_withrandomization_{iii}",
                    )

            # Test for equal variances of the moving averages and raw values
            self.perc_var_cell.append(perc_var_cell_val)
            self.perc_var_nuc.append(perc_var_nuc_val)
            self.perc_var_cyto.append(perc_var_cyto_val)
            self.perc_var_mt.append(perc_var_mt_val)

            curr_ab_norm = (
                curr_ab_cell_norm
                if self.wp_iscell[i]
                else curr_ab_nuc_norm
                if self.wp_isnuc[i]
                else curr_ab_cyto_norm
            )
            mvavg_comp = (
                mvavg_cell
                if self.wp_iscell[i]
                else mvavg_nuc
                if self.wp_isnuc[i]
                else mvavg_cyto
            )
            self.mvavgs_comp.append(mvavg_comp)
            self.mvavgs_x.append(mvavg_xvals)
            self.curr_pols.append(curr_pol)
            self.curr_ab_norms.append(curr_ab_norm)
            self.curr_freds.append(curr_fred)
            self.curr_fgreens.append(curr_fgreen)
            self.curr_mockbulk_phases.append(curr_mockbulk_phase)
            self.curr_area_cell.append(
                self.fucci_polar_coords.pol_sort_area_cell[curr_well_inds]
            )
            self.curr_area_cell.append(
                self.fucci_polar_coords.pol_sort_area_nuc[curr_well_inds]
            )
            self.curr_well_plate_imgnb.append(
                self.fucci_polar_coords.pol_sort_well_plate_imgnb[curr_well_inds]
            )
            cell_counts.append(len(curr_pol))

            # Plotting bimodals and microtubules over pseudotime
            if self.do_plotting:
                # Make the plots for each protein (takes 10 mins)
                windows = np.asarray(
                    [
                        np.arange(start, start + WINDOW)
                        for start in np.arange(len(curr_pol) - WINDOW + 1)
                    ]
                )
                mvperc_comp = MovingAverages.mvpercentiles(curr_ab_norm[windows])
                self.mvperc_comps.append(mvperc_comp)
                MovingAverages.temporal_mov_avg_protein(
                    curr_pol,
                    curr_ab_norm,
                    mvavg_xvals,
                    mvavg_comp,
                    mvperc_comp,
                    None,
                    self.folder,
                    fileprefixes[i],
                )

                # Make the plots for microtubules (takes 10 mins for all, so just do the arbitrary one in the Fig 1)
                if well == "C07_55405991":
                    mvperc_mt = MovingAverages.mvpercentiles(curr_mt_cell_norm[windows])
                    MovingAverages.temporal_mov_avg_protein(
                        curr_pol,
                        curr_mt_cell_norm,
                        mvavg_xvals,
                        mvavg_mt,
                        mvperc_mt,
                        None,
                        folder_mt,
                        f"{fileprefixes[i]}_mt",
                    )
                    if (
                        well == "C07_55405991"
                    ):  # keep here in case making all pseudotime plots
                        pd.DataFrame(
                            {
                                "ENSG": self.wp_ensg[i],
                                "Antibody": self.wp_ab[i],
                                "Compartment": "Cell",
                                "CCD": "No",
                                "cell_pseudotime": [
                                    ",".join([str(ppp) for ppp in pp])
                                    for pp in [curr_pol]
                                ],
                                "cell_intensity": [
                                    ",".join([str(yyy) for yyy in yy])
                                    for yy in [curr_mt_cell_norm]
                                ],
                                "mvavg_x": [
                                    ",".join([str(xxx) for xxx in xx])
                                    for xx in [mvavg_xvals]
                                ],
                                "mvavg_y": [
                                    ",".join([str(yyy) for yyy in yy])
                                    for yy in [mvavg_mt]
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
                                "phase": [",".join(pp) for pp in [curr_mockbulk_phase]],
                                "WellPlate": self.u_well_plates[i],
                            }
                        ).to_csv("output/mtplottingline.tsv", index=False, sep="\t")

                # Make the plots for each bimodal protein
                if clusters is not None:
                    MovingAverages.temporal_mov_avg_protein(
                        curr_pol,
                        curr_ab_norm,
                        mvavg_xvals,
                        mvavg_comp,
                        mvperc_comp,
                        clusters,
                        self.folder,
                        fileprefixes[i] + "_clust1&2",
                    )

        # percent variance attributed to cell cycle (mean POI intensities)
        self.perc_var_cell = np.array(self.perc_var_cell)
        self.perc_var_nuc = np.array(self.perc_var_nuc)
        self.perc_var_cyto = np.array(self.perc_var_cyto)
        self.perc_var_mt = np.array(self.perc_var_mt)

        self.perc_var_mt_rng = np.array(self.perc_var_mt_rng)
        self.perc_var_mt_rngperc_var_comp_rng = np.array(self.perc_var_comp_rng)

        # Let's check out which percent variances are greater than the permuted values
        self.perc_var_comp = utils.values_comp(
            self.perc_var_cell,
            self.perc_var_nuc,
            self.perc_var_cyto,
            self.wp_iscell,
            self.wp_isnuc,
            self.wp_iscyto,
        )
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
        ccd_var_comp_rng_wilcoxp_withbimodal = (
            np.apply_along_axis(
                scipy.stats.wilcoxon,
                1,
                (perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T,
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
                (self.perc_var_mt - self.perc_var_mt_rng.T).T,
                None,
                "wilcox",
                False,
                "greater",
            )
            .T[1]
            .T
        )
        mean_diff_from_rng_mt = np.mean(
            (self.perc_var_mt - self.perc_var_mt_rng.T).T, 1
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
            (perc_var_comp_withbimodal - perc_var_comp_rng_withbimodal.T).T, 1
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
            self.protein_bimodality.wp_comp_kruskal_gaussccd_adj <= alpha_ccd
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

        wp_ccd_unibimodal = (
            wp_comp_ccd_difffromrng | wp_comp_ccd_clust1 | wp_comp_ccd_clust2
        )

        self.wp_comp_ccd_difffromrng = wp_comp_ccd_difffromrng
        self.mean_diff_from_rng_mt = mean_diff_from_rng_mt
        self.wp_comp_ccd_clust1 = wp_comp_ccd_clust1
        self.wp_comp_ccd_clust2 = wp_comp_ccd_clust2
        self.wp_ccd_unibimodal = wp_ccd_unibimodal
        self.wp_comp_ccd_gauss = wp_comp_ccd_gauss
        self.mean_diff_from_rng = mean_diff_from_rng
        self.wp_comp_eq_percvar_adj = wp_comp_eq_percvar_adj
        self.mean_diff_from_rng_clust1 = mean_diff_from_rng_clust1
        self.wp_comp_eq_percvar_adj_clust1 = wp_comp_eq_percvar_adj_clust1
        self.mean_diff_from_rng_clust2 = mean_diff_from_rng_clust2
        self.wp_comp_eq_percvar_adj_clust2 = wp_comp_eq_percvar_adj_clust2

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

    def analyze_ccd_variation_protein(self):
        """Analyze the cell cycle dependence for all the proteins"""
        ccd_comp = (
            self.wp_comp_ccd_difffromrng
            | self.wp_comp_ccd_clust1
            | self.wp_comp_ccd_clust2
        )
        nonccd_comp = ~ccd_comp
        bioccd = np.genfromtxt(
            "input/ProteinData/BiologicallyDefinedCCD.txt", dtype="str"
        )  # from mitotic structures

        # pickle the results (removed ones passing in only one replicate)
        utils.np_save_overwriting("output/pickles/ccd_comp.npy", ccd_comp)
        utils.np_save_overwriting("output/pickles/nonccd_comp.npy", nonccd_comp)

        self.ccd_comp = ccd_comp
        self.nonccd_comp = nonccd_comp
        self.bioccd = bioccd
