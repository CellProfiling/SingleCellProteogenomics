# -*- coding: utf-8 -*-
"""
Loading the raw immunofluorescence data for single-cell protein expression.
Filtering the data for:
    - Sample with low cell count
    - Out of focus images removed
    - Antibodies that failed validation are removed
    - Small nuclei ar removed (mitotic cells) prior to this analylsis
    - Large nuclei are removed (segmentation artificats)
    - Only proteins displaying variability are kept for single-cell variability analysis (using manual annotations)

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SingleCellProteogenomics import utils

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

# EMPTYWELLS: These wells on the last plate didn't have cells; the segmentation algorithm still annotated some, so remove them
EMPTYWELLS = set(
    [
        "B11_6745",
        "C11_6745",
        "D11_6745",
        "E11_6745",
        "F11_6745",
        "G11_6745",
        "H11_6745",
        "A12_6745",
        "B12_6745",
        "C12_6745",
        "D12_6745",
        "E12_6745",
        "F12_6745",
        "G12_6745",
    ]
)

# Minimum number of cells per sample required for cell cycle analysis with pseudotime
MIN_CELL_COUNT = 60

# 0: use mean, (testing intensity that's already normalized for cell size)
# 1: use integrated, (testing integrated because it reflects that small cells are often brighter because they have rounded up and are thicker)
# 2: use integrated / nucleus area (testing a normalization by nucleus size, which may be more consistent in segmentation)
INTENSITY_SWITCH = 0  # cell size does increase into G2 and that has a substantial effect, so use the mean intensity, which is well behaved


def plot_areas(areas, title):
    """Histogram for areas of cell/nuc/cytoplasm"""
    bins = plt.hist(areas, bins=100, alpha=0.5)
    plt.vlines(np.mean(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) - 2 * np.std(areas), 0, np.max(bins[0]))
    plt.vlines(np.mean(areas) + 2 * np.std(areas), 0, np.max(bins[0]))
    plt.title(title)
    plt.ylabel("Count")
    plt.xlabel("Area")
    plt.savefig(f"figures/areas{title}.png")
    plt.close()


class ProteinData:
    """Loading in protein data"""

    def __init__(self, do_plotting):
        self.do_plotting = do_plotting

    def read_raw_data(self):
        """Read in the raw protein IF data"""
        print("reading raw protein IF data")
        my_df1 = pd.read_csv("input/ProteinData/FucciDataFirstPlates.csv.gz")
        my_df2 = pd.read_csv("input/ProteinData/FucciDataSecondPlates.csv.gz")
        self.raw_protein_data = pd.concat((my_df1, my_df2), sort=True)
        print("loaded raw data")

    def read_sample_info(self, df):
        """Get the metadata for all the samples"""
        self.plate = np.asarray(df.plate)
        self.u_plate = np.unique(self.plate)
        self.well_plate = np.asarray(df.well_plate)
        self.imgnb = np.asarray(df.ImageNumber)
        self.well_plate_imgnb = np.asarray(
            [f"{wp}_{self.imgnb[i]}" for i, wp in enumerate(self.well_plate)]
        )
        self.u_well_plates = np.unique(self.well_plate)
        self.ab_objnum = np.asarray(df.ObjectNumber)
        self.well_plate_imgnb_objnb = np.asarray(
            [
                f"{wp}_{self.imgnb[i]}_{self.ab_objnum[i]}"
                for i, wp in enumerate(self.well_plate)
            ]
        )
        self.area_cell = np.asarray(df.Area_cell)
        self.area_nuc = np.asarray(df.AreaShape_Area)
        self.area_cyto = np.asarray(df.Area_cyto)

        name_df = pd.read_csv("input/ProteinData/FucciStainingSummaryFirstPlates.csv")
        wppp1 = list(name_df["well_plate"])
        ensggg1 = list(name_df["ENSG"])
        abbb1 = list(name_df["Antibody"])
        rrrr = list(name_df["Results_final_update"])
        cccc1 = list(name_df["Compartment"])

        name_df2 = pd.read_csv("input/ProteinData/FucciStainingSummarySecondPlates.csv")
        wppp2 = list(name_df2["well_plate"])
        ensggg2 = list(name_df2["ENSG"])
        abbb2 = list(name_df2["Antibody"])
        cccc2 = list(name_df2["Compartment"])

        wppp, ensggg, abbb, cccc = (
            wppp1 + wppp2,
            ensggg1 + ensggg2,
            abbb1 + abbb2,
            cccc1 + cccc2,
        )

        self.ensg_dict = dict([(wppp[i], ensggg[i]) for i in range(len(wppp))])
        self.ab_dict = dict([(wppp[i], abbb[i]) for i in range(len(wppp))])
        self.result_dict = dict(
            [(wppp1[i], rrrr[i]) for i in range(len(wppp1))]
        )  # Only available for first plates
        self.compartment_dict = dict([(wppp[i], cccc[i]) for i in range(len(wppp))])
        self.ENSG = np.asarray(
            [
                self.ensg_dict[wp] if wp in self.ensg_dict else ""
                for wp in self.well_plate
            ]
        )
        self.antibody = np.asarray(
            [self.ab_dict[wp] if wp in self.ab_dict else "" for wp in self.well_plate]
        )
        self.result = np.asarray(
            [
                self.result_dict[wp] if wp in self.result_dict else ""
                for wp in self.well_plate
            ]
        )
        self.compartment = np.asarray(
            [
                self.compartment_dict[wp] if wp in self.compartment_dict else ""
                for wp in self.well_plate
            ]
        )

        # Pickle the results
        if not os.path.exists("output/"):
            os.mkdir("output/")
        if not os.path.exists("output/pickles/"):
            os.mkdir("output/pickles/")
        if not os.path.exists("figures/"):
            os.mkdir("figures/")
        utils.np_save_overwriting("output/pickles/plate.npy", self.plate)
        utils.np_save_overwriting("output/pickles/u_plate.npy", self.u_plate)
        utils.np_save_overwriting(
            "output/pickles/u_well_plates.npy", self.u_well_plates
        )
        utils.np_save_overwriting("output/pickles/area_cell.npy", self.area_cell)
        utils.np_save_overwriting("output/pickles/area_nuc.npy", self.area_nuc)
        utils.np_save_overwriting("output/pickles/area_cyto.npy", self.area_cyto)
        utils.np_save_overwriting("output/pickles/well_plate.npy", self.well_plate)
        utils.np_save_overwriting(
            "output/pickles/well_plate_imgnb.npy", self.well_plate_imgnb
        )
        utils.np_save_overwriting(
            "output/pickles/well_plate_imgnb_objnb.npy", self.well_plate_imgnb_objnb
        )

    def previous_results(self):
        """Process the results metadata into lists of previously annotated CCD proteins"""
        wp_ensg = np.asarray(
            [
                self.ensg_dict[wp] if wp in self.ensg_dict else ""
                for wp in self.u_well_plates
            ]
        )
        wp_ab = np.asarray(
            [
                self.ab_dict[wp] if wp in self.ab_dict else ""
                for wp in self.u_well_plates
            ]
        )
        wp_prev_ccd = np.asarray(
            [
                wp in self.result_dict and self.result_dict[wp].startswith("ccd")
                for wp in self.u_well_plates
            ]
        )
        wp_prev_notccd = np.asarray(
            [
                wp in self.result_dict and self.result_dict[wp].startswith("notccd")
                for wp in self.u_well_plates
            ]
        )
        wp_prev_negative = np.asarray(
            [
                wp in self.result_dict and self.result_dict[wp].startswith("negative")
                for wp in self.u_well_plates
            ]
        )
        prev_ccd_ensg = wp_ensg[wp_prev_ccd]
        prev_notccd_ensg = wp_ensg[wp_prev_notccd]
        prev_negative_ensg = wp_ensg[wp_prev_negative]

        # Pickle the results
        utils.np_save_overwriting("output/pickles/wp_ensg.npy", wp_ensg)
        utils.np_save_overwriting("output/pickles/wp_ab.npy", wp_ab)
        utils.np_save_overwriting("output/pickles/wp_prev_ccd.npy", wp_prev_ccd)
        utils.np_save_overwriting("output/pickles/wp_prev_notccd.npy", wp_prev_notccd)
        utils.np_save_overwriting(
            "output/pickles/wp_prev_negative.npy", wp_prev_negative
        )
        utils.np_save_overwriting("output/pickles/prev_ccd_ensg.npy", prev_ccd_ensg)
        utils.np_save_overwriting(
            "output/pickles/prev_notccd_ensg.npy", prev_notccd_ensg
        )
        utils.np_save_overwriting(
            "output/pickles/prev_negative_ensg.npy", prev_negative_ensg
        )

        self.wp_ensg = wp_ensg
        self.wp_ab = wp_ab
        self.wp_prev_ccd = wp_prev_ccd
        self.wp_prev_notccd = wp_prev_notccd
        self.wp_prev_negative = wp_prev_negative
        self.prev_ccd_ensg = prev_ccd_ensg
        self.prev_notccd_ensg = prev_notccd_ensg
        self.prev_negative_ensg = prev_negative_ensg

    def apply_manual_filtering(self):
        """Filter raw data based on manual annotations"""
        my_df = self.raw_protein_data
        result_dict = self.result_dict
        ab_dict = self.ab_dict

        # filter some wells in the last plate didn't have anything.
        before = f"{len(my_df)}: number of cells before filtering empty wells"
        print(before)
        my_df = my_df[~my_df.well_plate.isin(EMPTYWELLS)]
        after = f"{len(my_df)}: number of cells after filtering empty wells"
        print(after)

        print()
        print("filtering out of focus")
        my_df_filtered = my_df
        oof = pd.read_csv("input/ProteinData/OutOfFocusImages.txt", header=None)[0]
        well_plate = np.asarray(my_df_filtered.well_plate)
        imgnb = np.asarray(my_df_filtered.ImageNumber)
        well_plate_imgnb = np.asarray(
            [f"{wp}_{imgnb[i]}" for i, wp in enumerate(well_plate)]
        )
        before = f"{len(my_df_filtered)}: number of cells before filtering out of focus images"
        print(before)
        my_df_filtered = my_df_filtered[~np.isin(well_plate_imgnb, oof)]
        after = f"{len(my_df_filtered)}: number of cells after filtering out of focus images"
        print(f"{after}\nfinished filtering")

        print()
        print("filtering negative staining")
        new_data_or_nonnegative_stain = [
            wp not in result_dict
            or (
                not result_dict[wp].lower().startswith("negative")
                and not wp.startswith("H12")
            )
            for wp in my_df_filtered.well_plate
        ]
        before = f"{len(my_df_filtered)}: number of cells before filtering negative staining from first batch"
        print(before)
        my_df_filtered = my_df_filtered[new_data_or_nonnegative_stain]
        after = f"{len(my_df_filtered)}: number of cells after filtering negative staining from first batch"
        print(f"{after}\nfinished filtering")

        print()
        print("filtering bad fields of view (negative staining, unspecific, etc)")
        filterthese = pd.read_csv("input/ProteinData/FOV_ImgNum_Lookup.csv")
        badfov = filterthese["well_plate_imgnb"][(filterthese["UseImage"] == 0)]
        well_plate = np.asarray(my_df_filtered.well_plate)
        imgnb = np.asarray(my_df_filtered.ImageNumber)
        well_plate_imgnb = np.asarray(
            [f"{wp}_{imgnb[i]}" for i, wp in enumerate(well_plate)]
        )
        negative_controls = np.asarray([wp.startswith("H12") for wp in well_plate])
        before = f"{len(my_df_filtered)}: number of cells before filtering out of focus images"
        print(before)
        my_df_filtered = my_df_filtered[
            ~np.isin(well_plate_imgnb, badfov) & ~negative_controls
        ]
        after = f"{len(my_df_filtered)}: number of cells after filtering out of focus images"
        print(f"{after}\nfinished filtering")

        print()
        print("filtering failed antibodies")
        failedab = np.genfromtxt(
            "input/ProteinData/RecentlyFailedAntibodies.txt", dtype="str"
        )
        before = f"{len(my_df_filtered)}: number of cells before filtering antibodies failed in HPAv19"
        print(before)
        my_df_filtered = my_df_filtered[
            ~np.isin([ab_dict[wp] for wp in my_df_filtered.well_plate], failedab)
        ]
        after = f"{len(my_df_filtered)}: number of cells after filtering antibodies failed in HPAv19"
        print(f"{after}\nfinished filtering")

        print()
        print("filtering mitotic proteins")
        mitoticab = np.genfromtxt(
            "input/ProteinData/RemoveMitoticAndMicrotubules.txt", dtype="str"
        )
        before = f"{len(my_df_filtered)}: number of cells before filtering mitotic/microtubule proteins"
        print(before)
        my_df_filtered = my_df_filtered[
            ~np.isin([ab_dict[wp] for wp in my_df_filtered.well_plate], mitoticab)
        ]
        after = f"{len(my_df_filtered)}: number of cells after filtering mitotic/microtubule proteins"
        print(f"{after}\nfinished filtering")

        self.filtered_protein_data = my_df_filtered

    def apply_big_nucleus_filter(self):
        """Filter the super big nuclei"""
        my_df = self.filtered_protein_data

        area_cell = my_df.Area_cell
        area_nuc = my_df.AreaShape_Area
        area_cyto = my_df.Area_cyto

        if self.do_plotting:
            plot_areas(area_cell, "area_cell")
            plot_areas(area_nuc, "area_nuc")
            plot_areas(area_cyto, "area_cyto")

        upper_nucleus_cutoff = np.mean(area_nuc) + 2 * np.std(area_nuc)

        print("filtering super big nuclei")
        my_df_filtered = my_df
        cell_passes_nucleus_filter = (
            my_df_filtered.AreaShape_Area < upper_nucleus_cutoff
        )
        before = f"{len(my_df_filtered)}: number of cells before filtering out super big nuclei"
        print(before)
        my_df_filtered = my_df_filtered[cell_passes_nucleus_filter]
        after = f"{len(my_df_filtered)}: number of cells after filtering out super big nuclei"
        print(f"{after}\nfinished filtering on nuclei")

        area_cell_filtered = my_df_filtered.Area_cell
        area_nuc_filtered = my_df_filtered.AreaShape_Area
        area_cyto_filtered = my_df_filtered.Area_cyto

        if self.do_plotting:
            plot_areas(area_cell_filtered, "area_cell_filtered")
            plot_areas(area_nuc_filtered, "area_nuc_filtered")
            plot_areas(area_cyto_filtered, "area_cyto_filtered")

        self.filtered_protein_data = my_df_filtered

    def apply_cell_count_filter(self):
        """Filter low cell counts per sample"""
        my_df_filtered = self.filtered_protein_data
        well_plate = np.asarray(my_df_filtered.well_plate)
        u_well_plates = np.unique(my_df_filtered.well_plate)
        cell_count_dict = {}
        for wp in well_plate:
            if wp in cell_count_dict:
                cell_count_dict[wp] += 1
            else:
                cell_count_dict[wp] = 1
        cell_counts = np.array([cell_count_dict[wp] for wp in well_plate])
        print("filtering low cell counts")
        my_df_filtered = my_df_filtered[cell_counts >= MIN_CELL_COUNT]
        before = f"{len(self.raw_protein_data)}: number of cells before filtering out samples with < {MIN_CELL_COUNT} cells"
        after = f"{len(my_df_filtered)}: number of cells after filtering out samples with < {MIN_CELL_COUNT} cells"
        print(f"{before}\n{after}\nfinished filtering on cell count")
        self.filtered_protein_data = my_df_filtered

    def apply_variation_filter(self, output_dataframes=False):
        """Separate the varying and nonvarying samples"""
        my_df_filtered = self.filtered_protein_data
        result_dict, unfiltered_df = self.result_dict, self.raw_protein_data
        variable_firstbatch = np.asarray(
            [
                wp in result_dict
                and not result_dict[wp].replace(" ", "").startswith("novariation")
                for wp in my_df_filtered.well_plate
            ]
        )

        varann_secondbatch = pd.read_csv(
            "input/ProteinData/SecondBatchVariableLookup.csv"
        )
        variable_ann_secondbatch = np.asarray(
            [
                str(vv).lower().startswith("yes")
                for vv in varann_secondbatch["IsVariable"]
            ]
        )
        variable_wp_secondbatch = np.asarray(
            varann_secondbatch["well_plate"][variable_ann_secondbatch]
        )
        variable_secondbatch = np.isin(
            my_df_filtered.well_plate, variable_wp_secondbatch
        )

        self.my_df_filtered_variation = my_df_filtered[
            variable_firstbatch | variable_secondbatch
        ]
        self.my_df_filtered_novariation = my_df_filtered[
            ~(variable_firstbatch | variable_secondbatch)
        ]
        before = f"{len(unfiltered_df)}: number of cells before filtering for variation"
        withvar = f"{len(self.my_df_filtered_variation)}: number of cells in samples with variation"
        withoutvar = f"{len(self.my_df_filtered_novariation)}: number of cells in samples without variation"
        print(f"{before}\n{withvar}\n{withoutvar}")

        if output_dataframes:
            """output these dataframes (used for skewness / kurtosis analysis)"""
            self.my_df_filtered_variation.to_csv(
                "output/nuc_predicted_prob_phases_filtered_variation.csv"
            )
            self.my_df_filtered_novariation.to_csv(
                "output/nuc_predicted_prob_phases_filtered_novariation.csv"
            )

    def metacompartments(self, my_df_filtered_variation):
        """Get the compartments for the unique wellplates"""
        u_well_plates, compartment_dict = self.u_well_plates, self.compartment_dict
        wp_iscell = np.asarray(
            [
                compartment_dict[wp].lower().startswith("cell")
                if wp in compartment_dict
                else False
                for wp in u_well_plates
            ]
        )
        wp_isnuc = np.asarray(
            [
                compartment_dict[wp].lower().startswith("nuc")
                if wp in compartment_dict
                else False
                for wp in u_well_plates
            ]
        )
        wp_iscyto = np.asarray(
            [
                compartment_dict[wp].lower().startswith("cyto")
                if wp in compartment_dict
                else False
                for wp in u_well_plates
            ]
        )

        # Pickle the results
        utils.np_save_overwriting("output/pickles/wp_iscell.npy", wp_iscell)
        utils.np_save_overwriting("output/pickles/wp_isnuc.npy", wp_isnuc)
        utils.np_save_overwriting("output/pickles/wp_iscyto.npy", wp_iscyto)

        wp_nocompartmentinfo = ~wp_iscell & ~wp_isnuc & ~wp_iscyto
        withoutcomp = f"{sum(wp_nocompartmentinfo)}: samples without compartment information; to be filtered since they're biologically defined as CCD and not included in the analysis"
        before = f"{len(my_df_filtered_variation)}: number of cells before filtering for compartment information"
        print(f"{withoutcomp}\n{before}")
        my_df_filtered_compartmentvariation = my_df_filtered_variation[
            ~np.isin(
                my_df_filtered_variation.well_plate, u_well_plates[wp_nocompartmentinfo]
            )
        ]
        after = f"{len(my_df_filtered_compartmentvariation)}: number of cells after filtering for compartment information"
        print(after)

        self.wp_iscell, self.wp_isnuc, self.wp_iscyto = wp_iscell, wp_isnuc, wp_iscyto
        self.my_df_filtered_compartmentvariation = my_df_filtered_compartmentvariation

    def read_sample_data(self, df):
        """Read antibody intensity data for each sample and save it to a file for later use."""
        # Antibody data (mean intensity)
        ab_nuc = np.asarray(
            [
                df.Intensity_MeanIntensity_ResizedAb,
                df.Intensity_IntegratedIntensity_ResizedAb,
                df.Intensity_IntegratedIntensity_ResizedAb / df.AreaShape_Area,
            ][INTENSITY_SWITCH]
        )
        ab_cyto = np.asarray(
            [
                df.Mean_ab_Cyto,
                df.Integrated_ab_cyto,
                df.Integrated_ab_cyto / df.AreaShape_Area,
            ][INTENSITY_SWITCH]
        )
        ab_cell = np.asarray(
            [
                df.Mean_ab_cell,
                df.Integrated_ab_cell,
                df.Integrated_ab_cell / df.AreaShape_Area,
            ][INTENSITY_SWITCH]
        )
        mt_cell = np.asarray(
            [
                df.Mean_mt_cell,
                df.Integrated_mt_cell,
                df.Integrated_mt_cell / df.AreaShape_Area,
            ][INTENSITY_SWITCH]
        )

        # Fucci data (mean intensity)
        green_fucci = np.asarray(df.Intensity_MeanIntensity_CorrResizedGreenFUCCI)
        red_fucci = np.asarray(df.Intensity_MeanIntensity_CorrResizedRedFUCCI)

        # Pickle the results
        utils.np_save_overwriting("output/pickles/ab_nuc.npy", ab_nuc)
        utils.np_save_overwriting("output/pickles/ab_cyto.npy", ab_cyto)
        utils.np_save_overwriting("output/pickles/ab_cell.npy", ab_cell)
        utils.np_save_overwriting("output/pickles/mt_cell.npy", mt_cell)
        utils.np_save_overwriting("output/pickles/green_fucci.npy", green_fucci)
        utils.np_save_overwriting("output/pickles/red_fucci.npy", red_fucci)

        (
            self.ab_nuc,
            self.ab_cyto,
            self.ab_cell,
            self.mt_cell,
            self.green_fucci,
            self.red_fucci,
        ) = (ab_nuc, ab_cyto, ab_cell, mt_cell, green_fucci, red_fucci)
