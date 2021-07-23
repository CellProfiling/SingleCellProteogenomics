# -*- coding: utf-8 -*-
"""
Includes modeling of the cell cycle using FUCCI cell cycle markers.
The FUCCI markers are modeled using a polar coordinate model of the CDT1 and GMNN log-intensities.
The markers were measured for each cell in protein measurements using confocal microscopy.
The markers were measured for each cell in RNA analysis using FACS intentensities.

@author: devinsullivan
@author: Anthony J. Cesnik, cesnik [at] kth.se
"""

from SingleCellProteogenomics import utils, stretch_time
from SingleCellProteogenomics.MovingAverages import mvpercentiles, mvavg
from SingleCellProteogenomics.FucciCellCycle import FucciCellCycle
import scipy.optimize
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import copy, shutil, decimal

# Make PDF text readable
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.dpi"] = 300

# Number of bins for polar coord calculation, arbitrary choice for now
NBINS_POLAR_COORD = 150

# Window used for visualizing FUCCI intensities over pseudotime
WINDOW_FUCCI_PSEUDOTIME = 100

# Object representing FUCCI cell cycle phase durations
fucci = FucciCellCycle()  

## POLAR COORDINATE FUNCTIONS
def calc_R(xc, yc, x, y):
    """Calculate the distance of each 2D points from the center (xc, yc)"""
    return np.sqrt((x - xc) ** 2 + (y - yc) ** 2)


def f_2(c, x, y):
    """Calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc)"""
    print(c)
    Ri = calc_R(c[0], c[1], x, y)
    return Ri - Ri.mean()


def cart2pol(x, y):
    """Convert cartesian coordinates to polar coordinates"""
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol_sort(inds, more_than_start, less_than_start, *args):
    """Sort data by polar coordinates and reorder based on the start position of the polar coordinate model"""
    return [
        np.concatenate((arr[inds][more_than_start], arr[inds][less_than_start]))
        for arr in args
    ]


def pol2cart(rho, phi):
    """Apply uniform radius (rho) and convert back"""
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


## PLOTTING HELPERS
def plot_annotate_time(R_2, start_phi, fraction):
    """Pseudotime annotation helper for point on plot"""
    pt = pol2cart(R_2, start_phi + (1 - fraction) * 2 * np.pi)
    plt.scatter(pt[0], pt[1], c="c", linewidths=4)
    plt.annotate(f"  {round(fraction * fucci.TOT_LEN, 2)} hrs", (pt[0], pt[1]))


def drange(x, y, jump):
    """Increment `x` by a decimal `jump` until greater than `y`"""
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)


def get_phase_colormap():
    """Get colormap used for cell cycle phase illustrations"""
    return {"G1": "blue", "G2M": "orange", "S-ph": "green"}


class FucciPolarCoords:
    def __init__(self):
        pass

    def fucci_hist2d(
        self,
        analysis_title,
        nbins=200,
        show_gmnn=False,
        pol_sort_well_plate=[],
        pol_sort_ab_nuc=[],
        pol_sort_centered_data0=[],
        pol_sort_centered_data1=[],
    ):
        """
        Visualize the log-FUCCI intensities and phase transitions.
        If `show_gmnn` is true, generate an overlay of the GMNN antibody staining intensities.
        """
        centered_data = self.centered_data
        cart_data_ur = self.cart_data_ur
        start_pt = self.start_pt
        g1_end_pt = self.g1_end_pt
        g1s_end_pt = self.g1s_end_pt
        R_2 = self.R_2
        start_phi = self.start_phi

        fig, ax1 = plt.subplots(figsize=(10, 10))
        mycmap = copy.copy(plt.cm.get_cmap("gray_r"))
        mycmap.set_under(color="w", alpha=None)
        ax1.hist2d(
            centered_data[:, 0], centered_data[:, 1], bins=nbins, alpha=1, cmap=mycmap
        )
        hist, xbins, ybins = np.histogram2d(
            cart_data_ur[0], cart_data_ur[1], bins=nbins, density=True
        )
        extent = [xbins.min(), xbins.max(), ybins.min(), ybins.max()]
        im = ax1.imshow(
            np.ma.masked_where(hist == 0, hist).T,
            interpolation="nearest",
            origin="lower",
            extent=extent,
            cmap="plasma",
        )
        if (
            show_gmnn
        ):  # GMNN was tagged in the FUCCI cells and stained with antibodies; visualize the agreement
            gmnn = "H05_55405991"
            gmnn_well_inds = pol_sort_well_plate == gmnn
            gmnn_ab_nuc = pol_sort_ab_nuc[gmnn_well_inds]
            im = ax1.scatter(
                pol_sort_centered_data0[gmnn_well_inds],
                pol_sort_centered_data1[gmnn_well_inds],
                c=gmnn_ab_nuc,
            )
            fig.colorbar(im, ax=ax1)
        else:
            plt.scatter(start_pt[0], start_pt[1], c="c", linewidths=4)
            plt.scatter(g1_end_pt[0], g1_end_pt[1], c="c", linewidths=4)
            plt.scatter(g1s_end_pt[0], g1s_end_pt[1], c="c", linewidths=4)
            plt.scatter(0, 0, c="m", linewidths=4)
            plt.annotate("  0 hrs (start)", (start_pt[0], start_pt[1]))
            plt.annotate(
                f"  {fucci.G1_LEN} hrs (end of G1)", (g1_end_pt[0], g1_end_pt[1])
            )
            plt.annotate(
                f"  {fucci.G1_LEN + fucci.G1_S_TRANS} hrs (end of S)",
                (g1s_end_pt[0], g1s_end_pt[1]),
            )
            for yeah in list(drange(decimal.Decimal(0.1), 0.9, "0.1")):
                plot_annotate_time(R_2, start_phi, yeah)
        plt.xlabel(r"$\propto log_{10}(GMNN_{fucci})$", size=20)
        plt.ylabel(r"$\propto log_{10}(CDT1_{fucci})$", size=20)
        plt.tight_layout()
        if show_gmnn:
            plt.savefig("figures/GMNN_FUCCI_plot.pdf", transparent=True)
        else:
            plt.savefig(
                f"figures/masked_polar_hist_{analysis_title}.pdf", transparent=True
            )
        plt.close()

    def fucci_polar_coords(self, x, y, analysis_title):
        """
        Calculate the polar coordinate position of each cell based on the FUCCI intensities (x, y).
        """
        fucci_data = np.column_stack([x, y])
        center_est_xy = np.mean(x), np.mean(y)
        center_est2_xy = scipy.optimize.least_squares(f_2, center_est_xy, args=(x, y))
        xc_2, yc_2 = center_est2_xy.x
        Ri_2 = calc_R(*center_est2_xy.x, x, y)
        R_2 = Ri_2.mean()
        residu_2 = sum((Ri_2 - R_2) ** 2)

        # Center data
        centered_data = fucci_data - center_est2_xy.x

        pol_data = cart2pol(centered_data[:, 0], centered_data[:, 1])
        pol_sort_inds = np.argsort(pol_data[1])
        pol_sort_rho = pol_data[0][pol_sort_inds]
        pol_sort_phi = pol_data[1][pol_sort_inds]
        centered_data_sort0 = centered_data[pol_sort_inds, 0]
        centered_data_sort1 = centered_data[pol_sort_inds, 1]

        # Rezero to minimum --resoning, cells disappear during mitosis, so we should have the fewest detected cells there
        bins = plt.hist(pol_sort_phi, NBINS_POLAR_COORD)
        start_phi = bins[1][np.argmin(bins[0])]

        # Move those points to the other side
        more_than_start = np.greater(pol_sort_phi, start_phi)
        less_than_start = np.less_equal(pol_sort_phi, start_phi)
        pol_sort_rho_reorder = np.concatenate(
            (pol_sort_rho[more_than_start], pol_sort_rho[less_than_start])
        )
        pol_sort_inds_reorder = np.concatenate(
            (pol_sort_inds[more_than_start], pol_sort_inds[less_than_start])
        )
        pol_sort_phi_reorder = np.concatenate(
            (pol_sort_phi[more_than_start], pol_sort_phi[less_than_start] + np.pi * 2)
        )
        pol_sort_centered_data0 = np.concatenate(
            (centered_data_sort0[more_than_start], centered_data_sort0[less_than_start])
        )
        pol_sort_centered_data1 = np.concatenate(
            (centered_data_sort1[more_than_start], centered_data_sort1[less_than_start])
        )
        pol_sort_shift = pol_sort_phi_reorder + np.abs(np.min(pol_sort_phi_reorder))

        # Shift and re-scale "time"
        # reverse "time" since the cycle goes counter-clockwise wrt the fucci plot
        pol_sort_norm = pol_sort_shift / np.max(pol_sort_shift)
        pol_sort_norm_rev = 1 - pol_sort_norm
        pol_sort_norm_rev = stretch_time.stretch_time(pol_sort_norm_rev)
        plt.tight_layout()
        plt.savefig(f"figures/FucciAllPseudotimeHist_{analysis_title}.png")
        plt.close()

        # visualize that result
        start_pt = pol2cart(R_2, start_phi)
        g1_end_pt = pol2cart(R_2, start_phi + (1 - fucci.G1_PROP) * 2 * np.pi)
        g1s_end_pt = pol2cart(R_2, start_phi + (1 - fucci.G1_S_PROP) * 2 * np.pi)
        cart_data_ur = pol2cart(np.repeat(R_2, len(centered_data)), pol_data[1])

        self.pol_sort_norm_rev = pol_sort_norm_rev
        self.centered_data = centered_data
        self.pol_sort_centered_data0 = pol_sort_centered_data0
        self.pol_sort_centered_data1 = pol_sort_centered_data1
        self.pol_sort_phi = pol_sort_phi
        self.pol_sort_inds = pol_sort_inds
        self.pol_sort_inds_reorder = pol_sort_inds_reorder
        self.pol_sort_phi_reorder = pol_sort_phi_reorder
        self.more_than_start = more_than_start
        self.less_than_start = less_than_start
        self.start_pt = start_pt
        self.g1_end_pt = g1_end_pt
        self.g1s_end_pt = g1s_end_pt
        self.cart_data_ur = cart_data_ur
        self.R_2 = R_2
        self.start_phi = start_phi

        # plot it
        if self.do_plotting:
            self.fucci_hist2d(analysis_title)

    def plot_fucci_intensities_on_pseudotime(self):
        """Visualize FUCCI intensities over pseudotime"""
        pol_sort_norm_rev = self.pol_sort_norm_rev
        pol_sort_centered_data1 = self.pol_sort_centered_data1
        pol_sort_centered_data0 = self.pol_sort_centered_data0

        plt.figure(figsize=(5, 5))
        WINDOW_FUCCI_PSEUDOTIMEs = np.asarray(
            [
                np.arange(start, start + WINDOW_FUCCI_PSEUDOTIME)
                for start in np.arange(
                    len(pol_sort_norm_rev) - WINDOW_FUCCI_PSEUDOTIME + 1
                )
            ]
        )
        mvperc_red = mvpercentiles(pol_sort_centered_data1[WINDOW_FUCCI_PSEUDOTIMEs])
        mvperc_green = mvpercentiles(pol_sort_centered_data0[WINDOW_FUCCI_PSEUDOTIMEs])
        mvavg_xvals = mvavg(pol_sort_norm_rev, WINDOW_FUCCI_PSEUDOTIME)
        plt.fill_between(
            mvavg_xvals * fucci.TOT_LEN,
            mvperc_green[1],
            mvperc_green[-2],
            color="lightgreen",
            label="25th & 75th Percentiles",
        )
        plt.fill_between(
            mvavg_xvals * fucci.TOT_LEN,
            mvperc_red[1],
            mvperc_red[-2],
            color="lightcoral",
            label="25th & 75th Percentiles",
        )

        mvavg_red = mvavg(pol_sort_centered_data1, WINDOW_FUCCI_PSEUDOTIME)
        mvavg_green = mvavg(pol_sort_centered_data0, WINDOW_FUCCI_PSEUDOTIME)
        plt.plot(
            mvavg_xvals * fucci.TOT_LEN, mvavg_red, color="r", label="Mean Intensity"
        )
        plt.plot(
            mvavg_xvals * fucci.TOT_LEN, mvavg_green, color="g", label="Mean Intensity"
        )
        plt.xlabel("Cell Cycle Time, hrs")
        plt.ylabel("Log10 Tagged CDT1 & GMNN Intensity")
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.tight_layout()
        plt.savefig("figures/FUCCIOverPseudotime.pdf")
        plt.savefig("figures/FUCCIOverPseudotime.png")
        plt.close()

        pd.DataFrame(
            {
                "cell_cycle_time_hrs": mvavg_xvals * fucci.TOT_LEN,
                "mvavgs_red": mvavg_red,
                "mvavgs_red_25p": mvperc_red[1],
                "mvavgs_red_75p": mvperc_red[-2],
                "mvavgs_green": mvavg_green,
                "mvavgs_green_25p": mvperc_green[1],
                "mvavgs_green_75p": mvperc_green[-2],
            }
        ).to_csv("output/FucciIntensitiesOnPseudotime.tsv", index=False, sep="\t")

    def sort_protein_results(self, protein_ccd):
        """Sort the protein data according to FUCCI pseudotime"""
        (
            self.pol_sort_well_plate,
            self.pol_sort_well_plate_imgnb,
            self.pol_sort_well_plate_imgnb_objnb,
            self.pol_sort_ab_nuc,
            self.pol_sort_ab_cyto,
            self.pol_sort_ab_cell,
            self.pol_sort_mt_cell,
            self.pol_sort_area_cell,
            self.pol_sort_area_nuc,
            self.pol_sort_fred,
            self.pol_sort_fgreen,
            self.pol_sort_mockbulk_phases,
        ) = pol_sort(
            self.pol_sort_inds,
            self.more_than_start,
            self.less_than_start,
            protein_ccd.well_plate,
            protein_ccd.well_plate_imgnb,
            protein_ccd.well_plate_imgnb_objnb,
            protein_ccd.ab_nuc,
            protein_ccd.ab_cyto,
            protein_ccd.ab_cell,
            protein_ccd.mt_cell,
            protein_ccd.area_cell,
            protein_ccd.area_nuc,
            protein_ccd.log_red_fucci_zeroc_rescale,
            protein_ccd.log_green_fucci_zeroc_rescale,
            protein_ccd.mockbulk_phases,
        )

    def save_protein_results(self):
        utils.np_save_overwriting(
            "output/pickles/pol_sort_well_plate.npy", self.pol_sort_well_plate
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_well_plate_imgnb.npy",
            self.pol_sort_well_plate_imgnb,
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_well_plate_imgnb_objnb.npy",
            self.pol_sort_well_plate_imgnb_objnb,
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_norm_rev.npy", self.pol_sort_norm_rev
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_ab_nuc.npy", self.pol_sort_ab_nuc
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_ab_cyto.npy", self.pol_sort_ab_cyto
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_ab_cell.npy", self.pol_sort_ab_cell
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_mt_cell.npy", self.pol_sort_mt_cell
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_area_cell.npy", self.pol_sort_area_cell
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_area_nuc.npy", self.pol_sort_area_nuc
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_fred.npy", self.pol_sort_fred
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_fgreen.npy", self.pol_sort_fgreen
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_centered_data0.npy", self.pol_sort_centered_data0
        )
        utils.np_save_overwriting(
            "output/pickles/pol_sort_centered_data1.npy", self.pol_sort_centered_data1
        )

    def pseudotime_protein(self, protein_ccd):
        """Generate a polar coordinate model of cell cycle progression based on the FUCCI intensities"""
        self.do_plotting = protein_ccd.do_plotting
        self.fucci_polar_coords(
            protein_ccd.fucci_data[:, 0], protein_ccd.fucci_data[:, 1], "Protein"
        )
        self.sort_protein_results(protein_ccd)
        self.save_protein_results()

        # Generate some plots
        if self.do_plotting:
            self.fucci_hist2d(
                "Protein",
                200,
                True,
                self.pol_sort_well_plate,
                self.pol_sort_ab_nuc,
                self.pol_sort_centered_data0,
                self.pol_sort_centered_data1,
            )
            self.plot_fucci_intensities_on_pseudotime()


# def pseudotime_rna(adata):
#     '''
#     Model the pseudotime of FUCCI markers measured by FACS for cell analyzed with single-cell RNA sequencing
#     Input: RNA-Seq data; cell cycle phase for each cell
#     Output: Cell cycle pseudotime for each cell; plots illustrating FUCCI intensities over pseudotime
#     '''
#     adataValidPhase = adata.obs["phase"] != "N/A"
#     utils.general_scatter_color(adata.obs["Green530"][adataValidPhase], adata.obs["Red585"][adataValidPhase], "log(Green530)", "log(Red585)",
#         adata.obs["phase"][adataValidPhase].apply(lambda x: get_phase_colormap()[x]), "Cell Cycle Phase", False, "", "figures/FucciPlotByPhase_RNA.png",
#         get_phase_colormap())

#     polar_coord_results = fucci_polar_coords(adata.obs["Green530"], adata.obs["Red585"], "RNA")
#     pol_sort_norm_rev, centered_data, pol_sort_centered_data0, pol_sort_centered_data1, pol_sort_phi, pol_sort_inds, pol_sort_inds_reorder, pol_sort_phi_reorder, more_than_start, less_than_start, start_pt, g1_end_pt, g1s_end_pt, cart_data_ur, R_2, start_phi = polar_coord_results

#     # Assign cells a pseudotime and visualize in fucci plot
#     pol_unsort = np.argsort(pol_sort_inds_reorder)
#     fucci_time = pol_sort_norm_rev[pol_unsort]
#     adata.obs["fucci_time"] = fucci_time

#     plt.figure(figsize=(6,5))
#     plt.scatter(adata.obs["Green530"], adata.obs["Red585"], c = adata.obs["fucci_time"], cmap="RdYlGn")
#     cbar = plt.colorbar()
#     cbar.set_label('Pseudotime',size=20)
#     cbar.ax.tick_params(labelsize=18)
#     plt.xlabel("log10(GMNN GFP Intensity)",size=20)
#     plt.ylabel("log10(CDT1 RFP Intensity)",size=20)
#     plt.tight_layout()
#     plt.savefig("figures/FucciAllFucciPseudotime.pdf")
#     plt.close()

#     # Save fucci times, so they can be used in other workbooks
#     fucci_time_inds = np.argsort(adata.obs["fucci_time"])
#     fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
#     cell_time_sort = pd.DataFrame({
#         "fucci_time" : fucci_time_sort,
#         "cell" : np.take(adata.obs["Well_Plate"], fucci_time_inds)})
#     cell_time_sort.to_csv("output/CellPseudotimes.csv", index=False)
#     cell_polar_coords = pd.DataFrame({
#         "cell" : np.take(adata.obs["Well_Plate"], fucci_time_inds),
#         "phase_by_facs_gating" : adata.obs["phase"][fucci_time_inds],
#         "raw_green530" : adata.obs["MeanGreen530"][fucci_time_inds],
#         "raw_red585" : adata.obs["MeanRed585"][fucci_time_inds],
#         "green530_lognorm_rescale" : adata.obs["Green530"][fucci_time_inds],
#         "red585_lognorm_rescale" : adata.obs["Red585"][fucci_time_inds],
#         "polar_coord_phi" : pol_sort_phi[pol_unsort][fucci_time_inds],
#         "fucci_time_hrs" : fucci_time_sort * fucci.TOT_LEN})
#     cell_polar_coords.to_csv("output/fucci_coords.csv", index=False)
#     pd.DataFrame({"fucci_time": fucci_time}).to_csv("output/fucci_time.csv", index=False)


def pseudotime_umap(adata, isIsoform=False):
    """
    Display FUCCI pseudotime on the UMAP created from the gene expression.
    Input: RNA-Seq data
    Output: UMAP diagram
    """
    nneighbors = [5, 10, 15, 30, 100]  # used nn=10 in the paper
    mindists = [0, 0.01, 0.05, 0.1, 0.5, 1]  # used 0.5 (default) in the paper
    for nn in nneighbors:
        for md in mindists:
            sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=40)
            sc.tl.umap(adata, min_dist=md)
            sc.pl.umap(adata, color=["fucci_time"], show=False, save=True)
            shutil.move(
                "figures/umap.pdf",
                f"figures/umapAllCellsSeqFucciPseudotime_nn{nn}_md{md}{'_Isoform' if isIsoform else ''}.pdf",
            )
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
