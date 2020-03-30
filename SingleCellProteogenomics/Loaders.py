# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:41:32 2020

@author: antho
"""

import numpy as np

def load_protein_fucci_pseudotime():
    return {
        "u_plate" :  np.load("output/pickles/u_plate.npy", allow_pickle=True),
        "well_plate" : np.load("output/pickles/well_plate.npy", allow_pickle=True),
        "well_plate_imgnb" : np.load("output/pickles/well_plate_imgnb.npy", allow_pickle=True),
        "u_well_plates" : np.load("output/pickles/u_well_plates.npy", allow_pickle=True),
        "ab_nuc" : np.load("output/pickles/ab_nuc.npy", allow_pickle=True),
        "ab_cyto" : np.load("output/pickles/ab_cyto.npy", allow_pickle=True),
        "ab_cell" : np.load("output/pickles/ab_cell.npy", allow_pickle=True),
        "mt_cell" : np.load("output/pickles/mt_cell.npy", allow_pickle=True),
        "area_cell" : np.load("output/pickles/area_cell.npy", allow_pickle=True),
        "area_nuc" : np.load("output/pickles/area_nuc.npy", allow_pickle=True),
       
        "wp_ensg" :  np.load("output/pickles/wp_ensg.npy", allow_pickle=True), 
        "wp_ab" :  np.load("output/pickles/wp_ab.npy", allow_pickle=True), 
       
        "green_fucci" : np.load("output/pickles/green_fucci.npy", allow_pickle=True),
        "red_fucci" : np.load("output/pickles/red_fucci.npy", allow_pickle=True),
        "log_green_fucci_zeroc" : np.load("output/pickles/log_green_fucci_zeroc.npy", allow_pickle=True),
        "log_red_fucci_zeroc" : np.load("output/pickles/log_red_fucci_zeroc.npy", allow_pickle=True),
        "log_green_fucci_zeroc_rescale" : np.load("output/pickles/log_green_fucci_zeroc_rescale.npy", allow_pickle=True),
        "log_red_fucci_zeroc_rescale" : np.load("output/pickles/log_red_fucci_zeroc_rescale.npy", allow_pickle=True),
        "wp_comp_kruskal_gaussccd_adj" :  np.load("output/pickles/wp_comp_kruskal_gaussccd_adj.npy", allow_pickle=True),
        "wp_pass_kruskal_gaussccd_bh_comp" :  np.load("output/pickles/wp_pass_kruskal_gaussccd_bh_comp.npy", allow_pickle=True),
        "fucci_data" :  np.load("output/pickles/fucci_data.npy", allow_pickle=True),

        # "pol_sort_well_plate" :  np.load("output/pickles/pol_sort_well_plate.npy", allow_pickle=True),
        # "pol_sort_norm_rev" :  np.load("output/pickles/pol_sort_norm_rev.npy", allow_pickle=True),
        # "pol_sort_ab_nuc" :  np.load("output/pickles/pol_sort_ab_nuc.npy", allow_pickle=True),
        # "pol_sort_ab_cyto" :  np.load("output/pickles/pol_sort_ab_cyto.npy", allow_pickle=True),
        # "pol_sort_ab_cell" :  np.load("output/pickles/pol_sort_ab_cell.npy", allow_pickle=True),
        # "pol_sort_mt_cell" :  np.load("output/pickles/pol_sort_mt_cell.npy", allow_pickle=True),
        # "pol_sort_area_cell" : np.load("output/pickles/pol_sort_area_cell.npy", allow_pickle=True),
        # "pol_sort_area_nuc" : np.load("output/pickles/pol_sort_area_nuc.npy", allow_pickle=True),
        # "pol_sort_fred" :  np.load("output/pickles/pol_sort_fred.npy", allow_pickle=True),
        # "pol_sort_fgreen" :  np.load("output/pickles/pol_sort_fgreen.npy", allow_pickle=True),
       
        "wp_iscell" :  np.load("output/pickles/wp_iscell.npy", allow_pickle=True),
        "wp_isnuc" :  np.load("output/pickles/wp_isnuc.npy", allow_pickle=True),
        "wp_iscyto" :  np.load("output/pickles/wp_iscyto.npy", allow_pickle=True)}
     
