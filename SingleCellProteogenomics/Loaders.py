# -*- coding: utf-8 -*-
"""
Methods for loading results of previous analysis files. 

Separating analysis files improved the speed of testing and rerunning analyses.

@author: Anthony J. Cesnik, cesnik@stanford.edu
"""

import numpy as np
import pandas as pd
from SingleCellProteogenomics import utils

def load_protein_fucci_pseudotime():
    '''Loads results of previous analysis files for Protein FUCCI pseudotime analysis'''
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
       
        "wp_iscell" :  np.load("output/pickles/wp_iscell.npy", allow_pickle=True),
        "wp_isnuc" :  np.load("output/pickles/wp_isnuc.npy", allow_pickle=True),
        "wp_iscyto" :  np.load("output/pickles/wp_iscyto.npy", allow_pickle=True),
        
        "curr_wp_phases" : np.load("output/pickles/curr_wp_phases.npy", allow_pickle=True),
        "mockbulk_phases" : np.load("output/pickles/mockbulk_phases.npy", allow_pickle=True)}

def load_temporal_delay():
    '''Loads results of previous analysis files for temporal delay analysis'''
    return_dict = load_protein_fucci_pseudotime()
    add_dict = {"ccd_comp" : np.load("output/pickles/ccd_comp.npy", allow_pickle=True),
        "ccdtranscript" : np.load("output/pickles/ccdtranscript.npy", allow_pickle=True),
        "var_comp" : np.load("output/pickles/var_comp.npy", allow_pickle=True),
        "gini_comp" : np.load("output/pickles/gini_comp.npy", allow_pickle=True),
        "cv_comp" : np.load("output/pickles/cv_comp.npy", allow_pickle=True),
        "var_cell" : np.load("output/pickles/var_cell.npy", allow_pickle=True),
        "gini_cell" : np.load("output/pickles/gini_cell.npy", allow_pickle=True),
        "cv_cell" : np.load("output/pickles/cv_cell.npy", allow_pickle=True),
        "pol_sort_well_plate" :  np.load("output/pickles/pol_sort_well_plate.npy", allow_pickle=True),
        "pol_sort_norm_rev" :  np.load("output/pickles/pol_sort_norm_rev.npy", allow_pickle=True),
        "pol_sort_ab_nuc" :  np.load("output/pickles/pol_sort_ab_nuc.npy", allow_pickle=True),
        "pol_sort_ab_cyto" :  np.load("output/pickles/pol_sort_ab_cyto.npy", allow_pickle=True),
        "pol_sort_ab_cell" :  np.load("output/pickles/pol_sort_ab_cell.npy", allow_pickle=True),
        "pol_sort_mt_cell" :  np.load("output/pickles/pol_sort_mt_cell.npy", allow_pickle=True),
        "pol_sort_area_cell" : np.load("output/pickles/pol_sort_area_cell.npy", allow_pickle=True),
        "pol_sort_area_nuc" : np.load("output/pickles/pol_sort_area_nuc.npy", allow_pickle=True),
        "pol_sort_fred" :  np.load("output/pickles/pol_sort_fred.npy", allow_pickle=True),
        "pol_sort_fgreen" :  np.load("output/pickles/pol_sort_fgreen.npy", allow_pickle=True)}
    for item in add_dict.items():
        return_dict[item[0]] = item[1]
    return return_dict

def load_ptm_and_stability(adata):
    '''Loads results of previous analysis files for PTM and protein stability analyses'''
    return_dict = load_temporal_delay()

    # Get the labels in terms of ENSG
    ccdtranscript = return_dict["ccdtranscript"]
    ccd_comp = np.load("output/pickles/ccd_comp.npy", allow_pickle=True)
    nonccd_comp = np.load("output/pickles/nonccd_comp.npy", allow_pickle=True)
    ccdprotein_transcript_regulated = np.load("output/pickles/ccdprotein_transcript_regulated.npy", allow_pickle=True)
    ccdprotein_nontranscript_regulated = np.load("output/pickles/ccdprotein_nontranscript_regulated.npy", allow_pickle=True)
    genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])
    ccd_regev_filtered, ccd_filtered, nonccd_filtered = utils.ccd_gene_lists(adata)

    # Convert them to gene names and store them as such
    ccdtranscript_names = set(utils.ccd_gene_names(adata.var_names[ccdtranscript]))
    nonccdtranscript_names = set(utils.ccd_gene_names(adata.var_names[~ccdtranscript]))
    ccdprotein_transcript_regulated_names = set(utils.ccd_gene_names(adata.var_names[ccdprotein_transcript_regulated]))
    ccdprotein_nontranscript_regulated_names = set(utils.ccd_gene_names(adata.var_names[ccdprotein_nontranscript_regulated]))
    genes_analyzed_names = set(utils.ccd_gene_names(genes_analyzed))
    ccd_regev_filtered_names = set(utils.ccd_gene_names(ccd_regev_filtered))
    ccd_filtered_names = set(utils.ccd_gene_names(ccd_filtered))
    nonccdprotein_names = set(utils.ccd_gene_names(return_dict["wp_ensg"][nonccd_comp]))
    ccdprotein_names = set(utils.ccd_gene_names(np.unique(return_dict["wp_ensg"][ccd_comp]))) # just pseudotime

    add_dict = {"wp_max_pol" : np.load("output/pickles/wp_max_pol.npy", allow_pickle=True),
        "ccdtranscript_names" : ccdtranscript_names,
        "nonccdtranscript_names" : nonccdtranscript_names,
        "ccdprotein_transcript_regulated_names" : ccdprotein_transcript_regulated_names,
        "ccdprotein_nontranscript_regulated_names" : ccdprotein_nontranscript_regulated_names,
        "genes_analyzed_names" : genes_analyzed_names,
        "ccd_regev_filtered_names" : ccd_regev_filtered_names,
        "ccd_filtered_names" : ccd_filtered_names,
        "nonccdprotein_names" : nonccdprotein_names,
        "ccdprotein_names" : ccdprotein_names}
    for item in add_dict.items():
        return_dict[item[0]] = item[1]
    return return_dict