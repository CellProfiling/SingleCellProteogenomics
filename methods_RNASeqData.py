from imports import *

def read_counts_and_phases():
    '''Read data into scanpy; Read phases and FACS intensities'''
    dd = "All" # can also be 355,356, or 357
    adata = sc.read_csv(f"input/{dd}CountsForScanpy.csv")
    adata.var_names_make_unique()
    # adata.raw = adata

    phases = pd.read_csv("input/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")
    phases_filt = phases[phases["Well_Plate"].isin(adata.obs_names)]
    phases_filt = phases_filt.reset_index(drop=True) # remove filtered indices

    # Assign phases and log intensities; require log intensity
    adata.obs["phase"] = np.array(phases_filt["Stage"])
    adata.obs["phase_ajc"] = np.array(phases_filt["StageAJC"])
    adata.obs["Green530"] = np.array(phases_filt["Green530"])
    adata.obs["Red585"] = np.array(phases_filt["Red585"])
    adata = adata[pd.notnull(adata.obs["Green530"]) & pd.notnull(adata.obs["Red585"])] # removes dark mitotic cells

    # Read in fucci pseudotime from previous analysis
    if os.path.isfile("output/fucci_time.csv"):
        adata.obs["fucci_time"] = np.array(pd.read_csv("output/fucci_time.csv")["fucci_time"])

    return adata, phases_filt

def qc_filtering(adata):
    '''QC and filtering'''
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=100)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

def ccd_gene_lists(adata):
    '''Read in the published CCD genes / Diana's CCD / Non-CCD genes'''
    ccd_regev=pd.read_csv("input/ccd_regev.txt")
    ccd_regev_filtered = [gene for gene in ccd_regev["gene"] if gene in adata.var_names]
    ccd=pd.read_csv("input/ccd_genes.txt")
    nonccd=pd.read_csv("input/nonccd_genes.txt")
    ccd_filtered = [gene for gene in ccd["gene"] if gene in adata.var_names]
    nonccd_filtered = [gene for gene in nonccd["gene"] if gene in adata.var_names]
    return ccd_regev_filtered, ccd_filtered, nonccd_filtered