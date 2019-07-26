from imports import *

def read_counts_and_phases(dd, count_or_rpkm, use_spike_ins, biotype_to_use):
    '''
    Read data into scanpy; Read phases and FACS intensities
    * dd: "All", "355", "356", "357"
    * count_or_rpkm: "Counts" or "Tpms"
    '''
    read_file = f"input/{count_or_rpkm}.csv" + (".ercc.csv" if use_spike_ins else "")
    if biotype_to_use != None and len(biotype_to_use) > 0:
        print(f"filtering for biotype: {biotype_to_use}")
        biotype_file = f"{read_file}.{biotype_to_use}.csv"
        if not os.path.exists(biotype_file):
            gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
            biotyped = gene_info[gene_info["biotype"] == biotype_to_use]["gene_id"]
            pd.read_csv(read_file)[biotyped ].to_csv(biotype_file)
        read_file = biotype_file

    adata = sc.read_csv(read_file)
    print(f"data shape: {adata.X.shape}")
    # adata.raw = adata

    phases = pd.read_csv("input/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")

    # Assign phases and log intensities; require log intensity
    adata.obs["phase"] = np.array(phases["Stage"])
    adata.obs["Green530"] = np.array(phases["Green530"])
    adata.obs["Red585"] = np.array(phases["Red585"])
    adata = adata[pd.notnull(adata.obs["Green530"]) & pd.notnull(adata.obs["Red585"])] # removes dark mitotic cells
    
    # Read in fucci pseudotime from previous analysis
    if os.path.isfile("output/fucci_time.csv"):
        adata.obs["fucci_time"] = np.array(pd.read_csv("output/fucci_time.csv")["fucci_time"])

    # Get info about the genes
    gene_info = pd.read_csv("input/IdsToNames.csv", header=None, names=["name", "biotype", "description"], index_col=0)
    adata.var["name"] = gene_info["name"]
    adata.var["biotype"] = gene_info["biotype"]
    adata.var["description"] = gene_info["description"]

    return adata, phases


def qc_filtering(adata, do_log_normalize):
    '''QC and filtering'''
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=100)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    if do_log_normalize: sc.pp.log1p(adata)
    phases = pd.read_csv("input/WellPlatePhasesLogNormIntensities.csv").sort_values(by="Well_Plate")
    phases_filt = phases[phases["Well_Plate"].isin(adata.obs_names)]
    phases_filt = phases_filt.reset_index(drop=True) # remove filtered indices
    print(f"data shape after filtering: {adata.X.shape}")
    return phases_filt


def ccd_gene_lists(adata):
    '''Read in the published CCD genes / Diana's CCD / Non-CCD genes'''
    gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    ccd_regev=pd.read_csv("input/ccd_regev.txt")    
    ccd=pd.read_csv("input/ccd_genes.txt")
    nonccd=pd.read_csv("input/nonccd_genes.txt")
    ccd_regev_filtered = list(gene_info[(gene_info["name"].isin(ccd_regev["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    ccd_filtered = list(gene_info[(gene_info["name"].isin(ccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    nonccd_filtered = list(gene_info[(gene_info["name"].isin(nonccd["gene"])) & (gene_info["gene_id"].isin(adata.var_names))]["gene_id"])
    return ccd_regev_filtered, ccd_filtered, nonccd_filtered

def ccd_gene_names(id_list_like):
    '''Convert gene ID list to gene name list'''
    gene_info = pd.read_csv("input/IdsToNames.csv", index_col=False, header=None, names=["gene_id", "name", "biotype", "description"])
    return gene_info[(gene_info["gene_id"].isin(id_list_like))]["name"]
