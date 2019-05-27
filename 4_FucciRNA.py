#%% [markdown]
# # Cell Cycle scRNA-Seq Analysis
# We've shown in single cell imaging data that there is variability that is correlated to the cell cycle, as well as a majority of proteins that vary outside of the cell cycle, which might be due to metabolic processes or other sources of variation.
# 
# Here, we've collected single-cell RNA-Seq (scRNA-Seq) data for these cells.
# * How many cells were analyzed?
# * How many reads per cell?
# * How many genes show variability in expression at the RNA level?
# * How many of the genes that are temporally regulated over the cell cycle, using the fucci colors to build the cell cycle trajectory?
# * How many of the genes that show variability not correlated to the cell cycle?

#%% Imports
from imports import *

#%% Read in RNA-Seq data again and the CCD gene lists
from methods_RNASeqData import read_counts_and_phases, qc_filtering, ccd_gene_lists

adata, phases_filt = read_counts_and_phases()
qc_filtering(adata)

ccd_regev_filtered, ccd_filtered, nonccd_filtered = ccd_gene_lists(adata)

#%% fucci pseudotime
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.diffmap(adata)

dd = "All"
plt.rcParams['figure.figsize'] = (10, 10)
sc.pl.diffmap(adata, color='fucci_time', projection="3d", show=True, save=True)
shutil.move("figures/diffmap.pdf", f"figures/diffmap{dd}CellsFucciPseudotime3d.pdf")
sc.pl.umap(adata, color=["fucci_time"], show=True, save=True)
shutil.move("figures/umap.pdf", f"figures/umap{dd}CellsSeqFucciPseudotime.pdf")

#%% Expression vs Pseudotime, uncomment to run again
def plot_expression_pseudotime(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        plt.scatter(adata.obs["fucci_time"], normalized_exp_data)
        plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
        plt.ylabel("RNA-Seq Counts, Normalized By Cell",size=36,fontname='Arial')
        plt.title(gene,size=36,fontname='Arial')
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

# plot_expression_pseudotime(ccd_regev_filtered, "figures/RegevGeneProfiles")
# plot_expression_pseudotime(ccd_filtered, "figures/DianaCcdGeneProfiles")
# plot_expression_pseudotime(nonccd_filtered, "figures/DianaNonCcdGeneProfiles")
# plot_expression_pseudotime(["ENO1"], "figures/OtherGeneProfiles")


#%% Expression Fucci, uncomment to run again
def plot_expression_facs(genelist, pppp, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        expression_data = np.exp(adata.X[:,list(adata.var_names).index(gene)]) - 1
        normalized_exp_data = expression_data / np.max(expression_data)
        plt.scatter(pppp["Green530"], pppp["Red585"], normalized_exp_data)
        plt.tight_layout()
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel("Log Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
        plt.title(gene,size=20,fontname='Arial')
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

phasesFilt = phases_filt[pd.notnull(phases_filt.Green530) & pd.notnull(phases_filt.Red585)] # stage may be null
phasesfiltfilt = phasesFilt[phasesFilt["Well_Plate"].isin(adata.obs_names)]

# plot_expression_facs(ccd_regev_filtered, phasesfiltfilt, "figures/RegevGeneFucci")
# plot_expression_facs(ccd_filtered, phasesfiltfilt, "figures/DianaCcdGeneFucci")
# plot_expression_facs(nonccd_filtered, phasesfiltfilt, "figures/DianaNonCcdGeneFucci")

#%% Expression-FUCCI facs of anillin by plate
def plot_expression_facs_plate(genelist, pppp, plate, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    plt.subplot(1, 3)
    for tt in ["355", "356", "357"]:
        for gene in genelist:
            pppp[gene] = adata.X[:,list(adata.var_names).index(gene)]
            pppplate = pppp[pd.notnull(pppp.Green530) & pd.notnull(pppp.Red585) & pd.notnull(pppp.Stage) & phasesfiltfilt.Well_Plate.str.endswith(plate)]
            plt.scatter(pppplate["Green530"], pppplate["Red585"], c = pppplate[gene])
            plt.tight_layout()
            cbar = plt.colorbar()
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel("Log Normalized RNA-Seq Counts", rotation=270,size=16,fontname='Arial')
            plt.title(gene,size=20,fontname='Arial')
            plt.savefig(f"{outfolder}/{gene}Plate{plate}.png")
            plt.close()
            pppp.drop(gene, 1)

    plot_expression_facs_plate(["ANLN"], phasesfiltfilt, tt, f"figures/FucciPlotByPlates")

        
#%% Moving avg Expression vs Pseudotime, uncomment to run again
def moving_average(a, n):
    '''A formula for the moving average of an array (a) with window size (n)'''
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

expression_data = np.exp(adata.X) - 1
normalized_exp_data = (expression_data.T / np.max(expression_data, axis=0)[:,None]).T
fucci_time_inds = np.argsort(adata.obs["fucci_time"])
fucci_time_sort = np.take(np.array(adata.obs["fucci_time"]), fucci_time_inds)
norm_exp_sort = np.take(normalized_exp_data, fucci_time_inds, axis=0)

def plot_expression_avg_pseudotime(genelist, outfolder):
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    for gene in genelist:
        nexp = norm_exp_sort[:,list(adata.var_names).index(gene)]
        df = pd.DataFrame({"fucci_time" : fucci_time_sort, gene : nexp})
        # plt.scatter(df["fucci_time"], df[gene], label="Normalized Expression")
        bin_size = 100
        plt.plot(df["fucci_time"], 
            df[gene].rolling(bin_size).mean(), 
            color="blue", 
            label=f"Moving Average by {bin_size} Cells")
        plt.fill_between(df["fucci_time"], 
            df[gene].rolling(bin_size).quantile(0.10),
            df[gene].rolling(bin_size).quantile(0.90), 
            color="lightsteelblue", 
            label="10th & 90th Percentiles")
        # plt.plot(df["fucci_time"], color="orange", label="Normalized Expression, 10th Percentile")
        # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() + 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
        # plt.plot(df["fucci_time"], df[gene].rolling(bin_size).mean() - 2 * df[gene].rolling(bin_size).std(), color="purple", label="Normalized Expression, 95% CI")
        plt.xlabel("Fucci Pseudotime",size=36,fontname='Arial')
        plt.ylabel("RNA-Seq Counts, Normalized By Cell",size=36,fontname='Arial')
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.title(gene,size=36,fontname='Arial')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{outfolder}/{gene}.png")
        plt.close()

plot_expression_avg_pseudotime(ccd_regev_filtered, "figures/RegevGeneProfilesAvgs")
# plot_expression_avg_pseudotime(ccd_filtered, "figures/DianaCcdGeneProfilesAvgs")
# plot_expression_avg_pseudotime(nonccd_filtered, "figures/DianaNonCcdGeneProfilesAvgs")


#%%
