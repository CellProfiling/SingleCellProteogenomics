#%%
from imports import *

#%% Read in the plate counts
if not os.path.isfile("input/AllCountsForScanpy.csv") or not os.path.isfile("input/355CountsForScanpy.csv"):
    counts355 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_355\\counts.tab", delimiter="\t", index_col=0)
    counts356 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_356\\counts.tab", delimiter="\t", index_col=0)
    counts357 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_357\\counts.tab", delimiter="\t", index_col=0)
    counts355.columns += "_355"
    counts356.columns += "_356"
    counts357.columns += "_357"

if not os.path.isfile("input/AllRpkmForScanpy.csv") or not os.path.isfile("input/355RpkmForScanpy.csv"):
    rpkms355 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_355\\rpkms.tab", delimiter="\t", index_col=0)
    rpkms356 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_356\\rpkms.tab", delimiter="\t", index_col=0)
    rpkms357 = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\ESCG_data\\SS2_18_357\\rpkms.tab", delimiter="\t", index_col=0)
    rpkms355.columns += "_355"
    rpkms356.columns += "_356"
    rpkms357.columns += "_357"

#%% Read in all data; 
# keep symbols that aren't accepted to keep all data through analysis

# Symbol lookup
def gene_symbol_lookup():
    symbols = pd.read_csv("C:\\Users\\antho\\Box\\ProjectData\\CellCycle\\gene_symbols.synonyms.hgnc.txt", delimiter="\t")
    accepted_symbols = set(symbols["Approved symbol"])
    accepted_lookup = {}
    for idx, row in symbols.iterrows():
        for sss in str(row["Synonyms"]).split(", ") + str(row["Previous symbols"]).split(", "):
            accepted_lookup[sss] = row["Approved symbol"]
    return accepted_symbols, accepted_lookup

accepted_symbols, accepted_lookup = gene_symbol_lookup()
pd.DataFrame({"accepted_symbols" : list(accepted_symbols)}).to_pickle("output/accepted_gene_symbols.pkl")

def correct_symbol(x):
    if x in accepted_symbols: return x
    elif x in accepted_lookup: return accepted_lookup[x]
    elif x[-2] == "-": return x[:-2] # transcript symbols?
    else: return x

if not os.path.isfile("input/AllCountsForScanpy.csv") or not os.path.isfile("input/355CountsForScanpy.csv"):
    counts = pd.concat([counts355,counts356,counts357], axis=1, sort=False)
    countsT = counts.T
    countsT_acceptedish = countsT.rename(index=str, columns = dict([(x, correct_symbol(x)) for x in list(counts.index)]))
    print(f"{len([x for x in list(countsT.columns) if x not in accepted_symbols])}: number of genes that had no accepted symbol")
    print(f"{len([x for x in list(countsT_acceptedish.columns) if x not in accepted_symbols])}: number of genes that still have no accepted symbol")
    countsT_acceptedish.sort_index().to_csv("input/AllCountsForScanpy.csv")
    counts355.T.sort_index().to_csv("input/355CountsForScanpy.csv")
    counts356.T.sort_index().to_csv("input/356CountsForScanpy.csv")
    counts357.T.sort_index().to_csv("input/357CountsForScanpy.csv")

if not os.path.isfile("input/AllRpkmForScanpy.csv") or not os.path.isfile("input/355RpkmForScanpy.csv"):
    rpkms = pd.concat([rpkms355,rpkms356,rpkms357], axis=1, sort=False)
    rpkmsT = rpkms.T
    rpkmsT_acceptedish = rpkmsT.rename(index=str, columns = dict([(x, correct_symbol(x)) for x in list(rpkms.index)]))
    print(f"{len([x for x in list(rpkmsT.columns) if x not in accepted_symbols])}: number of genes that had no accepted symbol")
    print(f"{len([x for x in list(rpkmsT_acceptedish.columns) if x not in accepted_symbols])}: number of genes that still have no accepted symbol")
    rpkmsT_acceptedish.sort_index().to_csv("input/AllRpkmsForScanpy.csv")
    rpkms355.T.sort_index().to_csv("input/355RpkmsForScanpy.csv")
    rpkms356.T.sort_index().to_csv("input/356RpkmsForScanpy.csv")
    rpkms357.T.sort_index().to_csv("input/357RpkmsForScanpy.csv")

#%%
