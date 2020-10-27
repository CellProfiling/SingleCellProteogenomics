# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 10:31:23 2020

@author: antho
"""

bioccd = np.genfromtxt("input/ProteinData/BiologicallyDefinedCCD.txt", dtype='str') # from mitotic structures
knownccd1 = np.genfromtxt("input/ProteinData/knownccd.txt", dtype='str') # from gene ontology, reactome, cyclebase 3.0, NCBI gene from mcm3
knownccd2 = np.genfromtxt("input/ProteinData/known_go_ccd.txt", dtype='str') # from GO cell cycle
knownccd3 = np.genfromtxt("input/ProteinData/known_go_proliferation.txt", dtype='str') # from GO proliferation
knownccd = np.concatenate((knownccd1, knownccd2, knownccd3))
ccdold = np.unique(np.concatenate((np.concatenate(np.asarray(pd.read_csv("c:/Users/antho/Desktop/knownCCD_old.csv"))), bioccd)))
oldNovel = ccdold[~np.isin(ccdold, knownccd)]
oldKnown = ccdold[np.isin(ccdold, knownccd)]
proteinresults = pd.read_csv("output/CellCycleVariationSummary.csv")

def compare(cutoff):
    newccd = np.unique(np.concatenate((proteinresults[proteinresults["mean_percvar_diff_from_random"] > cutoff]["ENSG"], bioccd)))
    print(f"COMPARISON MEANDIFF={cutoff}")
    print("NOVEL")
    newNovel = newccd[~np.isin(newccd, knownccd)]
    oldNovelOnly = ~np.isin(oldNovel, newNovel)
    bothNovel = np.isin(oldNovel, newNovel)
    newNovelOnly = ~np.isin(newNovel, oldNovel)
    print(f"{sum(oldNovelOnly)}: oldNovelOnly")
    print(f"{sum(bothNovel)}: bothNovel")
    print(f"{sum(newNovelOnly)}: newNovelOnly")
    print("KNOWN")
    newKnown = newccd[np.isin(newccd, knownccd)]
    oldKnownOnly = ~np.isin(oldKnown, newKnown)
    bothKnown = np.isin(oldKnown, newKnown)
    newKnownOnly = ~np.isin(newKnown, oldKnown)
    print(f"{sum(oldKnownOnly)}: oldKnownOnly")
    print(f"{sum(bothKnown)}: bothKnown")
    print(f"{sum(newKnownOnly)}: newKnownOnly")
    print("ALL")
    ccdTogether = np.unique(np.concatenate((ccdold, newccd)))
    oldCcdOnly = ~np.isin(ccdold, newccd)
    bothCcd = np.isin(ccdold, newccd)
    newCcdOnly = ~np.isin(newccd, ccdold)
    print(f"{sum(oldCcdOnly)}: oldCcdOnly")
    print(f"{sum(bothCcd)}: bothCcd")
    print(f"{sum(newCcdOnly)}: newCcdOnly")
    
    pd.DataFrame({"gene" : newNovel}).to_csv("output/novelCcdProteins{cutoff}.csv") #etc

compare(0.05)
compare(0.08)