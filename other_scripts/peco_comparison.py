# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 14:38:14 2020

@author: antho
"""

from SingleCellProteogenomics.utils import *
from SingleCellProteogenomics import utils
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable
import math
from sklearn.linear_model import LinearRegression

peco=pd.read_csv("C:/Users/antho/Dropbox/Projects/CellCycle/pseudotime.csv")
print(f"{list(peco["Unnamed: 0"]) == list(adata.obs["Well_Plate"])}: cells are still parallel")

# plot with geminin log TPM values as color
xx = adata.obs["fucci_time"]
yy = peco["cellcycle_peco"] / (2 * math.pi)
plt.scatter(xx, yy, c=adata.X[:,list(adata.var_names).index("ENSG00000112312")], alpha=0.5)
plt.show(); plt.close()

# plot with geminin intensity as color
plt.scatter(xx, yy, c=adata.obs["Green530"], alpha=0.5)
plt.xlabel("FUCCI Pseudotime")
plt.ylabel("Predicted Pseudotime from scRNA-Seq using peco")
cbar = plt.colorbar()
cbar.set_label('Log10 Green530 Intensity',fontname='Arial',size=12)
plt.savefig("figures/PecoPredictedPseudotime.png")
plt.show(); plt.close()

reg = LinearRegression().fit(np.hstack((xx[:-20], yy[:-20])), np.hstack((xx[-20:], yy[-20:])))
pred = reg.predict(xx[-20:], yy[-20:])

print(f"{np.array(adata.obs['Well_Plate'])[(np.array(xx) < 0.2) & (np.array(yy) > 0.8)][:5]}: top left")
print(f"{np.array(adata.obs['Well_Plate'])[(np.array(xx) < 0.2) & (np.array(yy) < 0.2)][:5]}: bottom left")
print(f"{np.array(adata.obs['Well_Plate'])[(np.array(xx) > 0.8) & (np.array(yy) > 0.8)][:5]}: top right")