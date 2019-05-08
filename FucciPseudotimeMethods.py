import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import operator

from sklearn.neighbors import RadiusNeighborsRegressor
from scipy.optimize import least_squares
from scipy.optimize import minimize_scalar

###CONSTANTS
NBINS = 150 #number of bins, arbitrary choice for now
DO_PLOTS = True #flag of whether to plot each well and save the plot
PSIN_INIT = [np.nan,1,1,np.nan]
WINDOW = 20 #Number of points for moving average window, arbitrary choice
PERCVAR_CUT = 0.1 #Min % of the variance explained by the cell cycle.
FDR_CUT = 0.05 #False discovery rate we will tolerate

##MT CONSTANTS -- Calculated by running in mt mode, here we use cytosol
MT_GINI = 0.06600124369521103
MT_PERC_VAR = 0.09086819
MT_FDR = 0.3783993#FROM DIANA: cyto 0.3783993  and cell 0.0756249

