# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:20:24 2020

@author: antho
"""

from utils import *
import utils
import numpy as np
import sklearn.mixture
import seaborn as sbn
plt.rcParams['pdf.fonttype'], plt.rcParams['ps.fonttype'], plt.rcParams['savefig.dpi'] = 42, 42, 300 #Make PDF text readable

NBINS = 150 #number of bins, arbitrary choice for now

def calc_R(xc, yc, x, y):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c,x,y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    print(c)
    Ri = calc_R(c[0],c[1],x,y)
    return Ri - Ri.mean()

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol_sort(inds, nuc, cyto, cell, mt):
    '''Sort data by polar coordinates'''
    return nuc[inds], cyto[inds], cell[inds], mt[inds]
