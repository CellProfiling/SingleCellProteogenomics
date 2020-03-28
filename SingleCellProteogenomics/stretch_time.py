#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 16:09:46 2018

@author: devinsullivan
"""

import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
from copy import deepcopy

def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i = round(i+step,14)

def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def histedges_equalA(x, nbin):
    pow = 0.5
    dx = np.diff(np.sort(x))
    tmp = np.cumsum(dx ** pow)
    tmp = np.pad(tmp, (1, 0), 'constant')
    return np.interp(np.linspace(0, tmp.max(), nbin + 1),
                     tmp,
                     np.sort(x))

def stretch_time(time_data,nbins=1000):
    #This function is supposed to create uniform density space

    n, bins, patches = plt.hist(time_data, histedges_equalN(time_data, nbins), normed=True)
    #data_hist = plt.hist(time_data,nbins)

    tmp_time_data = deepcopy(time_data)
#    ndecimals = np.ceil(np.log10(nbins))
#    rnd_time_data = np.round(time_data,decimals=int(ndecimals))

    trans_time = np.zeros([len(time_data)])
    for i,c_bin in enumerate(bins[1:]):
        #get curr bin indexs
        c_inds = np.argwhere(tmp_time_data<c_bin)
        trans_time[c_inds] = i/nbins
        tmp_time_data[c_inds] = np.inf
#        print(len(c_inds))

    return trans_time
