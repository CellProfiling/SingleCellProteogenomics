#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates a uniform density across pseudotime for polar coordinate calculations.

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
        i = round(i + step, 14)

def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x))

def histedges_equalA(x, nbin):
    pow = 0.5
    dx = np.diff(np.sort(x))
    tmp = np.cumsum(dx ** pow)
    tmp = np.pad(tmp, (1, 0), 'constant')
    return np.interp(np.linspace(0, tmp.max(), nbin + 1), tmp, np.sort(x))

def stretch_time(time_data,nbins=1000):
    '''This function is supposed to create uniform density space'''
    n, bins, patches = plt.hist(time_data, histedges_equalN(time_data, nbins), density=True)
    tmp_time_data = deepcopy(time_data)
    trans_time = np.zeros([len(time_data)])
    # Get bin indexes
    for i,c_bin in enumerate(bins[1:]):
        c_inds = np.argwhere(tmp_time_data<c_bin)
        trans_time[c_inds] = i/nbins
        tmp_time_data[c_inds] = np.inf
    return trans_time
