# get_ipython().system('pip install scanpy pandas')
# get_ipython().system('conda install -y -c vtraag python-igraph')
# get_ipython().system('conda install -y -c conda-forge statsmodel')

# Render our plots inline
# get_ipython().run_line_magic('matplotlib', 'inline')

#%% some setup
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpcolors
import matplotlib.patches as mpatches
import scanpy as sc
import os
import shutil
import scipy
import scipy.stats
