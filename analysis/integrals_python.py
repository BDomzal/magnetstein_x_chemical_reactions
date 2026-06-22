# ### Imports

experiment_name = 'Saccharose hydrolysis'

import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import pulp
import pickle
import time

import json

with open('../config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]


import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *

MIXTURE_PATH_INT = config["mixture_paths_int"][experiment_name]
#(it"s ok that we are using preprocessed_mixture for PMG 287 
#cause the only preprocessing here is zeros from edges i.e. taking region: (2.952301, 10.387961))
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PYTHON_INTEGRALS_PATH = config["python_integrals_paths"][experiment_name]
SUBSTANCES_NAMES = config["substances_names"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]


import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
from magnetstein import NMRSpectrum, estimate_proportions


# ### Data

# #### Mixture in time


mixture_time_data = load_mixture_time_data_v4(experiment_name, MIXTURE_PATH_INT, MIXTURE_SEPARATOR)


# Note that no baseline correction is needed for saccharose hydrolysis anymore!

# ### Integrals changing in time

data_cut_to_intervals = cut_mixture_time_data_to_regions_v3(experiment_name, mixture_time_data, INTEGRATION_INTERVALS)

python_integrals = compute_integrals_in_python(experiment_name, mixture_time_data, data_cut_to_intervals)

# ### Figures

for i in range(python_integrals.shape[1]-1):
    plt.plot(python_integrals[:,i] / python_integrals[:, :-1].sum(1), 'p')
plt.show()


# ### Saving results

save_python_integrals(python_integrals, INTEGRATION_INTERVALS, SUBSTANCES_NAMES, PYTHON_INTEGRALS_PATH, experiment_name)