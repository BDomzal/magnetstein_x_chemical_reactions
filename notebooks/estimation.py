
experiment_name = 'Saccharose hydrolysis'
kappa = 0.1
kappa_th = 0.1

import json

with open('../config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]


import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *

MIXTURE_PATH = config["mixture_paths"][experiment_name]
REAGENTS_PATH = config["reagents_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
REAGENTS_SEPARATOR = config["reagents_separators"][experiment_name]
RESULT_PATH = config["results_paths"][experiment_name]


mixture_time_data = load_mixture_time_data_v3(experiment_name, MIXTURE_PATH, MIXTURE_SEPARATOR)

reagents_spectra = load_and_preprocess_reagents_spectra(REAGENTS_PATH, REAGENTS_SEPARATOR)


proportions_in_times, noise_proportions_in_times, common_horizontal_axis, noise, noise_in_components = run_estimation_in_time(
                                                                                                                            mixture_time_data, 
                                                                                                                            reagents_spectra, 
                                                                                                                            what_to_compare='area', 
                                                                                                                            MTD=kappa,
                                                                                                                            MTD_th=kappa_th, 
                                                                                                                            results_path=RESULT_PATH
                                                                                                                                )

