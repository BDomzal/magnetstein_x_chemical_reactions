import json

with open('../../config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]

import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *


experiment_name = 'Saccharose hydrolysis'

MIXTURE_PATH = config["mixture_paths"][experiment_name]
SUCROSE_TIME_AXIS_PATH = config["sucrose_time_axis_path"]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
REAGENTS_PATH = config["reagents_paths"][experiment_name]
REAGENTS_SEPARATORS = config["reagents_separators"][experiment_name]
SATURATED_COLORS_FOR_COMPONENTS = config["saturated_colors_for_components"][experiment_name]
SUBSTANCES_NAMES = config["substances_names_for_visualisations"][experiment_name]
ADDITIONAL_ANALYSIS_RESULTS_PATHS = config["additional_analysis_results_paths_incomplete_library"][experiment_name]


time_range = get_time_range(experiment_name, '../' + SUCROSE_TIME_AXIS_PATH)


mixture_time_data = load_mixture_time_data_v3(experiment_name, '../' + MIXTURE_PATH, MIXTURE_SEPARATOR)


reagents_spectra = load_and_preprocess_reagents_spectra(['../' + reag_path for reag_path in REAGENTS_PATH], REAGENTS_SEPARATORS)

kappa = 0.1
kappa_th = 1


for nr in [None, 0, 1, 2, 3]:

    which_reagent_to_remove = nr


    reagents_spectra = load_and_preprocess_reagents_spectra(['../' + reag_path for reag_path in REAGENTS_PATH], REAGENTS_SEPARATORS)
    reagents_spectra = reduce_library(reagents_spectra, which_reagent_to_remove)


    proportions_in_times = []
    noise_proportions_in_times = []
    noise = []
    noise_in_components = []

    estimation = estimate_proportions_in_time(np.array(mixture_time_data), reagents_spectra, 
                                            what_to_compare='area', 
                                            solver=pulp.GUROBI(msg=False),
                                            MTD=kappa, MTD_th=kappa_th,
                                             verbose=True)

    proportions_in_times = estimation['proportions_in_time']
    noise_proportions_in_times = estimation['proportion_of_noise_in_reagents_in_time']
    noise = estimation['noise_in_mixture_in_time']
    noise_in_components = estimation['noise_in_reagents_in_time']
    common_horizontal_axis = estimation['common_horizontal_axis_in_time']


    create_incomplete_library_figure(proportions_in_times, time_range, SUBSTANCES_NAMES, which_reagent_to_remove,
                                        SATURATED_COLORS_FOR_COMPONENTS, kappa, kappa_th, results_path=ADDITIONAL_ANALYSIS_RESULTS_PATHS)

