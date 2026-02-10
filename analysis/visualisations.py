experiment_name = 'PMG 287 monitoring with full mixture spectrum'

import json

with open('../config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]

import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *
from visualisation_utils import *


SUCROSE_TIME_AXIS_PATH = config["sucrose_time_axis_path"]
MIXTURE_PATH = config["mixture_paths"][experiment_name]
REAGENTS_PATH = config["reagents_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
REAGENTS_SEPARATOR = config["reagents_separators"][experiment_name]
RESULT_PATH = config["results_paths"][experiment_name]
SUBSTANCES_NAMES = config["substances_names_for_visualisations"][experiment_name]
INTEGRALS_TO_COMPARE_WITH = config["integrals_to_compare_with"][experiment_name]
INTEGRALS_SEPARATORS = config["integrals_separators"][experiment_name]
PYTHON_INTEGRALS_TO_COMPARE_WITH = config["python_integrals_to_compare_with"][experiment_name]
PYTHON_INTEGRALS_SEPARATORS = config["python_integrals_separators"][experiment_name]
COLORS_FOR_COMPONENTS = config["colors_for_components"][experiment_name]
SATURATED_COLORS_FOR_COMPONENTS = config["saturated_colors_for_components"][experiment_name]


all_kappas, all_kappas_th = get_all_kappa_values(RESULT_PATH)
all_kappas, all_kappas_th = remove_kappa_values(all_kappas, all_kappas_th, experiment_name)



results_dict = create_results_dict(RESULT_PATH, all_kappas, all_kappas_th)


# ### Loading integrals computed in MNova


mnova_integrals_proportions = load_mnova_integrals(INTEGRALS_TO_COMPARE_WITH, INTEGRALS_SEPARATORS, experiment_name)


# ### Loading integrals computed in Python

python_integrals_proportions = load_python_integrals(PYTHON_INTEGRALS_TO_COMPARE_WITH, PYTHON_INTEGRALS_SEPARATORS, experiment_name)


# ### Setting time for particular experiments

time_range = get_time_range(experiment_name, SUCROSE_TIME_AXIS_PATH)

time_range_integrals_python = get_time_range_integrals_python(experiment_name, SUCROSE_TIME_AXIS_PATH)

time_range_integrals_mnova = get_time_range_integrals_mnova(experiment_name, SUCROSE_TIME_AXIS_PATH)


# ### Visualising results for different parameters, different components

component_nr = 1

plot_proportion_against_time_single_component(
                                            experiment_name,
                                            component_nr, 
                                            all_kappas, 
                                            all_kappas_th, 
                                            results_dict,
                                            SUBSTANCES_NAMES, 
                                            time_range,
                                            SATURATED_COLORS_FOR_COMPONENTS,
                                            path_to_save=None
                                            )


# ### Visualising results for different parameters, all components together

plot_proportion_against_time_all_components(    
                                            experiment_name,
                                            all_kappas, 
                                            all_kappas_th, 
                                            results_dict,
                                            SUBSTANCES_NAMES,
                                            time_range,
                                            SATURATED_COLORS_FOR_COMPONENTS,
                                            path_to_save=None
                                                )


# ### Visualising results for different parameters, all components together + integrals

plot_proportion_against_time_all_components_plus_integrals(
                                                                experiment_name,
                                                                all_kappas, 
                                                                all_kappas_th, 
                                                                results_dict,
                                                                SUBSTANCES_NAMES,
                                                                time_range,
                                                                SATURATED_COLORS_FOR_COMPONENTS,
                                                                COLORS_FOR_COMPONENTS,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                python_or_mnova = 'Mnova',
                                                                path_to_save=None
                                                                )


# ### Single plot magnetstein + integral (chosen values of parameters)

plot_proportion_and_integrals_against_time_chosen_kappas(
                                                        experiment_name,
                                                        RESULT_PATH,
                                                        0.5,
                                                        0.5,
                                                        SUBSTANCES_NAMES,
                                                        time_range,
                                                        SATURATED_COLORS_FOR_COMPONENTS,
                                                        COLORS_FOR_COMPONENTS,
                                                        mnova_integrals_proportions,
                                                        python_integrals_proportions,
                                                        time_range_integrals_mnova,
                                                        time_range_integrals_python,
                                                        python_or_mnova='Mnova',
                                                        include_ticks=True,
                                                        path_to_save=None)


# ### Single plot for magnetstein (for chosen values of parameters)

plot_proportion_against_time_chosen_kappas(
                                            experiment_name,
                                            RESULT_PATH,
                                            0.5,
                                            0.5,
                                            SUBSTANCES_NAMES,
                                            time_range,
                                            SATURATED_COLORS_FOR_COMPONENTS,
                                            include_ticks=True,
                                            path_to_save=None
                                            )


# ### Single plot magnetstein + MNova + Python (chosen values of parameters)

plot_proportion_and_all_integrals_against_time_chosen_kappas(
                                                                experiment_name,
                                                                RESULT_PATH,
                                                                0.5,
                                                                0.5,
                                                                SUBSTANCES_NAMES,
                                                                time_range,
                                                                SATURATED_COLORS_FOR_COMPONENTS,
                                                                COLORS_FOR_COMPONENTS,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                path_to_save=None
                                                                )

plot_integrals_against_time_chosen_kappas(
                                                experiment_name,
                                                SUBSTANCES_NAMES,
                                                COLORS_FOR_COMPONENTS,
                                                mnova_integrals_proportions,
                                                python_integrals_proportions,
                                                time_range_integrals_mnova,
                                                time_range_integrals_python,
                                                python_or_mnova='Python',
                                                path_to_save=None
                                                )


plot_proportion_against_time_all_components_added_plus_integrals(
                                                                experiment_name,
                                                                [0],
                                                                [1,2],
                                                                all_kappas, 
                                                                all_kappas_th, 
                                                                results_dict,
                                                                SUBSTANCES_NAMES,
                                                                time_range,
                                                                SATURATED_COLORS_FOR_COMPONENTS,
                                                                COLORS_FOR_COMPONENTS,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                python_or_mnova = 'Mnova',
                                                                path_to_save=None
                                                                )


components_numbers = [0]

plot_proportion_against_time_chosen_components_added(
                                                        experiment_name,
                                                        [0],
                                                        all_kappas, 
                                                        all_kappas_th, 
                                                        results_dict,
                                                        SUBSTANCES_NAMES,
                                                        time_range,
                                                        path_to_save=None
                                                        )


best_kappa = 0.5
best_kappa_th = 0.5

plot_proportion_against_time_single_component_chosen_kappa(
                                                                experiment_name,
                                                                0,
                                                                RESULT_PATH,
                                                                0.5,
                                                                0.5,
                                                                SUBSTANCES_NAMES,
                                                                time_range,
                                                                SATURATED_COLORS_FOR_COMPONENTS,
                                                                path_to_save=None
                                                                )



plot_proportion_against_time_components_added_together_chosen_kappa(
                                                                        experiment_name,
                                                                        [0,2],
                                                                        RESULT_PATH,
                                                                        0.5,
                                                                        0.5,
                                                                        SUBSTANCES_NAMES,
                                                                        time_range,
                                                                        path_to_save=None
                                                                        )



plot_noise_proportion_chosen_kappas(
                                        experiment_name,
                                        RESULT_PATH,
                                        0.5,
                                        0.5,
                                        time_range,
                                        path_to_save=None
                                        )



plot_noise_proportion_chosen_kappas_v2(
                                        experiment_name,
                                        RESULT_PATH,
                                        0.5,
                                        0.5,
                                        time_range,
                                        path_to_save=None
                                        )


plot_added_components_plus_integrals_against_time_chosen_kappas(
                                                                experiment_name,
                                                                [0],
                                                                [1,2],
                                                                RESULT_PATH,
                                                                0.5,
                                                                0.5,
                                                                time_range,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                colors_list = ['grey', 'purple'],
                                                                include_integrals = True,
                                                                path_to_save=None
                                                                )



plot_noise_proportion_in_components_chosen_kappas(
                                        experiment_name,
                                        RESULT_PATH,
                                        0.5,
                                        0.5,
                                        time_range,
                                        path_to_save=None
                                        )



plot_noise_proportion_against_time (
                                        experiment_name,
                                        all_kappas, 
                                        all_kappas_th, 
                                        SUBSTANCES_NAMES,
                                        results_dict,
                                        time_range,
                                        path_to_save=None
                                        )


visualise_signal_removed_from_mixture(
                                            experiment_name,
                                            0.5,
                                            0.5,
                                            RESULT_PATH,
                                            MIXTURE_PATH,
                                            MIXTURE_SEPARATOR,
                                            moment_of_time = 1,
                                            path_to_save=None
                                            )



visualise_signal_removed_from_mixture_plus_reagents(
                                                    experiment_name,
                                                    0.5,
                                                    0.5,
                                                    RESULT_PATH,
                                                    MIXTURE_PATH,
                                                    MIXTURE_SEPARATOR,
                                                    REAGENTS_PATH,
                                                    REAGENTS_SEPARATOR,
                                                    moment_of_time = 1,
                                                    path_to_save=None
                                                    )


visualise_signal_removed_from_components(
                                            experiment_name,
                                            [0,1],
                                            0.5,
                                            0.5,
                                            RESULT_PATH,
                                            REAGENTS_PATH,
                                            REAGENTS_SEPARATOR,
                                            SATURATED_COLORS_FOR_COMPONENTS,
                                            moment_of_time = 60,
                                            path_to_save=None
                                                    )



plot_mixture_and_reagents(
                                experiment_name,
                                RESULT_PATH,
                                MIXTURE_PATH,
                                MIXTURE_SEPARATOR,
                                REAGENTS_PATH,
                                REAGENTS_SEPARATOR,
                                moment_of_time = 1,
                                path_to_save=None
                                )


