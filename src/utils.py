import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from masserstein import NMRSpectrum, estimate_proportions, estimate_proportions_in_time

def load_spectrum(mixture_time_data, moment_of_time):
    ppm = mixture_time_data['ppm']
    intensity = mixture_time_data['t'+str(moment_of_time)]
    sp = NMRSpectrum(confs = list(zip(ppm, intensity)))
    return sp

def load_mixture_time_data(mixture_path, mixture_separator, last_columns_to_skip=2):

    mixture_time_data = pd.read_csv(mixture_path, sep = mixture_separator).iloc[:,:-last_columns_to_skip]
    names = ['ppm'] + ['t' + str(nb) for nb in range(1, mixture_time_data.shape[1])]
    mixture_time_data.columns = names

    return mixture_time_data

def load_mixture_time_data_v2(mixture_path, mixture_separator, last_columns_to_skip=2):

    mixture_time_data = pd.read_csv(mixture_path, sep = mixture_separator).iloc[:,:-last_columns_to_skip]
    mixture_time_data.dropna(inplace=True)
    names = ['ppm'] + ['t' + str(nb) for nb in range(1, mixture_time_data.shape[1])]
    mixture_time_data.columns = names

    return mixture_time_data


def create_reagent_spectrum(mixture_time_data, mixture_to_extract_reagents_spectra, reagent_intervals):

    reagents_sp = []

    mix_to_cut = pd.DataFrame(
                            data = load_spectrum(mixture_time_data, mixture_to_extract_reagents_spectra).confs,
                            columns = ['ppm', 't' + str(mixture_to_extract_reagents_spectra)]
                            )

    reagent = mix_to_cut[
                        mix_to_cut['ppm'].apply(
                                            lambda x: np.any([x>reagent_intervals[i][0] and x<reagent_intervals[i][1] for i in range(len(reagent_intervals))])
                                            )
                        ]

    ppm = reagent['ppm'].values
    ints = reagent['t' + str(mixture_to_extract_reagents_spectra)].values
    return NMRSpectrum(confs = list(zip(ppm, ints)))

def create_reagent_spectrum_v2(mixture_time_data, mixture_to_extract_reagents_spectra, reagent_intervals):

    reagents_sp = []

    mix_to_cut = pd.DataFrame(
                            data = load_spectrum(mixture_time_data, mixture_to_extract_reagents_spectra).confs,
                            columns = ['ppm', 't' + str(mixture_to_extract_reagents_spectra)]
                            )

    reagent = mix_to_cut[
                        mix_to_cut['ppm'].apply(
                                            lambda x: np.any([x>=reagent_intervals[i][0] and x<=reagent_intervals[i][1] for i in range(len(reagent_intervals))])
                                            )
                        ]

    return reagent

def cut_mixture_time_data_to_regions(mixture_time_data, intervals):

    res = mixture_time_data[mixture_time_data['ppm'].apply(
                                                            lambda x: np.any([x>intervals[i][0] and x<intervals[i][1] for i in range(len(intervals))])
                                                            )
                            ]
    return res.values

def cut_mixture_time_data_to_regions_v2(mixture_time_data, intervals):

    res = mixture_time_data[mixture_time_data['ppm'].apply(
                                                            lambda x: np.any([x>=intervals[i][0] and x<=intervals[i][1] for i in range(len(intervals))])
                                                            )
                            ]
    return res.values

def cut_ppm_axis(mixture_time_data, full_ppm_interval):

    res = mixture_time_data[mixture_time_data['ppm'].apply(
                                                            lambda x: x>full_ppm_interval[0] and x<full_ppm_interval[1]
                                                            )
                            ]
    return res

def preprocess_mnova_integrals(integrals_path, raw_integrals_separators, substances_names, skiprows=2, columns_to_use=[2,4,6,8]):
    
    integrals = pd.read_csv(integrals_path, sep = raw_integrals_separators, skiprows=skiprows).iloc[:,columns_to_use]
    integrals.fillna(0., inplace=True)
    integrals.columns = substances_names
    integrals[integrals < 0] = 0
    return integrals