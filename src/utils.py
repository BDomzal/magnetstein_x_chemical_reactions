import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from masserstein import NMRSpectrum, estimate_proportions, estimate_proportions_in_time
import pulp
import time
from textwrap import wrap
from scipy.stats import cauchy
import math

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

def get_time_range(experiment_name, sucrose_time_axis_path):

    if experiment_name.startswith('Saccharose hydrolysis'):
        return (1/60)*np.loadtxt(sucrose_time_axis_path)[::10]
    elif experiment_name.startswith('PMG 284 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1000,10))]
    elif experiment_name.startswith('PMG 287 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1000,10))]
    else:
        return

def get_time_range_integrals_python(experiment_name, sucrose_time_axis_path):

    if experiment_name.startswith('Saccharose hydrolysis'):
        return (1/60)*np.loadtxt(sucrose_time_axis_path)
    elif experiment_name.startswith('PMG 284 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1000))]
    elif experiment_name.startswith('PMG 287 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1000))]
    else:
        return

def get_time_range_integrals_mnova(experiment_name, sucrose_time_axis_path):
    if experiment_name.startswith('Saccharose hydrolysis'):
        return (1/60)*np.loadtxt(sucrose_time_axis_path)[:-1]
    elif experiment_name.startswith('PMG 284 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1001))]
    elif experiment_name.startswith('PMG 287 monitoring'):
        return [(1/60)*5*(el-1) for el in list(range(1,1001))]
    else:
        return

def load_mixture_time_data_v3(experiment_name, mixture_path, mixture_separator):

    mixture_time_data = pd.read_csv(mixture_path, sep = mixture_separator)

    if experiment_name.startswith('Saccharose hydrolysis'):
        ppm = mixture_time_data.iloc[:,0:1]
        every_10th = mixture_time_data.iloc[:,1:].iloc[:,::10]
        mixture_time_data = pd.concat((ppm, every_10th), axis=1)

    elif experiment_name.startswith('PMG 284 monitoring'):
        ppm = mixture_time_data.iloc[:,:-1].iloc[:,0:1]
        every_10th = mixture_time_data.iloc[:,:-1].iloc[:,1:].iloc[:,::10]
        mixture_time_data = pd.concat((ppm, every_10th), axis=1)

    elif experiment_name.startswith('PMG 287 monitoring'):
        ppm = mixture_time_data.iloc[:,:-1].iloc[:,0:1]
        every_10th = mixture_time_data.iloc[:,:-1].iloc[:,1:].iloc[:,::10]
        mixture_time_data = pd.concat((ppm, every_10th), axis=1)

    names = ['ppm'] + ['t' + str(nb) for nb in range(1, mixture_time_data.shape[1])]
    mixture_time_data.columns = names
    return mixture_time_data

def load_and_preprocess_reagents_spectra(reagents_path, reagents_separator):

    reagents_spectra = []
    for reagent in reagents_path:
        reag = pd.read_csv(reagent, sep=reagents_separator, header=None).iloc[:,:2]
        reagents_spectra.append(reag)

    reagents_spectra2 = []
    for reag in reagents_spectra:
        ppm = reag.iloc[:,0]
        ints = reag.iloc[:,1]
        sp = NMRSpectrum(confs = list(zip(ppm, ints)))
        reagents_spectra2.append(sp)

    reagents_spectra = reagents_spectra2

    for sp in reagents_spectra:
        sp.trim_negative_intensities()
        sp.normalize()

    return reagents_spectra

def reduce_library(reagents_spectra, which_reagent_to_remove):
    if which_reagent_to_remove is None:
        return reagents_spectra
    else:
        res = reagents_spectra[:which_reagent_to_remove] + reagents_spectra[(which_reagent_to_remove + 1):]
        return res

def get_names_and_colors(substances_names, saturated_colors_for_components, which_reagent_to_remove=None):

    color_list = saturated_colors_for_components

    if which_reagent_to_remove is None:
        return substances_names, color_list
    else:
        color_list = color_list[:which_reagent_to_remove] + color_list[(which_reagent_to_remove + 1):]
        names = substances_names[:which_reagent_to_remove] + substances_names[(which_reagent_to_remove + 1):]

    return names, color_list

def create_incomplete_library_figure(proportions_in_times, time_range, substances_names, which_reagent_to_remove,
                                    saturated_colors_for_components, kappa, kappa_th, results_path=None):


    names, color_list= get_names_and_colors(
                                           substances_names,
                                           saturated_colors_for_components,
                                           which_reagent_to_remove
                                            )
    res = np.array(proportions_in_times)
    _, nr_of_reagents = res.shape

    for reag_nr in range(nr_of_reagents):
        plt.plot(time_range, res[:,reag_nr], '.', color=color_list[reag_nr],
                label = names[reag_nr])

    if which_reagent_to_remove is None:
        plt.title('Magnetstein with complete library \n', size=20)
    else:
        plt.title('\n'.join(wrap('Magnetstein with incomplete library: missing ' + substances_names[which_reagent_to_remove],
                                40)
                           ),
                  size=20)
        plt.plot(time_range, 1-res.sum(axis=1), '.', color='#5e5e5e', label='removed signal')
        
    plt.xlabel('Time [min]')
    plt.ylabel('Proportion')
    leg = plt.legend(fontsize=15)
    for handle in leg.legend_handles:
        handle.set_markersize(15)

    plt.tight_layout()

    if results_path is None:
        plt.show()
    else:
        if which_reagent_to_remove is None:
            plt.savefig(results_path + 'all_components_kappa_'+str(kappa)+'_'+str(kappa_th)+'.png')
        else:
            plt.savefig(results_path + 'removed_component_'+str(which_reagent_to_remove)+'_kappa_'+str(kappa)+'_'+str(kappa_th)+'.png')
    plt.clf()


def add_in_proportion(main_peak, small_peak, chemical_shift_axis, small_peak_proportion):
    
    component_0 = NMRSpectrum(confs = list(zip(chemical_shift_axis, main_peak)))
    component_0.normalize()
    
    component_1 = NMRSpectrum(confs = list(zip(chemical_shift_axis, small_peak)))
    component_1.normalize(small_peak_proportion / (1 + small_peak_proportion))
    
    result_intensities = np.array(component_0.confs)[:,1] + np.array(component_1.confs)[:,1]
    result = NMRSpectrum(confs = list(zip(chemical_shift_axis, result_intensities)))
    result.normalize()
    
    return result


def second_peak_location(overlap_proportion, first_peak_position, scale):
    return 2*scale*np.tan((math.pi*0.5)*(1-overlap_proportion)) + first_peak_position


def create_second_peak(overlap_proportion, first_peak_position, scale, chemical_shift_axis):

    loc = second_peak_location(overlap_proportion, first_peak_position, scale)
    second_peak = cauchy.pdf(np.array(chemical_shift_axis), 
                                        loc = loc, 
                                        scale = scale)
    res = NMRSpectrum(confs=list(zip(chemical_shift_axis, second_peak)))
    res.normalize()
    return loc, res


def create_mixture(component_0, component_1, ground_truth_proportions, shift=0):
    
    assert np.all(np.array(component_0.confs)[:,0] == np.array(component_1.confs)[:,0])
    chemical_shift_axis = np.array(component_0.confs)[:,0] + shift
    
    component_0_intensities = np.array(component_0.confs)[:,1]
    component_1_intensities = np.array(component_1.confs)[:,1]
    mixture_intensities = ground_truth_proportions[0]*component_0_intensities + ground_truth_proportions[1]*component_1_intensities
    mixture = NMRSpectrum(confs=list(zip(chemical_shift_axis, mixture_intensities)))
    mixture.normalize()
    
    return mixture

def run_overlapping_simulations(overlap_proportion, shift, chemical_shift_axis, 
                                first_peak_position, small_peak_position,
                                small_peak_proportion, scale, shift_diff, 
                                ground_truth_proportions, results_path):

    abs_error_df = pd.DataFrame(index=overlap_proportion, columns=shift)

    for op in overlap_proportion:
        for sh in shift:
            print('------------------------------')
            print('Overlap proportion: ' + str(op))
            print('Shift: ' + str(sh))

            #Creating components
            first_peak = cauchy.pdf(np.array(chemical_shift_axis), 
                                            loc = first_peak_position, 
                                            scale = scale)
            small_peak = cauchy.pdf(np.array(chemical_shift_axis), 
                                            loc = small_peak_position, 
                                            scale = scale)

            component_0 = add_in_proportion(first_peak, small_peak, chemical_shift_axis, small_peak_proportion)

            second_peak_position, component_1 = create_second_peak(op, first_peak_position, scale, chemical_shift_axis)

            #Creating mixture
            mixture = create_mixture(component_0, component_1, ground_truth_proportions, sh)

            #Cutting out library from the mixture
            where_to_cut = (second_peak_position + first_peak_position)/2
            component_0_confs = [(pos, intensity) for pos, intensity in mixture.confs if pos<where_to_cut]
            component_1_confs = [(pos, intensity) for pos, intensity in mixture.confs if pos>where_to_cut]
            component_0 = NMRSpectrum(confs = component_0_confs)
            component_0.normalize()
            component_1 = NMRSpectrum(confs = component_1_confs)
            component_1.normalize()

            #Estimation
            estimation = estimate_proportions(mixture, [component_0, component_1], what_to_compare='area', 
                                          solver=pulp.GUROBI(msg=False),
                                         #MTD=sh+shift_diff, MTD_th=sh+shift_diff)
                                            MTD=0.25, MTD_th=0.22)
                                              
            abs_err = sum([abs(el-gt) for el, gt in zip(estimation['proportions'], ground_truth_proportions)])
            print(estimation['proportions'])
            
            abs_error_df.loc[op, sh] = abs_err

    abs_error_df = abs_error_df.apply(pd.to_numeric)

    if results_path is not None:
        name = results_path + 'estimation_error_small_peak_proportion_'+ str(small_peak_proportion) + \
                '_shift_diff_' + str(shift_diff) + '.csv'
        abs_error_df.to_csv(name)

    return abs_error_df


def visualise_overlap_and_shift(abs_error_df, results_path, small_peak_proportion, shift_diff):

    fig, ax = plt.subplots(figsize=(6, 5))

    # Display heatmap
    cax = ax.imshow(abs_error_df.values, cmap='RdYlGn_r', aspect='auto', vmax=1)

    # Add colorbar to the figure
    cbar = fig.colorbar(cax, ax=ax, label = 'Absolute error')
    cbar.set_label(label='Absolute error', fontsize=15)

    # Set axis ticks and labels
    ax.set_xticks(range(len(abs_error_df.columns)))
    ax.set_xticklabels([np.round(el, 2) for el in pd.to_numeric(abs_error_df.columns)], rotation=90)
    ax.set_yticks(range(len(abs_error_df.index))[::10])
    ax.set_yticklabels([np.round(el,3) for el in abs_error_df.index[::10]])
    ax.set_xlabel("Shift", size=15)
    ax.set_ylabel("Proportion of overlap", size=15)

    # Layout adjustment
    fig.tight_layout()

    plt.savefig(results_path + 'error_heatmap_small_peak_proportion_'+ str(small_peak_proportion) + \
                '_shift_' + str(shift_diff) + '.png')
