import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from masserstein import NMRSpectrum
from utils import load_mixture_time_data_v3, load_spectrum, load_and_normalize_spectrum


def convert(txt):
    k = float(txt)
    if k%1 ==0:
        return int(k)
    return k

def get_all_kappa_values(result_path):

    all_kappas = []
    all_kappas_th = []

    for file in os.listdir(result_path):
        if file.startswith('noise_proportions_in_times'):
            kappa = convert(file.split('_')[4])
            all_kappas.append(kappa)
            kappa_th = convert(file.split('_')[5][:-4])
            all_kappas_th.append(kappa_th)
            
    all_kappas = sorted(list(set(all_kappas)))
    all_kappas_th = sorted(list(set(all_kappas_th)))

    return all_kappas, all_kappas_th

def remove_kappa_values(all_kappas, all_kappas_th, experiment_name):

    if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
        kappas_to_remove = [0.05, 0.5, 0.75, 1, 1.2]
        kappas_th_to_remove = [0.05, 0.5, 0.75, 1, 1.2]

    elif experiment_name == 'PMG 284 monitoring with full mixture spectrum':
        kappas_to_remove = [0.5]
        kappas_th_to_remove = [0.5]

    else:
        kappas_to_remove = [1.2, 0.75]
        kappas_th_to_remove = [1.2, 0.75]

    all_kappas = sorted(list(set(all_kappas) - set(kappas_to_remove)))
    all_kappas_th = sorted(list(set(all_kappas_th) - set(kappas_th_to_remove)))

    return all_kappas, all_kappas_th

def create_results_dict(results_path, all_kappas, all_kappas_th, filename='proportions_in_times_'):

    results_dict = {}

    for kappa in all_kappas:
        for kappa_th in all_kappas_th:
            try:
                with open(results_path + filename +
                                      str(kappa)+'_'+str(kappa_th)+'.pkl', 'rb') as f:
                    y = pickle.load(f)
                    y = np.array(y)
                    results_dict[(kappa, kappa_th)] = y
            except FileNotFoundError:
                results_dict[(kappa, kappa_th)] = None

    return results_dict

def load_mnova_integrals(mnova_integrals_path, integrals_separators, experiment_name):

    integrals = pd.read_csv(mnova_integrals_path, sep=integrals_separators)

    if experiment_name == 'PMG 287 monitoring with full mixture spectrum':
        how_to_divide = 'whole spectrum'
    else:
        how_to_divide = 'sum'

    if experiment_name == 'Saccharose hydrolysis' or experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
        integrals = integrals.iloc[1:,:]
        integrals.columns = ['saccharose', 'alpha-glucose', 'beta-glucose', 'fructose']
        integrals = integrals.apply(pd.to_numeric)
        nominator = np.array(integrals[['saccharose', 'alpha-glucose', 'beta-glucose', 'fructose']])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['whole spectrum']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals[['saccharose', 'alpha-glucose', 
                                              'beta-glucose', 'fructose']]).sum(axis=1).reshape(-1,1)
        mnova_integrals_proportions = nominator/denominator
        
    if experiment_name == 'PMG 284 monitoring' or experiment_name == 'PMG 284 monitoring with full mixture spectrum':
        integrals = integrals.iloc[1:,[2,6,8,10]]
        integrals.rename(columns={'Y(X)': 'penten',
                                 'Y2(X)': 'silan',
                                 'Y3(X)': 'produkt',
                                 'Unnamed: 10': 'cale widmo'}, inplace=True)
        integrals = integrals[['penten', 'silan', 'produkt', 'cale widmo']]
        integrals = integrals.apply(pd.to_numeric)
        nominator = np.array(integrals[['penten', 'silan', 'produkt']])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['cale widmo']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals[['penten', 'silan', 'produkt']]).sum(axis=1).reshape(-1,1)
        mnova_integrals_proportions = nominator/denominator
        
    if experiment_name == 'PMG 287 monitoring' or experiment_name == 'PMG 287 monitoring with full mixture spectrum':
        integrals = integrals[['hexene sum', 'silane sum', 'product sum', 'all']].iloc[1:]
        for col in integrals.columns:
            integrals[col] = integrals[col].apply(lambda x: float(x.replace(',', '.')))
            
        nominator = np.array(integrals[['hexene sum', 'silane sum', 'product sum']])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['all']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals[['hexene sum', 'silane sum', 'product sum']]).sum(axis=1).reshape(-1,1)
        mnova_integrals_proportions = nominator/denominator

    return mnova_integrals_proportions


def load_python_integrals(python_integrals_path, python_integrals_separators, experiment_name):

    integrals = pd.read_csv(python_integrals_path, sep=python_integrals_separators)

    if experiment_name == 'PMG 287 monitoring with full mixture spectrum':
        how_to_divide = 'whole spectrum'
    else:
        how_to_divide = 'sum'

    if experiment_name == 'Saccharose hydrolysis' or experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
        integrals = integrals.apply(pd.to_numeric)
        nominator = np.array(integrals.iloc[:,:4])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['whole_spectrum']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals.iloc[:,:4].sum(axis=1)).reshape(-1,1)
        python_integrals_proportions = nominator/denominator
        
    if experiment_name == 'PMG 284 monitoring' or experiment_name == 'PMG 284 monitoring with full mixture spectrum':
        integrals = integrals.apply(pd.to_numeric)
        nominator = np.array(integrals.iloc[:,:3])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['whole_spectrum']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals.iloc[:,:3]).sum(axis=1).reshape(-1,1)
        python_integrals_proportions = nominator/denominator
        
    if experiment_name == 'PMG 287 monitoring' or experiment_name == 'PMG 287 monitoring with full mixture spectrum':
        integrals = integrals.apply(pd.to_numeric)
        nominator = np.array(integrals.iloc[:,:3])
        if how_to_divide == 'whole spectrum':
            denominator = np.array(integrals[['whole_spectrum']])
        elif how_to_divide == 'sum':
            denominator = np.array(integrals.iloc[:,:3]).sum(axis=1).reshape(-1,1)
        python_integrals_proportions = nominator/denominator

    return python_integrals_proportions

def plot_proportion_against_time_single_component(
                                                    experiment_name,
                                                    component_nr, 
                                                    all_kappas, 
                                                    all_kappas_th, 
                                                    results_dict,
                                                    substances_names, 
                                                    time_range,
                                                    saturated_colors_for_components,
                                                    path_to_save=None
                                                ):

    fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')
    fig.suptitle('Proportion of ' + substances_names[component_nr], fontsize=50)
    fig.set_size_inches(30, 18, forward=True)

    fig.text(0.5, 0.92, 'Kappa components', ha='center', size=25)
    fig.text(0.06, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=25)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]
            try:
                y = results_dict[(kappa, kappa_th)][:, component_nr]
                axx.plot(time_range, y, 'p', markersize=5, 
                         color=saturated_colors_for_components[component_nr])
            except TypeError:
                axx.plot(0,0)
            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]       

    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=20)
        ax.set_ylim(0, 1)
        
    for ax, col in zip(axs[-1], cols):
            ax.set_xlabel('Time [min]', size=10)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=20)

    if path_to_save is not None:
        fig.savefig(path_to_save + 'comparison_for_different_kappas_'+ substances_names[component_nr]+'.png')
    fig.show()

def plot_proportion_against_time_all_components(    
                                                experiment_name,
                                                all_kappas, 
                                                all_kappas_th, 
                                                results_dict,
                                                substances_names,
                                                time_range,
                                                saturated_colors_for_components,
                                                path_to_save=None
                                                ):

    if experiment_name == 'Saccharose hydrolysis with full mixture spectrum' or experiment_name == 'PMG 284 monitoring with full mixture spectrum':
        fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all')
    else:
        fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')

    #fig.suptitle('Proportion of substrates and product', fontsize=50)

    ylims_dict_sacch_hydr_full_mix = {0.005: 0.05, 0.01: 0.05, 0.1: 0.21, 0.25: 0.63, 0.5:1.}

    fig.text(0.5, 0.93, 'Kappa components', ha='center', size=23)
    fig.text(0.04, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=23)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]
            try:
                y = results_dict[(kappa, kappa_th)]
                for i in range(y.shape[1]):
                    if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
                        axx.axis(ymin=0, ymax=ylims_dict_sacch_hydr_full_mix[kappa])
                        axx.plot(time_range, y[:,i], 'p', markersize=1, 
                                 label=substances_names[i],
                                color=saturated_colors_for_components[i])
                    else:
                        axx.plot(time_range, y[:,i], 'p', markersize=1, 
                                 label=substances_names[i],
                                color=saturated_colors_for_components[i])
            except TypeError:
                axx.plot(0,0)
            except AttributeError:
                axx.plot(0,0)
            axx.tick_params(axis="y", labelsize=9)
            axx.tick_params(axis="x", labelsize=9)


            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]     


    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=15)
        if experiment_name == 'Saccharose hydrolysis with full mixture spectrum' or experiment_name == 'PMG 287 monitoring with full mixture spectrum' or experiment_name == 'PMG 284 monitoring with full mixture spectrum':
            pass
        else:
            ax.set_ylim(0,1)

    for ax, col in zip(axs[-1], cols):
        ax.set_xlabel('Time [min]', size=15)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=15)
    handles, labels = ax.get_legend_handles_labels()

    fig.legend(handles, labels, loc='upper right', markerscale=5, prop={'size': 8})

    plt.gcf().set_size_inches(11.7, 7.02)

    if path_to_save is not None:
        fig.savefig(path_to_save + 'comparison_for_different_kappas_components_separately.png', dpi=200)
    fig.show()


def plot_proportion_against_time_all_components_plus_integrals(
                                                                experiment_name,
                                                                all_kappas, 
                                                                all_kappas_th, 
                                                                results_dict,
                                                                substances_names,
                                                                time_range,
                                                                saturated_colors_for_components,
                                                                colors_for_components,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                python_or_mnova = 'Mnova',
                                                                path_to_save=None
                                                                ):

    ylims_dict_sacch_hydr_full_mix = {0.01: 0.05, 0.1: 0.21, 0.25: 0.63, 0.5:1.}

    if python_or_mnova == 'Mnova':
        integrals_proportions = mnova_integrals_proportions
        time_range_integrals = time_range_integrals_mnova
    elif python_or_mnova == 'Python':
        integrals_proportions = python_integrals_proportions
        time_range_integrals = time_range_integrals_python

    if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
        fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all')
    else:
        fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')
    fig.suptitle('Proportions of substrates and products', fontsize=80)
    fig.set_size_inches(40, 25, forward=True)

    fig.text(0.5, 0.92, 'Kappa components', ha='center', size=40)
    fig.text(0.06, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=40)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]
            y = results_dict[(kappa, kappa_th)]
            if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
                pass
            else:
                for i in range(integrals_proportions.shape[1]):
                    axx.plot(time_range_integrals, integrals_proportions[:,i], 'p', 
                             markersize=5, label=python_or_mnova + ' normalized integral for ' + substances_names[i],
                             color=colors_for_components[i]
                             #alpha=0.1
                            )
            for i in range(y.shape[1]):
                if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
                    axx.axis(ymin=0, ymax=ylims_dict_sacch_hydr_full_mix[kappa])
                    axx.plot(time_range, y[:,i], 'p', markersize=5, 
                            label=substances_names[i],
                            color=saturated_colors_for_components[i])
                else:
                    axx.plot(time_range, y[:,i], 'p', markersize=5, 
                             label=substances_names[i],
                            color=saturated_colors_for_components[i])
                axx.tick_params(axis='both', which='major', labelsize=20)


    #         except AttributeError:
    #             axx.plot(0,0)
            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]       

    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=40)
        if experiment_name == 'Saccharose hydrolysis with full mixture spectrum' or experiment_name == 'PMG 287 monitoring with full mixture spectrum':
            pass
        else:
            ax.set_ylim(0,1)
        
    for ax, col in zip(axs[-1], cols):
        ax.set_xlabel('Time [min]', size=20)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=40)
    handles, labels = ax.get_legend_handles_labels()

    legend = fig.legend(handles, labels, loc='upper right', prop={'size': 30}, markerscale=5)

    # for lh in legend.legendHandles[i+1:]: 
    #     lh.set_alpha(0.5)

    if path_to_save is not None:
        if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
            fig.savefig(path_to_save + 'comparison_for_different_kappas_components_separately.png')

        else:
            fig.savefig(path_to_save + 'comparison_for_different_kappas_components_separately_plus_integrals_' + python_or_mnova + '.png')
    
    fig.show()

def plot_proportion_and_integrals_against_time_chosen_kappas(
                                                            experiment_name,
                                                            results_path,
                                                            chosen_kappa,
                                                            chosen_kappa_th,
                                                            substances_names,
                                                            time_range,
                                                            saturated_colors_for_components,
                                                            colors_for_components,
                                                            mnova_integrals_proportions,
                                                            python_integrals_proportions,
                                                            time_range_integrals_mnova,
                                                            time_range_integrals_python,
                                                            include_ticks=True,
                                                            python_or_mnova='Mnova',
                                                            path_to_save=None
                                                            ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    ylims_dict_sacch_hydr_full_mix = {0.01: 0.05, 0.1: 0.21, 0.25: 0.63, 0.5:1.}

    python_or_mnova = 'Mnova'
    if python_or_mnova == 'Mnova':
        integrals_proportions = mnova_integrals_proportions
        time_range_integrals = time_range_integrals_mnova
    elif python_or_mnova == 'Python':
        integrals_proportions = python_integrals_proportions
        time_range_integrals = time_range_integrals_python

    include_integrals = True
    if include_integrals:
        for i in range(integrals_proportions.shape[1]):
            plt.plot(time_range_integrals, integrals_proportions[:,i], 'p', 
                     markersize=5, label=python_or_mnova + ' integration, ' + substances_names[i],
                     color = colors_for_components[i]
                     #alpha=0.1
                    )
    for i in range(y.shape[1]):
        if experiment_name == 'Saccharose hydrolysis with full mixture spectrum':
            plt.plot(time_range, y[:, i], 'p', 
                     label = 'Magnetstein, ' + substances_names[i],
                    color=saturated_colors_for_components[i])
            plt.ylim(0, ylims_dict_sacch_hydr_full_mix[chosen_kappa])
        else:
            plt.plot(time_range, y[:, i], 'p', 
                     label = 'Magnetstein, ' + substances_names[i],
                    color=saturated_colors_for_components[i])
            if experiment_name == 'PMG 287 monitoring with full mixture spectrum':
                plt.ylim(0, 0.4)
            else:
                plt.ylim(0,1)

    if include_ticks:
        plt.xlabel('Time [min]', size=15)
        plt.ylabel('Proportion', size=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
    else:
        plt.xlabel('Time', size=15)
        plt.xticks([])
        plt.yticks([])
        
    leg = plt.legend(title='Method, substance')    
    plt.gcf().set_size_inches(11.7, 5.85)

    # for lh in leg.legendHandles[:y.shape[1]]: 
    #     lh.set_alpha(0.4)

    if path_to_save is not None:
        if include_ticks:
            plt.savefig(path_to_save + 'all_components_on_one_plot_plus_integrals_'+
                                                                python_or_mnova +
                                                                '_kappa_' +
                                                            str(chosen_kappa)+'_kappa_th_' + 
                                                            str(chosen_kappa_th) +'.png', dpi=200)
        else:
            plt.savefig(path_to_save + 'all_components_on_one_plot_plus_integrals_'+
                                                                python_or_mnova +
                                                                '_kappa_' +
                                                            str(chosen_kappa)+'_kappa_th_' + 
                                                            str(chosen_kappa_th) +'_no_labels.png', dpi=200)
    plt.show()

def plot_proportion_against_time_chosen_kappas(
                                                experiment_name,
                                                results_path,
                                                chosen_kappa,
                                                chosen_kappa_th,
                                                substances_names,
                                                time_range,
                                                saturated_colors_for_components,
                                                include_ticks=True,
                                                path_to_save=None
                                                ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)


    for i in range(y.shape[1]):
        plt.plot(time_range, y[:, i], 'p', 
                 label = 'Magnetstein, ' + substances_names[i],
                color=saturated_colors_for_components[i])

    plt.xlabel('Time', size=15)
    plt.ylabel('Relative amount', size=15)
    # plt.xticks([])
    # plt.yticks([])
     
    plt.gcf().set_size_inches(10, 5)

    # for lh in leg.legendHandles[:y.shape[1]]: 
    #     lh.set_alpha(0.4)

    if path_to_save is not None:
        plt.savefig(path_to_save + 'all_components_on_one_plot_magnetstein_only'+
                                                            '_kappa_' +
                                                        str(chosen_kappa)+'_kappa_th_' + 
                                                        str(chosen_kappa_th) +'.png', dpi=300)
    plt.show()

def plot_proportion_and_all_integrals_against_time_chosen_kappas(
                                                                experiment_name,
                                                                results_path,
                                                                chosen_kappa,
                                                                chosen_kappa_th,
                                                                substances_names,
                                                                time_range,
                                                                saturated_colors_for_components,
                                                                colors_for_components,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                path_to_save=None
                                                                ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)


    for i in range(mnova_integrals_proportions.shape[1]):
        plt.plot(time_range_integrals_mnova, mnova_integrals_proportions[:,i], 'p', 
                 markersize=5, label='MNova integration, ' + substances_names[i],
                 color = colors_for_components[i]
                 #alpha=0.1
                )
    for i in range(python_integrals_proportions.shape[1]):
        plt.plot(time_range_integrals_python, python_integrals_proportions[:,i], 'v', 
                 markersize=5, label='Python integration, ' + substances_names[i],
                 color = colors_for_components[i]
                 #alpha=0.1
                )
    for i in range(y.shape[1]):
        plt.plot(time_range, y[:, i], 'p', 
                 label = 'Magnetstein, ' + substances_names[i],
                color=saturated_colors_for_components[i])

        plt.title('Proportions in time')
        
    plt.xlabel('Time [min]')
        
    if experiment_name == 'PMG 287 monitoring with full mixture spectrum':
        plt.ylim(0, 0.4)
    else:
        plt.ylim(0,1)
        
    leg = plt.legend(title='Method, substance')    
    # for lh in leg.legendHandles[:y.shape[1]]: 
    #     lh.set_alpha(0.4)
    plt.gcf().set_size_inches(10, 5)

    if path_to_save is not None:
        plt.savefig(path_to_save + 'all_components_on_one_plot_plus_integrals_MNova_Python'+
                                                            '_kappa_' +
                                                        str(chosen_kappa)+'_kappa_th_' + 
                                                        str(chosen_kappa_th) +'.png')
    plt.show()


def plot_integrals_against_time_chosen_kappas(
                                                experiment_name,
                                                substances_names,
                                                colors_for_components,
                                                mnova_integrals_proportions,
                                                python_integrals_proportions,
                                                time_range_integrals_mnova,
                                                time_range_integrals_python,
                                                python_or_mnova='Mnova',
                                                path_to_save=None
                                                ):


    if python_or_mnova == 'Mnova':
        integrals_proportions = mnova_integrals_proportions
        time_range_integrals = time_range_integrals_mnova
    elif python_or_mnova == 'Python':
        integrals_proportions = python_integrals_proportions
        time_range_integrals = time_range_integrals_python

        
    for i in range(integrals_proportions.shape[1]):
        plt.plot(time_range_integrals, integrals_proportions[:,i], 'p', 
                 markersize=5, label=python_or_mnova + ' normalized integral for ' + substances_names[i],
                 color = colors_for_components[i]#, alpha=0.03
                )

    plt.title('Proportions in time')

    plt.xlabel('Time [min]')
        
    #plt.ylim(0,1)
        
    leg = plt.legend()    

    plt.gcf().set_size_inches(10, 5)

    if path_to_save is not None:
        plt.savefig(path_to_save + 'all_components_integrals_on_one_plot_' +
                                                            python_or_mnova  + '.png')
    plt.show()


def plot_proportion_against_time_all_components_added_plus_integrals(
                                                                    experiment_name,
                                                                    substrates_numbers,
                                                                    product_numbers,
                                                                    all_kappas, 
                                                                    all_kappas_th, 
                                                                    results_dict,
                                                                    substances_names,
                                                                    time_range,
                                                                    saturated_colors_for_components,
                                                                    colors_for_components,
                                                                    mnova_integrals_proportions,
                                                                    python_integrals_proportions,
                                                                    time_range_integrals_mnova,
                                                                    time_range_integrals_python,
                                                                    python_or_mnova = 'Mnova',
                                                                    path_to_save=None
                                                                    ):


    if python_or_mnova == 'Mnova':
        integrals_proportions = mnova_integrals_proportions
        time_range_integrals = time_range_integrals_mnova
    elif python_or_mnova == 'Python':
        integrals_proportions = python_integrals_proportions
        time_range_integrals = time_range_integrals_python

    colors_list = ['grey', 'purple']

    fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')
    fig.suptitle('Proportion of substrates and product', fontsize=50)
    fig.set_size_inches(40, 25, forward=True)

    fig.text(0.5, 0.92, 'Kappa components', ha='center', size=25)
    fig.text(0.06, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=25)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]

            #substrates
            y = 0
            for nr in substrates_numbers:
                y = y + results_dict[(kappa, kappa_th)][:, nr]
            axx.plot(time_range, y, 'p', markersize=5,
                    label='substrates',
                        color = colors_list[0])
            #product
            y = 0
            for nr in product_numbers:
                y = y + results_dict[(kappa, kappa_th)][:, nr]
            axx.plot(time_range, y, 'p', markersize=5,
                    label='product',
                        color = colors_list[1])

            #integrals
            y = results_dict[(kappa, kappa_th)]
            axx.plot(time_range_integrals, integrals_proportions[:,substrates_numbers].sum(axis=1), 'p', 
                         markersize=5, label=python_or_mnova + ' normalized integral for substrates',
                         color = colors_list[0], alpha=0.03
                        )
            axx.plot(time_range_integrals, integrals_proportions[:,product_numbers].sum(axis=1), 'p', 
                         markersize=5, label=python_or_mnova + ' normalized integral for products',
                         color = colors_list[1], alpha=0.03
                        )
            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]       

    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=20)
        ax.set_ylim(0,1)
        
    for ax, col in zip(axs[-1], cols):
        ax.set_xlabel('Time [min]', size=10)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=20)
    handles, labels = ax.get_legend_handles_labels()

    legend = fig.legend(handles, labels, loc='upper right', prop={'size': 25}, markerscale=5)

    if path_to_save is not None:
        fig.savefig(path_to_save+
                        'comparison_for_different_kappas_components_together_plus_integrals_' +
                        python_or_mnova + 
                    '.png')
    fig.show()


def plot_proportion_against_time_chosen_components_added(
                                                        experiment_name,
                                                        components_to_add,
                                                        all_kappas, 
                                                        all_kappas_th, 
                                                        results_dict,
                                                        substances_names,
                                                        time_range,
                                                        path_to_save=None
                                                        ):

    components_numbers = components_to_add
    names = [substances_names[nr] for nr in components_numbers]

    fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')
    fig.suptitle('Added proportion of components: ' + str(names), fontsize=50)
    fig.set_size_inches(30, 18, forward=True)

    fig.text(0.5, 0.92, 'Kappa components', ha='center', size=25)
    fig.text(0.06, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=25)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]
            try:
                y = 0
                for nr in components_numbers:
                    y = y + results_dict[(kappa, kappa_th)][:, nr]
                axx.plot(time_range, y, 'p', markersize=5,
                        color='grey')

            except TypeError:
                axx.plot(0,0)
            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]       

    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=20)
        ax.set_ylim(0,1)
        
    for ax, col in zip(axs[-1], cols):
        ax.set_xlabel('Time [min]', size=10)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=20)

    if path_to_save is not None:
        fig.savefig(path_to_save + 'comparison_for_different_kappas_sum_of_components_' +
                                                      str(names) +'.png')
    fig.show()

def plot_proportion_against_time_single_component_chosen_kappa(
                                                                experiment_name,
                                                                component_nr,
                                                                results_path,
                                                                chosen_kappa,
                                                                chosen_kappa_th,
                                                                substances_names,
                                                                time_range,
                                                                saturated_colors_for_components,
                                                                path_to_save=None
                                                                ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    plt.plot(time_range, y[:, component_nr], 'p',
            color = saturated_colors_for_components[component_nr])

    #plt.title('Proportion of '+substances_names[experiment_name][component_nr])
    plt.title('Proportion of ' + substances_names[component_nr])

    plt.xlabel('Time [min]')
    if path_to_save is not None:
        fig.savefig(path_to_save + 'single_component_chosen_kappas_' +
                                                      substances_names[component_nr] +'.png')
    plt.show()


def plot_proportion_against_time_components_added_together_chosen_kappa(
                                                                        experiment_name,
                                                                        component_numbers,
                                                                        results_path,
                                                                        chosen_kappa,
                                                                        chosen_kappa_th,
                                                                        substances_names,
                                                                        time_range,
                                                                        path_to_save=None
                                                                        ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    plt.plot(time_range, y[:,component_numbers].sum(1), 'p',
            color='grey')

    plt.title('Proportion of components: ' + str(component_numbers))

    plt.xlabel('Time [min]')

    if path_to_save is not None:
        fig.savefig(path_to_save + 'components_added_together_chosen_kappas_kappa_' + str(chosen_kappa) + '_kappa_th_' + 
                                                      str(chosen_kappa_th) + '_components_' +
                                                      str(component_numbers) +'.png')
    plt.show()

def plot_noise_proportion_chosen_kappas(
                                        experiment_name,
                                        results_path,
                                        chosen_kappa,
                                        chosen_kappa_th,
                                        time_range,
                                        path_to_save=None
                                        ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    plt.plot(time_range, 1-y.sum(axis=1), 'p',
            color='grey')
    plt.title('Proportion of noise in mixture')

    plt.xlabel('Time [min]')

    if path_to_save is not None:
        fig.savefig(path_to_save + 'noise_proportion_chosen_kappas_kappa' + str(chosen_kappa) + '_kappa_th_' +
                                                      str(chosen_kappa_th) +'.png')
    plt.show()

def plot_noise_proportion_chosen_kappas_v2(
                                        experiment_name,
                                        results_path,
                                        chosen_kappa,
                                        chosen_kappa_th,
                                        time_range,
                                        path_to_save=None
                                        ):

    with open(results_path + 'noise_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = [sum(l) for l in y]

    plt.plot(time_range, y, 'p',
            color='grey')
    plt.title('Proportion of noise in mixture')
    plt.xlabel('Time [min]')

    if path_to_save is not None:
        fig.savefig(path_to_save + 'noise_proportion_v2_chosen_kappas_kappa' + str(chosen_kappa) + '_kappa_th_' + 
                                                      str(chosen_kappa_th) +'.png')
    plt.show()


def plot_added_components_plus_integrals_against_time_chosen_kappas(
                                                                experiment_name,
                                                                component_numbers,
                                                                product_numbers,
                                                                results_path,
                                                                chosen_kappa,
                                                                chosen_kappa_th,
                                                                time_range,
                                                                mnova_integrals_proportions,
                                                                python_integrals_proportions,
                                                                time_range_integrals_mnova,
                                                                time_range_integrals_python,
                                                                python_or_mnova='Mnova',
                                                                colors_list = ['grey', 'purple'],
                                                                include_integrals = True,
                                                                path_to_save=None
                                                                ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    #Version when we don't have spectrum of product (i.e. we assume that noise in mixture is product)

    if python_or_mnova == 'Mnova':
        integrals_proportions = mnova_integrals_proportions
        time_range_integrals = time_range_integrals_mnova
    elif python_or_mnova == 'Python':
        integrals_proportions = python_integrals_proportions
        time_range_integrals = time_range_integrals_python

    if include_integrals:
        plt.plot(time_range_integrals, integrals_proportions[:,component_numbers].sum(axis=1), 'p', 
                     markersize=5, label=python_or_mnova + ' normalized integral for substartes',
                     color = colors_list[0], alpha=0.03
                    )
        plt.plot(time_range_integrals, integrals_proportions[:,product_numbers].sum(axis=1), 'p', 
                     markersize=5, label=python_or_mnova + ' normalized integral for products',
                     color = colors_list[1], alpha=0.03
                    )

    plt.plot(time_range, y[:,component_numbers].sum(1), 'p', label='Substrates', color = colors_list[0])
    plt.plot(time_range, y[:,product_numbers].sum(1), 'p', label='Product', color = colors_list[1])

    plt.ylim(0,1)
    #plt.plot(list(range(y.shape[0]*10))[::10], 1-y.sum(axis=1), 'p', label='Product')

    leg = plt.legend(loc='upper left', prop = {"size": 7})    

    plt.gcf().set_size_inches(10, 5)

    plt.title('Proportion of substrates/product in time')

    plt.xlabel('Time [min]')
            
    if path_to_save is not None:
        plt.savefig(path_to_save +
                                    'substrates_and_product_together_plus_integrals_' + 
                                    python_or_mnova + 
                                    '_kappa_' + 
                                    str(best_kappa) +
                                    '_kappa_th_' +
                                    str(best_kappa_th)
                                    +'.png')
    plt.show()


def plot_noise_proportion_in_components_chosen_kappas(
                                        experiment_name,
                                        results_path,
                                        chosen_kappa,
                                        chosen_kappa_th,
                                        time_range,
                                        path_to_save=None
                                        ):

    with open(results_path + 'noise_proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    plt.plot(time_range, y, 'p', 
            color='grey')
    plt.title('Proportion of noise in components')

    plt.xlabel('Time [min]')

    if path_to_save is not None:
        fig.savefig(path_to_save + 'noise_proportion_in_components_chosen_kappas_kappa' + str(chosen_kappa) + '_kappa_th_' + 
                                                      str(chosen_kappa_th) +'.png')
    plt.show()


def plot_noise_proportion_against_time (
                                        experiment_name,
                                        all_kappas, 
                                        all_kappas_th, 
                                        substances_names,
                                        results_dict,
                                        time_range,
                                        path_to_save=None
                                        ):

    components_numbers = list(range(len(substances_names)))

    fig, axs = plt.subplots(len(all_kappas), len(all_kappas_th), sharex='all', sharey='all')
    fig.suptitle('Proportion of signal removed from mixture', fontsize=50)
    fig.set_size_inches(30, 18, forward=True)

    fig.text(0.5, 0.92, 'Kappa components', ha='center', size=25)
    fig.text(0.06, 0.5, 'Kappa mixture', va='center', rotation='vertical', size=25)

    for i, ax in enumerate(axs):
        kappa = all_kappas[i]
        for j, axx in enumerate(ax):
            kappa_th = all_kappas_th[j]
            try:
                y = 1
                for nr in components_numbers:
                    y = y - results_dict[(kappa, kappa_th)][:, nr]
                axx.plot(time_range, y, 'p', markersize=5,
                        color='grey')
            except TypeError:
                axx.plot(0,0)
            
    cols = [str(kappa_th) for kappa_th in all_kappas_th]
    rows = [str(kappa) for kappa in all_kappas]       

    for ax, col in zip(axs[0], cols):
        ax.set_title(col, size=20)
        ax.set_ylim(0,1)
        
    for ax, col in zip(axs[-1], cols):
        ax.set_xlabel('Time [min]', size=10)

    for ax, row in zip(axs[:,0], rows):
        ax.set_ylabel(row, rotation=90, size=20)

    if path_to_save is not None:
        fig.savefig(path_to_save + 'proportion_of_removed_signal.png')
    fig.show()


def visualise_signal_removed_from_mixture(
                                            experiment_name,
                                            chosen_kappa,
                                            chosen_kappa_th,
                                            results_path,
                                            mixture_path,
                                            mixture_separator,
                                            moment_of_time = 1,
                                            path_to_save=None
                                            ):

    with open(results_path + 'noise_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    with open(results_path + 'common_horizontal_axis_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        cha = pickle.load(f)
        cha = np.array(cha)

    noise = y[moment_of_time-1,:]
    mixture_time_data = load_mixture_time_data_v3(experiment_name, mixture_path, mixture_separator)
    if experiment_name == 'PMG 284 monitoring' or experiment_name == 'PMG 284 monitoring with full mixture spectrum':
        mixture_time_data.fillna(0, inplace=True)

    mix = load_spectrum(mixture_time_data, moment_of_time)
    bar_width = np.mean([cha[i]-cha[i-1] for i in range(1, len(cha))])

    plt.gca().invert_xaxis()
    plt.bar(cha, noise, alpha=0.7, label='signal removed from mixture', width=bar_width);
    plt.legend()
    NMRSpectrum.plot_all([mix], profile=True);

    if path_to_save is not None:
        plt.savefig(path_to_save + 'noise_removed_from_mixture.png')
    plt.show()

def visualise_signal_removed_from_mixture_plus_reagents(
                                                    experiment_name,
                                                    chosen_kappa,
                                                    chosen_kappa_th,
                                                    results_path,
                                                    mixture_path,
                                                    mixture_separator,
                                                    reagents_path,
                                                    reagents_separator,
                                                    moment_of_time = 1,
                                                    path_to_save=None
                                                    ):

    with open(results_path + 'noise_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)
    with open(results_path + 'common_horizontal_axis_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        cha = pickle.load(f)
        cha = np.array(cha)

    noise = y[moment_of_time-1,:]
    reagents = []
    for i, path in enumerate(reagents_path):
        reag = pd.read_csv(path, sep = reagents_separator)
        if experiment_name == 'PMG 284 monitoring' and i==2:
            reag.fillna(0, inplace=True)
        sp = NMRSpectrum(confs = list(zip(reag.iloc[:,0], reag.iloc[:,1])))
        sp.normalize()
        reagents.append(sp)

    if experiment_name == 'Saccharose hydrolysis':
        products = reagents[-2:]
    else:
        products = reagents[-1:]

    bar_width = np.mean([cha[i]-cha[i-1] for i in range(1, len(cha))])
    plt.gca().invert_xaxis()
    plt.bar(cha, noise, alpha=0.7, label='noise', width=bar_width);
    plt.legend()
    NMRSpectrum.plot_all(products, profile=True)

    if path_to_save is not None:
        plt.savefig(path_to_save + 'noise_removed_from_mixture_plus_reagents.png')
    plt.show()


def visualise_signal_removed_from_components(
                                            experiment_name,
                                            components_nrs,
                                            chosen_kappa,
                                            chosen_kappa_th,
                                            results_path,
                                            reagents_path,
                                            reagents_separator,
                                            saturated_colors_for_components,
                                            moment_of_time = 60,
                                            path_to_save=None
                                            ):

    with open(results_path + 'proportions_in_times_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)
        current_proportions = y[moment_of_time-1, components_nrs]

    with open(results_path + 'noise_in_components_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        y = pickle.load(f)
        y = np.array(y)

    with open(results_path + 'common_horizontal_axis_'+str(chosen_kappa)+'_'+str(chosen_kappa_th)+'.pkl', 'rb') as f:
        cha = pickle.load(f)
        cha = np.array(cha)

    noise_in_components = y[moment_of_time-1,:]

    reagents = []
    for i, path in enumerate(reagents_path):
        reag = pd.read_csv(path, sep = reagents_separator)
        if experiment_name == 'PMG 284 monitoring' and i==2:
            reag.fillna(0, inplace=True)
        sp = NMRSpectrum(confs = list(zip(reag.iloc[:,0], reag.iloc[:,1])))
        sp.normalize()
        reagents.append(sp)

    bar_width = np.mean([cha[i]-cha[i-1] for i in range(1, len(cha))])

    plt.gca().invert_xaxis()
    plt.bar(cha, noise_in_components, alpha=0.7, label='noise in components', width=bar_width, color='grey');
    for i in components_nrs:
        plt.plot(np.array(reagents[i].confs)[:,0], np.array(reagents[i].confs)[:,1],
                        color=saturated_colors_for_components[i]);
    plt.legend()
    if path_to_save is not None:
        plt.savefig(path_to_save + 'noise_removed_from_components.png')
    plt.show()

def plot_mixture_and_reagents(
                                experiment_name,
                                results_path,
                                mixture_path,
                                mixture_separator,
                                reagents_path,
                                reagents_separator,
                                moment_of_time = 1,
                                path_to_save=None
                                ):
    
    reagents = []
    for i, path in enumerate(reagents_path):
        reag = pd.read_csv(path, sep = reagents_separator)
        sp = NMRSpectrum(confs = list(zip(reag.iloc[:,0], reag.iloc[:,1])))
        sp.normalize()
        reagents.append(sp)

    mixture_time_data = load_mixture_time_data_v3(experiment_name, mixture_path, mixture_separator)
    if experiment_name == 'PMG 284 monitoring':
        mixture_time_data.fillna(0, inplace=True)

    plt.gca().invert_xaxis()
    NMRSpectrum.plot_all(reagents + [load_and_normalize_spectrum(mixture_time_data, moment_of_time)], profile=True)
    if path_to_save is not None:
        plt.savefig(path_to_save + 'mixture_and_reagents.png')
    plt.show()

