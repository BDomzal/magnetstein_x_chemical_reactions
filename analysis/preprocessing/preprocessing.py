# SETTINGS AND IMPORTS

import json

with open('../../config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]

import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *


# SACCHAROSE HYDROLYSIS

experiment_name = 'Saccharose hydrolysis'

MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]


mixture_time_data = load_mixture_time_data_v2(MIXTURE_PATH, MIXTURE_SEPARATOR, last_columns_to_skip=1)


# saccharose: (5.39173, 5.44305)
# alpha-glucose: (5.2178, 5.26134)
# beta-glucose: (4.62026, 4.67207)
# fructose: (3.97917, 4.01542)


reagent0 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[0]])
reagent1 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[1]])
reagent2 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[2]])
reagent3 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[3]])

mixture_time_data = cut_mixture_time_data_to_regions_v2(mixture_time_data, INTEGRATION_INTERVALS)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_saccharose.csv', reagent0, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_alpha_glucose.csv', reagent1, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_beta_glucose.csv', reagent2, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_fructose.csv', reagent3, delimiter = '\t')

np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture.csv', mixture_time_data, delimiter = '\t')


# SACCHAROSE HYDROLYSIS WITH FULL MIXTURE SPECTRUM 


experiment_name = 'Saccharose hydrolysis with full mixture spectrum'


MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]


mixture_time_data = load_mixture_time_data(MIXTURE_PATH, MIXTURE_SEPARATOR, last_columns_to_skip=1)


# saccharose: (5.39173, 5.44305)
# alpha-glucose: (5.2178, 5.26134)
# beta-glucose: (4.62026, 4.67207)
# fructose: (3.97917, 4.01542)


reagent0 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[0]])
reagent1 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[1]])
reagent2 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[2]])
reagent3 = create_reagent_spectrum_v2(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[3]])


#baseline correction
# for colname in ['t' + str(nb) for nb in range(1, mixture_time_data.shape[1])]:
#     mixture_time_data[colname] = mixture_time_data[colname].apply(lambda x: x-0.0395)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_saccharose.csv', reagent0, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_alpha_glucose.csv', reagent1, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_beta_glucose.csv', reagent2, delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_fructose.csv', reagent3, delimiter = '\t')

np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture.csv', mixture_time_data, delimiter = '\t')


# PMG 284 MONITORING


experiment_name = 'PMG 284 monitoring'


MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]


if experiment_name == 'PMG 284 monitoring':
    whole_mixture_or_regions = 'regions'
elif experiment_name == 'PMG 284 monitoring with full mixture spectrum':
    whole_mixture_or_regions = 'whole'


mixture_time_data = load_mixture_time_data(MIXTURE_PATH, MIXTURE_SEPARATOR)


reagent0 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[0]])
reagent1 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[1]])
reagent2 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[2]])
reagents_sp = [reagent0, reagent1, reagent2]


if whole_mixture_or_regions == 'whole':
    pass
elif whole_mixture_or_regions == 'regions':
    mixture_time_data = cut_mixture_time_data_to_regions(mixture_time_data, INTEGRATION_INTERVALS)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_pentene.csv', np.array(reagent0.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_triethylsilane.csv', np.array(reagent1.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_product.csv', np.array(reagent2.confs), delimiter = '\t')

if whole_mixture_or_regions == 'whole':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_whole.csv', mixture_time_data, delimiter = '\t')
    
elif whole_mixture_or_regions == 'regions':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture.csv', mixture_time_data, delimiter = '\t')


# PMG 284 MONITORING WITH FULL MIXTURE SPECTRUM


experiment_name = 'PMG 284 monitoring with full mixture spectrum'


MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]


if experiment_name == 'PMG 284 monitoring':
    whole_mixture_or_regions = 'regions'
elif experiment_name == 'PMG 284 monitoring with full mixture spectrum':
    whole_mixture_or_regions = 'whole'


mixture_time_data = load_mixture_time_data(MIXTURE_PATH, MIXTURE_SEPARATOR)


reagent0 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[0]])
reagent1 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[1]])
reagent2 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, [INTEGRATION_INTERVALS[2]])
reagents_sp = [reagent0, reagent1, reagent2]



if whole_mixture_or_regions == 'whole':
    pass
elif whole_mixture_or_regions == 'regions':
    mixture_time_data = cut_mixture_time_data_to_regions(mixture_time_data, INTEGRATION_INTERVALS)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_pentene.csv', np.array(reagent0.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_triethylsilane.csv', np.array(reagent1.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_product.csv', np.array(reagent2.confs), delimiter = '\t')

if whole_mixture_or_regions == 'whole':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_whole.csv', mixture_time_data, delimiter = '\t')
    
elif whole_mixture_or_regions == 'regions':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture.csv', mixture_time_data, delimiter = '\t')


# PMG 287 MONITORING


experiment_name = 'PMG 287 monitoring'


MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]
FULL_PPM_INTERVAL_PMG_287 = config["full_ppm_interval_PMG_287"]


if experiment_name == 'PMG 287 monitoring':
    whole_mixture_or_regions = 'regions'
elif experiment_name == 'PMG 287 monitoring with full mixture spectrum':
    whole_mixture_or_regions = 'whole'


mixture_time_data = load_mixture_time_data(MIXTURE_PATH, MIXTURE_SEPARATOR, last_columns_to_skip=1)


# Integral(3.343730, 3.256018) product
# Integral(3.722376, 3.681636) product
# Integral(3.756886, 3.722376) product
# Integral(4.078496, 3.994139) product

# Integral(3.406997, 3.343730) silane
# Integral(3.832615, 3.756886) silane
# Integral(6.564104, 6.422073) silane

# Integral(3.679897, 3.599926) hexene
# Integral(4.129027, 4.079348) hexene
# Integral(4.857709, 4.668059) hexene
# Integral(7.838666, 7.616258) hexene
# Integral(8.647265, 8.440374) hexene

# Integral(10.387961, 2.952301) all


reagent0 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[:5])
reagent1 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[5:8])
reagent2 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[8:])
reagents_sp = [reagent0, reagent1, reagent2]


if whole_mixture_or_regions == 'whole':
    
    mixture_time_data = cut_ppm_axis(mixture_time_data, FULL_PPM_INTERVAL_PMG_287)   
    
elif whole_mixture_or_regions == 'regions':
    
    mixture_time_data = cut_mixture_time_data_to_regions(mixture_time_data, INTEGRATION_INTERVALS)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_hexene.csv', np.array(reagent0.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_triethylsilane.csv', np.array(reagent1.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_product.csv', np.array(reagent2.confs), delimiter = '\t')

if whole_mixture_or_regions == 'whole':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_whole.csv', mixture_time_data, delimiter = '\t')
    
elif whole_mixture_or_regions == 'regions':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_cut_to_regions.csv', mixture_time_data, delimiter = '\t')


# PMG 287 MONITORING WITH FULL MIXTURE SPECTRUM


experiment_name = 'PMG 287 monitoring with full mixture spectrum'


MIXTURE_PATH = config["raw_mixture_paths"][experiment_name]
MIXTURE_SEPARATOR = config["mixture_separators"][experiment_name]
PREPROCESSING_OUTPUT_PATH = config["preprocessing_output_path"][experiment_name]
INTEGRATION_INTERVALS = config["integration_intervals"][experiment_name]
MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA = config["mixture_to_extract_reagents_spectra"][experiment_name]
FULL_PPM_INTERVAL_PMG_287 = config["full_ppm_interval_PMG_287"]


if experiment_name == 'PMG 287 monitoring':
    whole_mixture_or_regions = 'regions'
elif experiment_name == 'PMG 287 monitoring with full mixture spectrum':
    whole_mixture_or_regions = 'whole'


mixture_time_data = load_mixture_time_data(MIXTURE_PATH, MIXTURE_SEPARATOR, last_columns_to_skip=1)


# Integral(3.343730, 3.256018) product
# Integral(3.722376, 3.681636) product
# Integral(3.756886, 3.722376) product
# Integral(4.078496, 3.994139) product

# Integral(3.406997, 3.343730) silane
# Integral(3.832615, 3.756886) silane
# Integral(6.564104, 6.422073) silane

# Integral(3.679897, 3.599926) hexene
# Integral(4.129027, 4.079348) hexene
# Integral(4.857709, 4.668059) hexene
# Integral(7.838666, 7.616258) hexene
# Integral(8.647265, 8.440374) hexene

# Integral(10.387961, 2.952301) all


reagent0 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[:5])
reagent1 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[5:8])
reagent2 = create_reagent_spectrum(mixture_time_data, MIXTURE_TO_EXTRACT_REAGENTS_SPECTRA, INTEGRATION_INTERVALS[8:])
reagents_sp = [reagent0, reagent1, reagent2]


if whole_mixture_or_regions == 'whole':
    
    mixture_time_data = cut_ppm_axis(mixture_time_data, FULL_PPM_INTERVAL_PMG_287)   
    
elif whole_mixture_or_regions == 'regions':
    
    mixture_time_data = cut_mixture_time_data_to_regions(mixture_time_data, INTEGRATION_INTERVALS)


np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_hexene.csv', np.array(reagent0.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_triethylsilane.csv', np.array(reagent1.confs), delimiter = '\t')
np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_product.csv', np.array(reagent2.confs), delimiter = '\t')

if whole_mixture_or_regions == 'whole':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_whole.csv', mixture_time_data, delimiter = '\t')
    
elif whole_mixture_or_regions == 'regions':
    np.savetxt(PREPROCESSING_OUTPUT_PATH + 'preprocessed_mixture_cut_to_regions.csv', mixture_time_data, delimiter = '\t')


# INTEGRALS

experiment_name = 'Saccharose hydrolysis'

INTEGRALS_PATH = config["integrals_path"][experiment_name]
SUBSTANCES_NAMES = config["substances_names"][experiment_name]
RAW_INTEGRALS_SEPARATORS = config["raw_integrals_separators"][experiment_name]
OUTPUT_PATH = config["output_path"][experiment_name]

integrals = preprocess_mnova_integrals(INTEGRALS_PATH, RAW_INTEGRALS_SEPARATORS, SUBSTANCES_NAMES)

integrals.to_csv(OUTPUT_PATH + 'sacharoza_calki_nowe.csv', '\t', index=False)