# Chemical reactions analysis using Magnetstein
![main_workflow_final](https://github.com/BDomzal/magnetstein_x_chemical_reactions/blob/main/main_figure.png)

This repository contains the code for the analysis of chemical reaction kinetics using nuclear magnetic resonance (NMR) spectra and the Magnetstein algorithm. 

Magnetstein is available [here](https://github.com/BDomzal/magnetstein). To be able to use all the functionalities dedicated for the analysis of chemical reactions, make sure that you set branch to initial_proportions.

# Data

Due to the large size of the input data, they are not stored in this repository. The data are available [here](https://zenodo.org/records/14814657). The data folder should be stored as magnetstein_x_chemical_reactions/data to ensure that all the paths are correct.

# Code 

The code used for preprocessing, analysis and visualisations is available in notebooks/. To reproduce the analysis (without warm-start), use notebooks/estimation.ipynb. To run the analysis with warm-start, use notebooks/estimation_with_warm_start.ipynb.

# Results

The numerical results as well as figures are available in results/. Additional figures showing the comparison of results for different settings of parameters are available in kappa_panels/.

# Comparison with other tools

To compare the results from Magnetstein with other tools (Mnova and manual integration in Python), we use the data stored in mnova_integrals/ and python_integrals/.

# Citing 

If you use tools from this package, please cite:

Domżał, B., Nawrocka, E.K., Gołowicz, D., Ciach, M.A., Miasojedow, B., Kazimierczuk, K., & Gambin, A. (2023). Magnetstein: An Open-Source Tool for Quantitative NMR Mixture Analysis Robust to Low Resolution, Distorted Lineshapes, and Peak Shifts. _Analytical Chemistry_. DOI: [10.1021/acs.analchem.3c03594](https://doi.org/10.1021/acs.analchem.3c03594).
