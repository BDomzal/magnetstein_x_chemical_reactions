# Chemical reactions analysis using Magnetstein
![main_workflow_final](https://github.com/BDomzal/magnetstein_x_chemical_reactions/blob/main/main_figure.png)

This repository contains the code for the analysis of chemical reaction kinetics using nuclear magnetic resonance (NMR) spectra and the Magnetstein algorithm. 

Magnetstein is available [here](https://github.com/BDomzal/magnetstein). To be able to use all the functionalities dedicated for the analysis of chemical reactions, make sure that you set branch to initial_proportions.

## About the method

#### When to use Magnetstein instead of traditional integration methods

The method should be used when:

- peaks move throughout the reaction,
- peaks from different reagents overlap,
- lineshapes are distorted,
- there are contaminations and/or noise.
  
In such circumstances, the Magnetstein algorithm gains an advantage over the traditional integration methods due to the special properties of the Wasserstein metric that is the core concept of the method. The metric makes the algorithm robust to changes in peaks' locations and shapes, and to ambiguity in assigning signal to particular reagents due to overlap. The additional refinements make it possible for the algorithm to remove noise from the data.

#### What to use as an input

The input should consist of the two crucial parts:

- mixture, i.e. a series of NMR spectra of reaction mixture measured in the consecutive moments of time,
- library, i.e. a set of spectra of individual reagents expected to be present in the mixture (no need for a series here, just a single spectrum for each reagent).

#### How to interpret and set the values of the parameters

The user needs to define the values of two parameters: $\kappa_{mixture}$ and $\kappa_{components}$. These are so-called *denoising penalties* that can be interpreted as soft tolerance thresholds for shifting of the signal along the horizontal axis for the spectrum of the reaction mixture and for the spectra in the library, respectively. Another interpretation is viewing $\kappa_{mixture}$ as a certain measure of reliability of the mixture's spectrum. The higher its value, the less likely the algorithm is to remove noise from the spectrum. Similarly, $\kappa_{components}$ reflects our confidence in the purity of the spectra in the library.

If you are unsure about how to set the parameters, we recommend using the default values, i.e. $\kappa_{mixture}=0.25$ and $\kappa_{components}=0.22$. We checked experimentally that such settings produced accurate results for many datasets.

#### How to interpret the output

The main output of the algorithm is the series of vectors indexed by time:

$$p_t = \Big(p_{1,t}, p_{2,t}, \dots, p_{k,t}\Big),$$

where $p_{i,t}$ is the *proportion* of the $i$-th reagent in the reaction mixture. By *proportions* here we mean the relative amounts of reagents. Note that $p_{1,t} + p_{2,t} + \dots + p_{k,t}$ does not necessarily equal to 1 due to the presence of noise and contamination. The quantity $p_{0,t} := 1 - p_{1,t} - p_{2,t} - \dots - p_{k,t}$ is the Magnetstein's estimation of the relative amount of the signal coming from the contamination in the reaction mixture's spectrum. Similarly, the quantity $p_{0,t}'$, also returned by the algorithm, is the Magnestein's estimation of the relative amount of the signal coming from the contamination in the library.

#### How to report the results

In order to make the results reproducable, one needs to provide: 1) input data, i.e. series of spectra of the reaction mixture indexed by time, and the library; 2) information about the parameters settings, i.e. chosen values of $\kappa_{mixture}$ and $\kappa_{components}$. This is enough to rerun the analysis.

## About the repository

#### Prerequisites

To be able to run the code from this repository, set your environment as in requirements.txt file. We also strongly recommend installing [Gurobi](https://www.gurobi.com/) and using it as a solver in the analysis.

#### Data

Due to the large size of the input data, they are not stored in this repository. The data are available [here](https://zenodo.org/records/14814657). The data folder should be stored as magnetstein_x_chemical_reactions/data to ensure that all the paths are correct.

#### Code 

The code used for preprocessing, analysis and visualisations is available in notebooks/. To reproduce the analysis (without warm-start), use notebooks/estimation.ipynb. To run the analysis with warm-start, use notebooks/estimation_with_warm_start.ipynb.

#### Results

The numerical results as well as figures are available in results/. Additional figures showing the comparison of results for different settings of parameters are available in kappa_panels/.

#### Comparison with other tools

To compare the results from Magnetstein with other tools (Mnova and manual integration in Python), we use the data stored in mnova_integrals/ and python_integrals/.

## Citing 

If you use tools from this package, please cite:

Domżał, B., Nawrocka, E.K., Gołowicz, D., Ciach, M.A., Miasojedow, B., Kazimierczuk, K., & Gambin, A. (2023). Magnetstein: An Open-Source Tool for Quantitative NMR Mixture Analysis Robust to Low Resolution, Distorted Lineshapes, and Peak Shifts. _Analytical Chemistry_. DOI: [10.1021/acs.analchem.3c03594](https://doi.org/10.1021/acs.analchem.3c03594).
