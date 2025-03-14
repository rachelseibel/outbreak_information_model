# Outbreak Information Model

This repository contains the code for the outbreak information model (OIM) described in the pre-print "A model for outbreak information management in public health" by [R. Seibel](https://warwick.ac.uk/fac/sci/mathsys/people/students/mathsysii/seibel/), [M. Tildesley](https://warwick.ac.uk/fac/sci/lifesci/people/mtildesley/) and [E. Hill](https://edmhill.github.io) (2024). The model is implemented in Python 3.11.7 and can be used to deterministically simulate the spread of a human infectious disease in a population given real-time awareness of intervention effectiveness. The model is described in detail in the pre-print, which can be found [on medRxiv](https://www.medrxiv.org/content/10.1101/2024.01.17.24301344v1).

Archived code associated with this study is available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15024878.svg)] [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10519963.svg)](https://doi.org/10.5281/zenodo.10519963)

## To run the model:
1. Clone the repository
2. Install the required packages (pip install using requirements.txt or activate conda environment from cloned repository)
3. Determine which population structure you want to use (homogeneous or heterogeneous) and which scenario (1 or 2) you wish to run. Note that some scenarios are separated into multiple files (e.g. homogeneous scenario 1 is split into 1a and 1b).
4. Open the corresponding model file (e.g. src/main_homogeneous.py) and set the file path for the parameter variable (e.g. "results/homogeneous/inputs/parameters_scenario1a.csv")
5. Execute the cells in the corresponding Jupyter notebook (e.g. Homogeneous.ipynb). Save the results to the "outputs" folder for the corresponding population structure and scenario.

## To interpret the results:
The model outputs a number of files, which are saved in the "outputs" folder for the correponding population structure and scenario. The files are:
* `results.csv`: Contains the cumulative results of the model run. Each row corresponds to a single simulation and contains the following columns:
    * `pathogen`: The pathogen used in the simulation
    * `split`: For the heterogeneous simulations only, the split between the three subpopulations
    * `memory_window`: The memory window used in the simulation
    * `behaviour_function`: The behaviour function used in the simulation
    * `r`: The vaccine opinion(s) used in the simulation
    * `alpha`: The information sensitivity used in the simulation
    * `vaccine_efficacy`: The vaccine efficacy (as a proportion) used in the simulation
    * `final_cases`: The number of cases at the end of the simulation
    * `final_deaths`: The number of deaths at the end of the simulation
    * `final_vaccinated`: The number of vaccinated individuals at the end of the simulation
    * `epidemic_duration`: The duration of the epidemic in days
* `temporal.csv`: Contains the temporal results of the model run. Each row corresponds to a single day in a single simulation and contains the following columns:
    * `pathogen`: The pathogen used in the simulation
    * `split`: For the heterogeneous simulations only, the split between the three subpopulations
    * `memory_window`: The memory window used in the simulation
    * `behaviour_function`: The behaviour function used in the simulation
    * `r`: The vaccine opinion(s) used in the simulation
    * `alpha`: The information sensitivity used in the simulation
    * `vaccine_efficacy`: The vaccine efficacy (as a proportion) used in the simulation
    * `t`: The day of the simulation
    * `S`: The number of susceptible, unvaccinated individuals on day t
    * `E`: The number of exposed, unvaccinated individuals on day t
    * `I`: The number of infectious, unvaccinated individuals on day t
    * `R`: The number of recovered, unvaccinated individuals on day t
    * `PD`: The number of pre-death, unvaccinated individuals on day t
    * `Sv`: The number of susceptible, vaccinated individuals on day t
    * `Ev`: The number of exposed, vaccinated individuals on day t
    * `Iv`: The number of infectious, vaccinated individuals on day t
    * `Rv`: The number of recovered, vaccinated individuals on day t
    * `PDv`: The number of pre-death, vaccinated individuals on day t
    * `C`: The number of cumulative cases in unvaccinated individuals on day t
    * `Cv`: The number of cumulative cases in vaccinated individuals on day t
