# Robust climate mitigation entails net negative emissions for centuries

Supporting material for Thomas Gasser, Armon Rezai, Côme Cheritel, Artem Baklanov, and Michael Obersteiner, “Robust climate mitigation entails net negative emissions for centuries.”

Choses importantes à signaler : 
- besoin d'une lcience gams
- besoin d'un cluster de calcul (temps de calculs pouvant être très long)
- ordres des codes
- Gams capricieux : parfois, il faut commenter et bien choisir le guess de départ qui doit être le plus proche possible de ce qui est utilisé
- description du solver utilisé (IPOPT)
- parler de OSCAR (lien vers le repo est le papier

# Description


## Structure

This repository includes the elements to reproduce the results presented in Gasser et al. "Robust climate mitigation entails net negative emissions for centuries.

This repository includes three different sections:
1. **Input data** from the Pathfinder model
2. **Codes** for the different robust simulation models whose results are presented
3. **Results tables** in excel spreadsheet format

## User suitability

Please note that the programming language used is GAMS (https://www.gams.com/). Since the size of the various models we run far exceeds the 1000 limit, a GAMS licence is required.

Furthermore note that running the simulation codes is extremely computationnally intensive since the model consists of the optimisation of about 20 state variables over 100 periods in 600 different states of the world. The power requirement can only be granted by a computing cluster.


## Folders

The folders in this repository are broadly consistent with the steps outlined above:

0_data_cleaning/ - Code for cleaning and constructing the dataset used to estimate the mortality-temperature relationship.

1_estimation/ - Code for estimating and plotting all mortality-temperature regression models present in the paper.

2_projection/ - Code for running future projections using Climate Impact Lab projection tools, and extracting, summarizing, and plotting the projection output.

3_valuation/ - Code for calculating the VSL based on various assumptions and applying those values to our projected impacts of climate change on mortality risk.

4_damage_function/ - Code for estimating empirical damage functions based upon monetized damages and GMST anomalies.

5_scc/ - Code for applying a CO2 pulse from the FAIR simple climate model to global damage functions, and summing damages over time to calculate mortality partial SCCs.

For run instructions on each step of the analysis, refer to the README files located within the corresponding directories.
Setup
Requirements For Using Code In This Repo

    You need to have python, Stata, and R programming capabilities, or at least environments to run code in these languages, on your computer.

    We use conda to manage python environments, so we recommend installing conda if you haven't already done so following these instructions.

## Known issues

    The model requires a high number of substeps to remain stable under high CO2 (because of the ocean carbon cycle). This can be set using the nt argument when calling run_xarray.

    The temperature-driven mode (Tdriven) is extremely sensitive to its forcings: it can be very difficult to make it transition smoothly from historical to projections. This is unavoidable because mathematically it requires the second derivative of T and the first derivative of ERFx as input.

    Unclear whether the my_AR1 class from cls_calib is actually needed.


# Setup Instructions

    Clone the following repos to a chosen directory, which we'll call yourREPO from now onwards, with the following commands:

cd <yourREPO>
git clone https://github.com/ClimateImpactLab/carleton_mortality_2022.git

    Install the conda environment included in this repo by running the following commands under the root of this repo:

cd <yourREPO>/carleton_mortality_2022
conda env create -f mortalityverse.yml

Try activating the environment:

conda activate mortalityverse

Please remember that you will need to activate this environment whenever you run python scripts in this repo, including the pip install -e . commands in the following section.

Also, you need to install Jupyter for the scc calculation code

conda install -c conda-forge jupyterlab

    Download data either from Zenodo or from the QJE Dataverse and unzip it somewhere on your machine with at least 85 GB of space. Let's call this location yourDATA.

    Set up a few environment variables so that all the code runs smoothly.

Append the following lines to your ~/.bash_profile.

First, run:

nano ~/.bash_profile

Then, point the variable DB in the yourDATA dierctory in the downloaded data, and do the same for OUTPUT. Point the REPO variable to yourREPO path used above containing this repo and other repos by adding the following lines to .bash_profile:

export REPO=<yourREPO>
export DB=<yourDATA>/data
export OUTPUT=<yourDATA>/output
export LOG=<yourDATA>/log

Save and exit. Then, run source ~/.bash_profile to load the changes we just made.

    Setup for the whole repo is complete! Please follow the READMEs in each subdirectory to run each part of the analysis. In general, each directory will contain one or more staging files where individual analysis or output producing scripts can be run from in one go. Before running, it is recommended that users review and set the TRUE/FALSE toggles to produce the desired set of outputs. More detail is available in the section READMEs.
