# Robust climate mitigation entails net negative emissions for centuries

Supporting material for Thomas Gasser, Armon Rezai, Côme Cheritel, Artem Baklanov, and Michael Obersteiner, “Robust climate mitigation entails net negative emissions for centuries.”

Choses importantes à signaler : 
- besoin d'une licence gams
- besoin d'un cluster de calcul (temps de calculs pouvant être très long)
- ordres des codes
- Gams capricieux : parfois, il faut commenter et bien choisir le guess de départ qui doit être le plus proche possible de ce qui est utilisé
- description du solver utilisé [CONOPT4](https://www.gams.com/latest/docs/S_CONOPT4.html)
- parler de pathfinder ([lien vers le repo](https://github.com/tgasser/Pathfinder) et [le papier](https://gmd.copernicus.org/articles/15/8831/2022/))


# Description


## Structure

This repository includes the elements to reproduce the results presented in Gasser et al. "Robust climate mitigation entails net negative emissions for centuries.

This repository includes three different sections:
1. **Input data** from the Pathfinder model
2. **Codes** for the different robust simulation models whose results are presented
3. **Results tables** in excel spreadsheet format

## User suitability

Please note that the programming language used is `GAMS` (https://www.gams.com/). Since the size of the various models we run far exceeds the 1000 limit, a GAMS licence is required.

Furthermore note that running the simulation codes is extremely computationnally intensive since the model consists of the optimisation of about 20 state variables over 100 periods in 600 different states of the world. The power requirement can only be granted by a computing cluster.


## Folders

The folders in this repository are broadly consistent with the steps outlined above:

`0_input_data/` - The several spreadsheets required to run the simulations.

`1_simulation/` - Code for estimating robust cost-benefit and cost-effective pathways.

`2_results/` - Excel spreadsheets containing the results of the GAMS simulation.


## Known issues

GAMS is a capricious software. Sometimes, with a given set of parameters, a model will be "infeasible", i.e. the feasibility step in the optimisation will not converge (the solver does not find a feasible solution from which to start the optimization). In this case, it may be relevant to exchange the rows of some states of the world with each other in the Excel input files (be careful to exchange the rows identically in the four input files, otherwise you change the underlying distribution of states of the world). More details in the following sections.

Similarly, if the initial guess is not close enough to the outcome, the model may not be able to find a feasible solution. This is particularly the case with some cost-effective simulations when switching from one probabilistic target to another. In this case, the easiest thing to do is:
1. to determine which simulation whose results have already been obtained appears to be the closest to the one sought
2. comment out model simulations with probabilistic targets preceding the target under consideration (i.e. the one whose result has already been released) using the `$ontext` and `$offtext` command lines
3. write the `resultsxxx.gdx` file corresponding to the simulation result selected in **1.** to the `$gdxin resultsxxx.gdx` and `execute_loadpoint "resultsxxx.gdx"` code lines
4. rerun the code and wait to see if the new guess results leads to a feasible solution. If the solution is not feasible, select another simulation result file close to the one sought and go back to step **3.**.

# How-to

Requirements for using the code in this repo

    You need to have `GAMS` as well as a `GAMS` licence too.

## Basic instructions

First you should download the input data files (folder `0_input_data/`) and the simulation GAMS code files (folder `1_simulation/`) and put all files in the same folder.

You may open the GAMS IDE project file `Robust_Pathfinder_ncc_parameters_log_600.gpr` with the GAMS IDE, and open the GAMS file `Robust_Pathfinder_ncc_parameters_log_600.gms`.

This first file performs:
- the cost-benefit Monte Carlo simulation,
- the cost-benefit robust simulation,
- the cost-effective robust probabilistic simulations with a constraint limiting warming to 2 degrees with a probability of 95%, 90%, 60%, 50%, 33% and 10%,
- the cost-effective robust probabilistic simulations with a constraint limiting warming to 1.5 degrees with a probability of 95%, 90%, 60%, 50%, 33% and 10% (note that the 95%, 90% and 66% trajectories are not achievable),
- the cost-effective Monte Carlo simulation with a constraint limiting warming to 2 degrees,
- the cost-effective Monte Carlo simulation with a constraint limiting warming to 1.5 degrees.

Once these simulations have run, it is possible to run the codes corresponding to the probabilistic targets and the Monte Carlo simulations.
Here is a list of the simulations performed:
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_03.gms` estimates the optimal trajectories to limit the permafrost thawing to 30% with a probability of 95%, 90%, 60%, 50%, 33% and 10% (note that the 95% and 90% trajectories are not achievable)
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_04.gms` estimates the optimal trajectories to limit the permafrost thawing to 40% with a probability of 95%, 90%, 60%, 50%, 33% and 10% (note that the 95% trajectory is not achievable)
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_05.gms` estimates the optimal trajectories to limit the permafrost thawing to 50% with a probability of 95%, 90%, 60%, 50%, 33% and 10%
- `Robust_Pathfinder_ncc_parameters_log_600_pH_02.gms` estimates the optimal trajectories to limit the ocean acidification to 0.2 pH points compared to the preindustrial era equilibrium with a probability of 95%, 90%, 60%, 50%, 33% and 10%
- `Robust_Pathfinder_ncc_parameters_log_600_SLR_dot_50.gms` estimates the optimal trajectories to limit the sea-level rise speed to 5.0 mm.yr-1 with a probability of 95%, 90%, 60%, 50%, 33% and 10%
- `Robust_Pathfinder_ncc_parameters_log_600_SLR_dot_75.gms` estimates the optimal trajectories to limit the sea-level rise speed to 7.5 mm.yr-1 with a probability of 95%, 90%, 60%, 50%, 33% and 10%
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_03_MC.gms` estimates the optimal Monte Carlo trajectories limiting the permafrost thawing to 30%.
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_04_MC.gms` estimates the optimal Monte Carlo trajectories limiting the permafrost thawing to 40%.
- `Robust_Pathfinder_ncc_parameters_log_600_permafrost_05_MC.gms` estimates the optimal Monte Carlo trajectories limiting the permafrost thawing to 50%.
- `Robust_Pathfinder_ncc_parameters_log_600_pH_02_MC.gms` estimates the optimal Monte Carlo trajectories limiting the ocean acidification to 0.2 pH points compared to the preindustrial era equilibrium.
- `Robust_Pathfinder_ncc_parameters_log_600_SLR_dot_50_MC.gms` estimates the optimal Monte Carlo trajectories limiting the sea-level rise speed to 5.0 mm.yr-1.
- `Robust_Pathfinder_ncc_parameters_log_600_SLR_dot_75_MC.gms` estimates the optimal Monte Carlo trajectories limiting the sea-level rise speed to 7.5 mm.yr-1.

## Input data

Each simulation requires a set of parameters for each state of the world. Each sets of parameters are distributed over 4 files produced using the [Pathfinder model](https://doi.org/10.5194/gmd-15-8831-2022):
- `par_v1` denotes the set of parameters of the climate model in each state of the world,
- `ini_2015_v1` denotes the initial value of the state variables for each state of the world,
- `ERFx_for_dice` denotes the exogenous radiative forcing scenario for each state of the world,
- `Eluc_for_dice` denotes the land-use change emissions scenario for each state of the world.

Original outputs of the [Pathfinder model](https://doi.org/10.5194/gmd-15-8831-2022) are given by the files `par_v1_original.xlsx`, `ini_2015_v1_original.xlsx`, `ERFx_for_dice.xlsx`and `Eluc_for_dice.xlsx`. Morover, for matters of practicalities, files `ERFx_for_dice.xlsx` and `Eluc_for_dice.xlsx` need to be transposed as files `ERFx_for_dice_2_original.xlsx` and `Eluc_for_dice_2_original.xlsx`. As explained above, GAMS might be a capricious software. In order to obtain a feasible model, some lines of the four models are interchanged (in the same way between the four files) which allows to obtain the files `par_v1.xlsx`, `ini_2015_v1.xlsx`, `ERFx_for_dice_2.xlsx`and `Eluc_for_dice_2.xlsx` that are used as input for the simulation file `Robust_Pathfinder_ncc_parameters_log_600.gms` as well as for the cost-effective probabilistic targets simulations files. Regarding the Monte Carlo simulation files they require again some interchanging of lines to obtain a feasible model, so the input files are:
-  `par_v1_03A_MC.xlsx`, `ini_2015_v1_03A_MC.xlsx`, `ERFx_for_dice_2_03A_MC.xlsx`and `Eluc_for_dice_2_03A_MC.xlsx for the simulation file `Robust_Pathfinder_ncc_parameters_log_600_permafrost_03_MC.gms`,
-  `par_v1_04A_MC.xlsx`, `ini_2015_v1_04A_MC.xlsx`, `ERFx_for_dice_2_04A_MC.xlsx`and `Eluc_for_dice_2_04A_MC.xlsx for the simulation file `Robust_Pathfinder_ncc_parameters_log_600_permafrost_04_MC.gms`,
-  `par_v1.xlsx`, `ini_2015_v1.xlsx`, `ERFx_for_dice_2.xlsx`and `Eluc_for_dice_2.xlsx
-  `par_v1.xlsx`, `ini_2015_v1.xlsx`, `ERFx_for_dice_2.xlsx`and `Eluc_for_dice_2.xlsx
-  `par_v1.xlsx`, `ini_2015_v1.xlsx`, `ERFx_for_dice_2.xlsx`and `Eluc_for_dice_2.xlsx
-  `par_v1.xlsx`, `ini_2015_v1.xlsx`, `ERFx_for_dice_2.xlsx`and `Eluc_for_dice_2.xlsx


# References
    
The Pathfinder model is described in Bossy, T., Gasser, T., and Ciais, P.: Pathfinder v1.0.1: a Bayesian-inferred simple carbon–climate model to explore climate change scenarios, Geosci. Model Dev., 15, 8831–8868, https://doi.org/10.5194/gmd-15-8831-2022, 2022. 
The repository of the code is available [here](https://github.com/tgasser/Pathfinder).

The GAMS solver used is a non-linear optimization solver named [CONOPT4](https://www.gams.com/latest/docs/S_CONOPT4.html)
    
