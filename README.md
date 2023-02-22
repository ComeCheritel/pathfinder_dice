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
2. comment on model simulations with probabilistic targets previous to the code under consideration using the `$ontext' and `$offtext' command lines
3. load the `resultsxxx.gdx` file corresponding to the simulation result selected in **1.** to the `$gdxin resultsxxx.gdx` and `execute_loadpoint "resultsxxx.gdx"` code lines

De la même manière, si le guess au départ n'est pas assez proche du résultat, il se peut que le modèle ne parvienne pas à trouver de solution feasible. C'est en particulier ce qui arrive dans le cas de certaines simulations cost-effective lors du passage d'une cible probabilisée à une autre. Dans ce cas, le plus simple est:
1. de déterminer quelle simulation dont les résultats sont déjà obtenus apparait comme la plus proche de celle recherchée
2. de commenter les simuations de modèle avec des cibles probabilistes prcédentes du code considéré en utilisant les lisges de commande `$ontext` et `$offtext`
3. de charger le fichier `resultsxxx.gdx` correspondant au résultat de simulation selectionnné en **1.** aux lignes du code `$gdxin resultsxxx.gdx` et `execute_loadpoint "resultsxxx.gdx"  ;`

# How-to

Requirements for using the code in this repo

    You need to have `GAMS` as well as a `GAMS` licence too.

## Basic instruction

First you should download the input data files (folder `0_input_data/`) and the simulation GAMS code files (folder `1_simulation/`) and put all files in the same folder.

You may open the GAMS IDE project file `Robust_Pathfinder_ncc_parameters_log_600.gpr` with the GAMS IDE, and open the GAMS file `Robust_Pathfinder_ncc_parameters_log_600.gms`.

This first file performs:
- the cost-benefit Monte Carlo simulation,
- the cost-benefit robust simulation,
- the cost-effective robust probabilistic simulations with a constraint limiting warming to 2 degrees,
- the cost-effective robust probabilistic simulations with a constraint limiting warming to 1.5 degrees,
- the cost-effective Monte Carlo simulation with a constraint limiting warming to 2 degrees,
- the cost-effective Monte Carlo simulation with a constraint limiting warming to 1.5 degrees.

Once this simulation has run

# References
    
The Pathfinder model is described in Bossy, T., Gasser, T., and Ciais, P.: Pathfinder v1.0.1: a Bayesian-inferred simple carbon–climate model to explore climate change scenarios, Geosci. Model Dev., 15, 8831–8868, https://doi.org/10.5194/gmd-15-8831-2022, 2022. 
The repository of the code is available [here](https://github.com/tgasser/Pathfinder).

The GAMS solver used is a non-linear optimization solver named [CONOPT4](https://www.gams.com/latest/docs/S_CONOPT4.html)
    
