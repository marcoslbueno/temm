# temm
R code for the temporal exceptional model mining paper

Simulated data are provided separately in case you want to use it with other algorithms. For using the TEMM algorithm, no data need to be downloaded as the script will automatically generate the data.

Two steps are needed to run simulations. 

1. Learn the exceptional subgroups from simulated data

This will generate the data based on several parameters, see ``sd_simuldata.R`` for more info. 

To proceed, type in the R console:

``source("sd_simuldata.R")``

This will create several files in your computer containing the results (R format).

2. Compute classification measures (AUROC, precision, recall, etc)

The subgroups found will be evaluated now. Type in the console:

``source("sd_simul_processresults.R")``

Note that the parameters you set in the learning phase are required in this script as well.

The final results will be written to disk in two csv files, one for the unitary subgroups and another one for the specialized subgroups.
