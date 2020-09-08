# temm
R code for the temporal exceptional model mining paper

Simulated data are provided separately in case you want to use it with other algorithms. For using the TEMM algorithm, no data need to be downloaded as the script will automatically generate the data.

Two steps are needed to run simulations. 

1. Learn the exceptional subgroups from simulated data
This will generate the data based on several parameters, see ``sd_simuldata.R`` for more info. To proceed, run in R:
``source("sd_simuldata.R)``
This will create several files in your computer containing the results (R format).

2. Compute classification measures (AUROC, precision, recall, etc)
Type in R:
``source("sd_simul_processresults.R)``
Note that the parameters you set in the learning phase are required in this script as well.

The final results will be written to disk in two csv files.
