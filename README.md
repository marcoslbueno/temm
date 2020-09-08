# temm
R code for the paper "Temporal Exceptional Model Mining using Dynamic Bayesian Networks", published at AALTD workshop at ECML 2020. For more details, see: https://project.inria.fr/aaltd20/files/2020/08/AALTD_20_paper_Bueno.pdf

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

The final results will be written to disk in two csv files, one for the unitary subgroups and another one for the specialized subgroups. Examples:

``results-dsimul-unitary-ndynamic17-numseq10.csv``

``results-dsimul-specialized-ndynamic17-numseq10.csv``
