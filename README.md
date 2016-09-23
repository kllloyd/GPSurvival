# GPSurvival
Code from *Gaussian processes for survival data with right censoring*.

This code generates/extracts the data and plots for the paper. 
All code is written in R.

To run, all files are required.
Each plot has an associated script:
-runSyntheticExp1
-runSyntheticExp2
-runSyntheticExp3
-runRealYuanEtAl
-runRealTothillEtAl

Upon completion, each script will save results into a folder named "Runs" in the working directory, which is required to already exist.

Required libraries are:
-curatedOvarianData
-fields
-foreach
-gbm
-glmnet
-impute
-ipred
-MASS
-Matrix
-mice
-nlme
-NORMT3
-pdist
-randomForestSRC
-rgl
-rms
-survcomp
-survival
-synapseClient
-zoo
