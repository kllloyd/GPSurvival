# GPSurvival
Code from *Gaussian processes for survival data with right censoring*.

This code generates/extracts the data and plots for the paper. 
All code is written in R.

To run, all files are required.
Each plot has an associated script:
- runSyntheticExp1.R
- runSyntheticExp2.R
- runSyntheticExp3.R
- runRealYuanEtAl.R
- runRealTothillEtAl.R

Upon completion, each script will save results into a folder named "Runs" in the working directory, which is required to already exist.

Both Gaussian process regression and the Gaussian process for survival data models, GPS1, GPS2 and GPS3, may be run using the ApplyGP.R function, as applied in the scripts above.

For the real data experiments, data and/or code were used as supplied by the original authors where possible. 
For the Yuan et al. (2014) data, code was obtained from Synapse, synapse id:syn1720423. Data is extracted from the same source using the Synapse R client. This requires a login.
The Tothill et al. (2008) data was extracted using the R package curatedOvarianData, as the data set GSE9891.

Required libraries are:
- curatedOvarianData
- fields
- foreach
- gbm
- glmnet
- impute
- ipred
- MASS
- Matrix
- mice
- nlme
- NORMT3
- pdist
- randomForestSRC
- rgl
- rms
- survcomp
- survival
- synapseClient
- zoo
