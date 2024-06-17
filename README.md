# CBW_JSDM
Joint species distribution model for the Chesapeake Bay Watershed.

R scripts to run cross-validation analysis performed for the paper 54 species model presented in the paper:

McLaughlin, P., Krause, K. P., Maloney, K. O., Woods, T., & Wagner, T. (2024). Evaluating the effectiveness of joint species distribution modeling for fresh water fish communities within large watersheds. Canadian Journal of Fisheries and Aquatic Sciences, (ja).

-kFolds_groupAssigments.R: produces the “folds_54Species.rds” file which specifies which sites are in the holdout or fitting dataset during the cross validation analysis. The scipt will produce 100 different group assignments so the cross validation analysis can be run 100 times per species.

-fitModel_54Species.R: can be run on a cluster to fit the GJAM to the 100 different group assignments provided in the “folds_54Species.rds” file. The script will produce 100 output files (‘out(JOB ID).rds’) which are used during the cross validation analysis so the model will not have to be re-fit for each species. This script also performs the marginal predictions for the cross-validation analysis. 

-conPreds_54Species.R: can be run on a cluster 5400 times (100 times per species) to perform the cross-validation analysis for 54 species. The script will make conditional predictions given an increasing amount of conditional data and calculate optimal cutoffs for making presence absence predictions from the predicted presence probabilities.

-analyzeResults_54Species.R: analyzes results from the cross validation analysis obtained by running the “condPreds_54Species.R” script described above, and produces plots and tables presented in the main text and supplementary material.
