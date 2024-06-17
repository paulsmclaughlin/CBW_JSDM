#==================================================
# This script will run cross validation for 
# predictions on the cluster for one species 
# conditionally given another species or set of 
# species. This should be run after the 
# "modelFit_54Species.R" script has been run on 
# the cluster so the results from that can be used 
# without re-fitting the model for each species.
#==================================================

#get job ID from cluster
#jobid=16
jobid = as.integer(Sys.getenv("INPUT"))
set.seed(jobid)

#load GJAM and data
library(gjam)
library(InformationValue) #use package for choosing optimal cutoffs for predicting 1's and 0's
load("runData_54Species.RData")

#load in results from previous run so we don't need to re-run the model
outFileNumber <- rep(1:100,54)[jobid]
outIn <- readRDS(paste0('outFiles/out',outFileNumber,'.rds'))
out <- outIn$out
marginalCutoff<-out$marginalCutoffs
removeIndex <- outIn$removeIndex
ydat_pa <- ifelse(ydat>0,1,0)
xx <- orig_xdat[removeIndex,]
yy <- ydat_pa[removeIndex,]
eff <- list( columns = 1:ncol(yy), values = log(as.numeric(unlist(ef_stime)[removeIndex])))
newdata   <- list(xdata = xx,nsim = 2, efforts=eff )

#data.frame for conditional predictions
predYMu_cond <- data.frame(matrix(NA,nrow=nrow(yy),ncol=ncol(yy)))

#vector to store optimal cutoffs to make presence-absence predictions
#from probability predictions (there are multiple cuttoffs because
#we're considering different cutoffs for the varying amounts of
#conditional species given)
conditionalCutoffs <- rep(NA,length(givenSpecies))

#construct matrix of species ids sorted by highest to lowest residual correlations
speciesID <- rep(1:ncol(ydat_pa),each=100)[jobid]
sorted_corr <- sort(abs(out$parameters$corMu[,speciesID]),decreasing = T)
sorted_corr  <- sorted_corr[-1]#remove current species
givenSpecies <- match(names(sorted_corr),gsub(" ", "", colnames(ydat_pa)))



#----------------------------------------------------------------------------------
# Make conditional predictions for holdout set. This will first make conditional
# predictions for the larger fitting dataset to calculate optimal cutoffs to make
# presence-absence predictions from probability of presence predictions, then
# makes conditional probability of presence predictions given an increasing number
# of species to condition on.
#----------------------------------------------------------------------------------
#number of species given
ins <- c(1,2,3,4,5,6,7,8,9,10,20,30,53)
for(i in 1:length(ins)){
  
  #conditional predictions for fitting dataset (just predicting these to get optimal cutoffs)
  id <- givenSpecies[1:ins[i]]
  condDat <- as.data.frame(ydat_pa[-removeIndex,id])
  names(condDat) <- colnames(ydat_pa)[id]
  newdata   <- list(ydataCond=condDat, xdata = orig_xdat[-removeIndex,], nsim = 200, efforts=log(as.numeric(unlist(ef_stime)[-removeIndex])))
  preds <- gjamPredict(out, newdata = newdata,FULL=FALSE)
  fullCondPreds <- preds$sdList$yMu[,speciesID]
  
  #get cutoffs from conditional predictions for full dataset
  cond <- optimalCutoff(actuals=ydat_pa[-removeIndex,speciesID],
                        predictedScores=fullCondPreds,
                        optimiseFor="Both",returnDiagnostics=T)
  conditionalCutoffs[i] <- cond$optimalCutoff
  
  #make conditional predictions for out of sample data
  condDat <- as.data.frame(ydat_pa[removeIndex,id])
  names(condDat) <- colnames(ydat)[id]
  newdata   <- list(ydataCond=condDat, xdata = xx, nsim = 200, efforts=eff)
  preds <- gjamPredict(out, newdata = newdata,FULL=FALSE)
  predYMu_cond[,i] <- preds$sdList$yMu[,speciesID]
}

#save
saveRDS(list(preds=outIn$preds,condPreds=predYMu_cond,removeIndex=removeIndex,margCutoff=marginalCutoff,conCutoff=conditionalCutoffs,sortedCor=sorted_corr),file=paste0('preds/',colnames(ydat_pa)[speciesID],'_conPreds',jobid,'.rds'))

