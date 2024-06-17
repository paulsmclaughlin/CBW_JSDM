#======================================================
# This script performs one run of the cross validation
# analysis, and thus can be run multiple times on
# a cluster. The script does the following:
#
#  1) GJAM model will be fit to 90% of the data,
#     and hold out 10% of the data for testing.
#
#  2) Marginal probability predictions will be 
#     made for the held out 10% of the data for a post 
#     script cross validation analysis performed in the
#     analyze_crossVal.R script using the cutoffs 
#     calculated in step (2).
#
#  3) Optimal cutoffs for converting probability
#     predictions into 1's and 0's will be
#     calculated from marginal predictions made 
#     for the 90% dataset.
#===================================================

#Get job ID from cluster
jobid = as.integer(Sys.getenv("INPUT"))
set.seed(jobid)

#Load data and packages 
library(gjam)
library(InformationValue) #package to choose optimal cutoffs for 1's and 0's
load("runData_54Species.RData")

#convert abundance data to presence-absence
ydat_pa <- ifelse(ydat>0,1,0)

#choose which fold current run is in to choose hold out and fitting data
removeIndex.df <- readRDS("folds_54Species_randomStrat.rds")
removeIndex.df.index <- rep(1:100,54)[jobid]
removeIndex <- as.numeric(removeIndex.df[removeIndex.df.index,])#sample(1:nrow(ydat),112,replace = F)#folds.read[[run.jobid]][[fold.jobid]]
xdat_run <- orig_xdat[-removeIndex,]
ydat_pa_run <- ydat_pa[-removeIndex,]
ef_run <- list(columns = 1:ncol(ydat), values = log(as.numeric(unlist(ef_stime)))[-removeIndex])
ml <-list(FULL=F, PREDICTX = T, ng=5000,burnin=1000,typeNames=rep("PA",ncol(ydat_pa_run)),
          REDUCT=F,effort=ef_run,PLOTALLY = F)

#make HUC8 random
ml$random <- 'HUC8'


#--------------------------------------------------------------
# Fit GJAM model
#--------------------------------------------------------------
out = gjam(~ TotDASqKM + DevAll08 + AgrAll08 + DamNrmStorWs +
             TOTSTREAMSLOPE + TOTOLSONCAO + PrecipWs +
             TmeanWs + TOTFRESHWATERWD +
             TOTBFI + TOTTWI, xdata=xdat_run, ydata=ydat_pa_run,
           modelList=ml)


#Calculate optimal marginal cutoffs to use for post script cross validation analysis
margCuts <- rep(NA,ncol(ydat))
for(i in 1:ncol(ydat)){
  marg <- optimalCutoff(actuals=ydat_pa_run[,i],
                        predictedScores=out$prediction$ypredMu[,i],
                        optimiseFor='Both',returnDiagnostics=T)
  margCuts[i] <- marg$optimalCutoff
}


#---------------------------------------------------------------
#Out-of-sample unconditional (marginal) probability predictions
#(1's and 0's predictions made later in post script analysis)
#---------------------------------------------------------------
xx <- orig_xdat[removeIndex,]
yy <- ydat_pa[removeIndex,]
eff <- list( columns = 1:ncol(yy), values = log(as.numeric(unlist(ef_stime)[removeIndex])))
newdata   <- list(xdata = xx,nsim = 200, efforts=eff )
preds <- gjamPredict(out, newdata = newdata,FULL=FALSE)
predYMu <- preds$sdList$yMu

#save model run and marginal predictios
saveRDS(list(out=out,preds=predYMu,removeIndex=removeIndex,marginalCutoffs=margCuts),file=paste0('out',jobid,'.rds'))
