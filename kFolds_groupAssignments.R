#===========================================================
# This script randomly assigns sites to holdout or fitting
# dataset, and stratifies this assignment to ensure that
# no species has fewer than five observations in the 
# holdout dataset. The script will create random
# assignments for 100 cross validations runs.
#===========================================================


#Set working directory to parent dirctory of current R file is (where runData file is stored)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dirname(getwd()))

#load data and construct presence-absence data from abundance dataset
load('runData_54Species.RData')
ydat_pa <- ifelse(ydat>0,1,0)

#-----------------------------------------------------------
# randomly assign sites to hold out group or fitting group
# and ensure there's at least 5 observations for each species
# in the hold out set
#-----------------------------------------------------------
nRemove <- 112
removeIndex.df <- data.frame(matrix(NA,nrow=100,ncol=nRemove))
colSums.df <- data.frame(matrix(NA,nrow=100,ncol=54))
colSums_fit.df <- data.frame(matrix(NA,nrow=100,ncol=54))
nSamps <- 0
while(nSamps < 100){
  
  removeIndex <- sample(1:nrow(ydat_pa),nRemove,replace = F)
  
  #check that there's at least 5 observations for all the species
  if(TRUE %in% c(c(0:4) %in% colSums((ydat_pa[removeIndex,])))){
    #do nothing
  }else{
    removeIndex.df[nSamps+1,] <- removeIndex
    colSums.df[nSamps+1,] <- colSums(ydat_pa[removeIndex,])
    colSums_fit.df[nSamps+1,] <- colSums(ydat_pa[c(1:nrow(ydat_pa))[-removeIndex],])
    nSamps <- nSamps + 1
  }
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
saveRDS(removeIndex.df,"folds_54Species.rds")
