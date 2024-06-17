#===================================================
# this script will analyze results from running the
# conPreds_54Species.R script on the cluster, and
# produce plots and tables given in the main text
# and supplementary material.
#===================================================

#================================================
# load data and packages
#================================================
#Set working directory to parent dirctory of current R file is (where runData file is stored)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dirname(getwd()))

#load data and construct presence-absence data from abundance dataset
load('runData_54Species.RData')
f <- list(ydata=ifelse(ydat>0,1,0))

library(InformationValue) #used for picking optimal cutoff for 1's vs. 0's
library(gjam)
library(pROC) #for AUC calculation
library(vcd) #for Cohen's Kappa (K) calculation (Kappa function, not base kappa function)
library(ggplot2)
library(ggrepel)#to have non-overlapping labels in ggplot

#convert abundance data to presence-absence
ydat_pa <- ifelse(ydat>0,1,0)

#lower case names for plots (names stored from the model are capitalized)
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub',
               'Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout',
               'Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub',
               'Cutlips minnow','East. blacknose dace','East. creek chubsuck.','Eastern mosquitofish',
               'Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish',
               'Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace',
               'Margined madtom','Mount. redbelly dace','North. hog sucker','Pirate perch',
               'Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel',
               'River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner',
               'Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner',
               'Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')


#initial matrices and lists for storing results
AUC <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
RMSE <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
K <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
AUC_c <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
RMSE_c <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
K_c <- matrix(NA,nrow=100,ncol=ncol(f$ydata))
AUC_c.list <- list()
RMSE_c.list <- list()
K_c.list <- list()


setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/results'))
dataFiles_all <- list.files()



#=================================================================
# function from stackoverflow comment by AGS (September 28, 2013)
# https://stackoverflow.com/questions/4903092/calculate-auc-in-r
#=================================================================
getROC_AUC = function(probs, true_Y){
  probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
  val = unlist(probsSort$x)
  idx = unlist(probsSort$ix)  
  
  roc_y = true_Y[idx];
  stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
  stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    
  
  auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
  return(auc=auc)
}


#-----------------------------------------------------------------------
# function to calculate TP, FP, TN, FN to construct a confusion table
# to calculate Cohen's Kappa
#-----------------------------------------------------------------------
make_confusionMatrix <- function(actual,pred){
  result <- ifelse(actual==pred & actual==1,'tp',
              ifelse(actual!=pred & actual==0,'fp',
                 ifelse(actual!=pred & actual==1, 'fn',
                    ifelse(actual==pred & actual==0, 'tn', 'error'))))
  res.tbl <- table(result)
  
  return(matrix(c(res.tbl[4],
                  res.tbl[2],
                  res.tbl[1],
                  res.tbl[3]),nrow=2,ncol=2))
}


#======================================================================
# Load in and store results, and calculate AUC, RMSE and Cohen's Kappa
#======================================================================
nSpecies <- ncol(ydat)
auc.df <- data.frame(matrix(NA,nrow=nSpecies,ncol=2))
dataFileLengths <- data.frame()
nObs <- data.frame(matrix(NA,ncol=nSpecies,nrow=100))
ins <- c(1,2,3,4,5,6,7,8,9,10,20,30,nSpecies-1)
unSorted_cors.list <- list()
for(speciesID in c(1:nSpecies)){
  
  print(speciesID)
  
  dataFiles <- dataFiles_all[which(gsub("\\_.*","",dataFiles_all)==colnames(ydat)[speciesID])]
  unsorted_cors.df <- data.frame(matrix(NA,ncol=nSpecies-1,nrow=length(dataFiles)))
  
  if(length(dataFiles)==0){
    #do nothing
  }else{
    
    #==================================
    # Predictions for cross validation
    #==================================
    margCuts <- rep(NA,length(dataFiles))
    nZero <- 0
    nLessThan5 <- 0
    for(i in 1:length(dataFiles)){
      
      preds <- readRDS(dataFiles[i])
      predYMu.i <- preds$preds[,speciesID]
      predYMu_cond.i <- preds$condPreds
      removeIndex.i <- preds$removeIndex
      margCut <- preds$margCutoff
      conCut <- preds$conCutoff
      
      margCuts[i] <- margCut
      
      #re-order sorted correlation matrix to match the order
      #of the names of y_dat and store in a matrix so we can
      #find largest average correlation for each species
      sortedCor <- preds$sortedCor
      sortedCor_nameOrder <- match(gsub(" ","",colnames(ydat)),names(sortedCor))
      unsorted_cors.df[i,] <- na.omit(sortedCor[sortedCor_nameOrder])
      unsorted_names <- names(na.omit(sortedCor[sortedCor_nameOrder]))
      
      #marginal AUC
      actual <- f$ydata[removeIndex.i,speciesID]
      
      #store the number of observations for current dataset and species
      nObs[i,speciesID] <- sum(actual)
      
      #keep track of the numbers of observed data with less than 5
      #observations or zero observations
      if(sum(actual)<5){
        nLessThan5 <- nLessThan5 + 1
        
        if((sum(actual)==0)){
          nZero <- nZero + 1
        }
      }
      
      if(1 %in% actual){
        
        #roc_object_marg <- roc(actual, predYMu.i)
        AUC[i,speciesID] <- getROC_AUC(predYMu.i,actual) #auc(roc_object_marg)

        preds.i <- ifelse(predYMu.i > margCut, 1, 0)
        confusMat <- make_confusionMatrix(actual, preds.i)
        K[i,speciesID] <- Kappa(confusMat)$Unweighted[1]
        RMSE[i,speciesID] <- sqrt(mean((na.omit(preds.i-actual))^2)) #auc(roc_object_marg)
        
        
        for(j in 1:length(ins)){
          
          #conditional statistics
          actual <- f$ydata[removeIndex.i,speciesID]
          AUC_c[i,j] <- getROC_AUC(predYMu_cond.i[,j],actual)
          
          preds.ij <- ifelse(predYMu_cond.i[,j] > conCut[j], 1, 0)
          confusMat <- make_confusionMatrix(actual, preds.ij)
          K_c[i,j] <- Kappa(confusMat)$Unweighted[1]
          RMSE_c[i,j] <- sqrt(mean((na.omit(preds.ij-actual))^2))
        }
      }
      
    }
    #put everything into a data.frame
    dataFileLengths <- rbind(dataFileLengths,data.frame(species=colnames(ydat)[speciesID],
                                                        speciesID=speciesID,
                                                        length = length(dataFiles),
                                                        nZero = nZero,
                                                        nLessThan5=nLessThan5))
  }
  
  #store results for current species in result lists
  AUC_c.list[[speciesID]] <- AUC_c[,1:length(ins)]
  RMSE_c.list[[speciesID]] <- RMSE_c[,1:length(ins)]
  K_c.list[[speciesID]] <- K_c[,1:length(ins)]
  
  colnames(unsorted_cors.df) <- unsorted_names 
  unSorted_cors.list[[speciesID]] <- unsorted_cors.df
}


#determine species with highest average residual correlation
corNames <- matrix(nrow=nSpecies,ncol=nSpecies-1)
corVals <- matrix(nrow=nSpecies,ncol=nSpecies-1)
corIndexes <- matrix(nrow=nSpecies,ncol=nSpecies-1)
for(i in 1:nSpecies){
  corrMat.i <- unSorted_cors.list[[i]]
  colMeans.i <- sort(colMeans(corrMat.i),decreasing = T)
  corNames[i,] <- names(colMeans.i)
  corVals[i,] <- as.numeric(colMeans.i)
  corIndexes[i,] <- match(corNames[i,],gsub(" ","",colnames(ydat)))
}



#---------------------------------------
# Manipulate AUC, RMSE and Kappa plots
#---------------------------------------
library(manipulate)
mf <- function(i,metric){
  speciesID <- i
  par(mfrow=c(1,1))
  clrs <- c('black',rep('white',59))
  
  if(metric=='AUC'){
    boxDat <- cbind(AUC[,speciesID],AUC_c.list[[speciesID]])
    lineHeight <- median(na.omit(AUC[,speciesID]))
  }else if(metric=='RMSE'){
    boxDat <- cbind(RMSE[,speciesID],RMSE_c.list[[speciesID]])
    lineHeight <- median(na.omit(RMSE[,speciesID]))
  }else if(metric=='K'){
    boxDat <- cbind(K[,speciesID],K_c.list[[speciesID]])
    lineHeight <- median(na.omit(K[,speciesID]))
  }
  
  boxplot(boxDat,col = clrs, main=paste0(colnames(ydat)[speciesID],' cross validation analysis ',metric),
          xlim=c(1,length(ins)+1),
          xaxt = "n",
          ylab=metric,
          xlab='Number of conditional species given',
          whisklty=1)
  axis(1,at=c(1:14),
       labels=c(0,1,2,3,4,5,6,7,8,9,10,20,30,59))
  abline(h=lineHeight,lty='dashed')
  lines(c(.58,1.42),c(lineHeight,lineHeight),col='white',lty='dashed')
}
manipulate(mf(i,metric),i=slider(1,60),metric=picker('AUC','RMSE','K'))





#======================================
# make all the boxplots for sup
#======================================
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','Eastern blacknose dace','Eastern creek chubsucker','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mountain redbelly dace','Northern hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')
metric<-'AUC'
par(mfrow=c(2,2), lwd=.3, cex.axis = 1,cex.main = 1.3, cex.lab = 1.2,
    oma=c(.2, 1, 1, 1), mar=c(4,2.7,1,1),mgp=c(1.5,.5,0))
for(i in dataFileLengths$speciesID){
  speciesID <- i
  clrs <- c('orange',rep('light blue',59))
  
  boxDat <- cbind(AUC[,speciesID],AUC_c.list[[speciesID]])
  lineHeight <- median(na.omit(AUC[,speciesID]))
  
  boxplot(boxDat,col = clrs, main='',
          xlim=c(1,length(ins)+1),
          xaxt = "n",
          ylab=metric,
          xlab='Number of species conditioned on',
          whisklty=1,pch=16,cex=.3,
          boxwex=.65)
  title(plotNames[speciesID], line = .3)
  axis(1,at=c(1:14),
       labels=c(0,1,2,3,4,5,6,7,8,9,10,20,30,53))
  abline(h=lineHeight,lty='dashed')
}



#=========================================
# make one set of boxplots for main text
#=========================================
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','East. blacknose dace','East. creek chubsuck.','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mount. redbelly dace','North. hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')

metric<-'AUC'
par(mfrow=c(3,2), lwd=.3, cex.axis = .7,cex.main = .7, cex.lab = .7,
    oma=c(.2, 1, .5, 1), mar=c(3,2,1,1),mgp=c(1.35,.5,0))
for(i in c(49,30,21,19,48)){
  speciesID <- i
  clrs <- c('black',rep('white',54))
  
  boxDat <- cbind(AUC[,speciesID],AUC_c.list[[speciesID]])
  lineHeight <- median(na.omit(AUC[,speciesID]))
  
  boxplot(boxDat,col = clrs, main=paste0(plotNames[speciesID]),
          xlim=c(1,length(ins)+1),
          xaxt = "n",
          ylab=metric,
          xlab='Number of conditional species given',
          whisklty=1,pch=16,cex=.3,
          boxwex=.65)
  text(x = 1:14,
       y = par("usr")[3] - .1*(par("usr")[4]-par("usr")[3]),
       labels = c(0:10,20,30,53),
       xpd = NA,
       #Rotate the labels by 35 degrees.
       srt = 35,
       cex = 0.6)
  abline(h=lineHeight,lty='dashed')
  lines(c(.58,1.42),c(lineHeight,lineHeight),col='white',lty='dashed')
  
  #draw ticks
  for(i in 1:14){
    x.i <- rep(1:14)[i]
    y.i <- 
  lines(x=c(x.i,x.i),
        y=c(par("usr")[3], par("usr")[3] - .05*(par("usr")[4]-par("usr")[3])),col='black',
        xpd=NA)
  }
}



#-----------------------------------------------------------------------------------
# boxplots of difference between marginal and most improved conditional predictions
#-----------------------------------------------------------------------------------
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','East. blacknose dace','East. creek chubsuck.','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mount. redbelly dace','North. hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')
metric<-'AUC'
par(mfrow=c(1,1), lwd=.3, cex.axis = .7,cex.main = .7, cex.lab = .7,
    oma=c(2.3, 1, .5, 1), mar=c(3,2,1,1), mgp=c(1.35,.5,0))
largestMedianImprovement <- data.frame(matrix(NA,ncol=nSpecies,nrow=nrow(AUC)))
colMeans <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  percentChange <- na.omit(percentChange)
  percentChange <- percentChange[!is.infinite(rowSums(percentChange)),]
  medPercentChange <- apply(FUN=median,X=percentChange,MARGIN=2)
  i.index <- which.max(medPercentChange)
  largestMedianImprovement[,i] <- 100*(AUC_c.list[[i]][,i.index]-AUC[,i])/AUC[,i]
  colMeans[i] <- max(medPercentChange)
}
colMeans <- data.frame(cbind(mean=as.numeric(colMeans),id=1:nSpecies))
colMeans <- colMeans[order(colMeans$mean),]
clrs <- ifelse(colMeans$mean>.05,'gray','white')

boxplot(largestMedianImprovement[,colMeans$id], outline=FALSE, col=clrs,
        xaxt = "n",
        #xlab='Species ID',
        ylab='Percent AUC difference',
        #main='Percent difference between best conditional and marginal AUCs',
        whisklty=1,
        at=seq(1,nSpecies*1.5,by=1.5),
        boxwex=1.1)

text(x = seq(1,nSpecies*1.5,by=1.5)-2.2,
     y = par("usr")[3] - .03*(par("usr")[4]-par("usr")[3]),
     labels = plotNames[colMeans$id],
     xpd = NA,
     #Rotate the labels by 35 degrees.
     srt = 280,
     cex = 0.38,pos=4)
abline(h=0,lty=2)

#draw ticks
for(i in 1:nSpecies){
  x.i <- seq(1,nSpecies*1.5,by=1.5)[i]
  y.i <- lines(x=c(x.i,x.i),
               y=c(par("usr")[3], par("usr")[3] - .02*(par("usr")[4]-par("usr")[3])),col='black',
               xpd=NA)
}


#==================================================================================================
# see how many species needed to be conditionally given to be within 10% of the maximum improvement
#==================================================================================================
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','East. blacknose dace','East. creek chubsuck.','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mount. redbelly dace','North. hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')
minMax <- rep(NA,nSpecies)
percentThresh <- 0.20
maxIndex <- rep(NA,nSpecies)
percentChanges <- rep(NA,nSpecies)
cMedians <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  maxIndex[i] <- which.max(apply(X=AUC_c.list[[i]],FUN=median,MARGIN=2))
  cMedian <- median(na.omit(AUC_c.list[[i]][,maxIndex[i]]))
  mMedian <- median(na.omit(AUC[,i]))
  percentChange <- (cMedian-mMedian)/mMedian
  percentChanges[i] <- percentChange
  cMedians[i] <- cMedian
  minMax.i <- c()
  for(j in 1:length(ins)){
    percentChange.j <- (median(na.omit(AUC_c.list[[i]][,j]))-median(na.omit(AUC[,i])))/median(na.omit(AUC[,i]))
    if(abs(percentChange.j-percentChange)<(percentThresh*percentChange)){
      minMax.i <- c(minMax.i,j)
    }
  }
  minMax[i] <- min(minMax.i)
}
#combine common and scientific names
sNames <- c("Anguilla rostrata", "Etheostoma zonale", "Cottus caeruleomentum", "Lepomis macrochirus", "Nocomis leptocephalus", "Enneacanthus gloriosus", "Pimephales notatus", "Salvelinus fontinalis", "Ameiurus nebulosus", "Salmo trutta", "Campostoma anomalum", "Esox niger", "Notropis amoenus", "Luxilus cornutus", "Semotilus atromaculatus", "Exoglossum maxillingua", "Rhinichthys atratulus", "Erimyzon oblongus", "Gambusia holbrooki", "Umbra pygmaea", "Semotilus corporalis", "Etheostoma flabellare", "Notemigonus crysoleucas", "Lepomis cyanellus", "Etheostoma blennioides", "Etheostoma nigrum", "Micropterus salmoides", "Lethenteron appendix", "Rhinichthys cataractae", "Noturus insignis", "Chrosomus oreas", "Hypentelium nigricans", "Aphredoderus sayanus", "Cottus girardi", "Lepomis gibbosus", "Oncorhynchus mykiss", "Lepomis auritus", "Esox americanus americanus", "Nocomis micropogon", "Ambloplites rupestris", "Lythrurus ardens", "Notropis rubellus", "Clinostomus funduloides", "Notropis analostanus", "Petromyzon marinus", "Percina peltata", "Micropterus dolomieu", "Cyprinella spiloptera", "Notropis hudsonius", "Notropis procne", "Etheostoma olmstedi", "Thoburnia rhothoeca", "Catostomus commersonii", "Ameiurus natalis")
combinedNames <- rep(NA,nSpecies)
for(i in 1:nSpecies){
  combinedNames[i] <- paste0(colnames(ydat[i])," (",sNames[i],")")
}
#consruct, order and save data.frame
minMaxDF <- data.frame(species=plotNames,minMax=ins[minMax],
                       trueMax=ins[maxIndex],improvement=round(100*percentChanges,2),
                       median_auc=round(cMedians,2),
                       top_cor=corNames[,1:3])
                       
minMaxDF <- minMaxDF[order(minMaxDF$improvement,decreasing=T),]
minMaxDF





#==============================================================
# Scatterplots of marginal and highest conditional median AUC
#==============================================================
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','East. blacknose dace','East. creek chubsuck.','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mount. redbelly dace','North. hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')
metric<-'AUC'
par(mfrow=c(1,1), lwd=.3, cex.axis = .7,cex.main = .7, cex.lab = .7,
    oma=c(2.3, 1, .5, 1), mar=c(3,2,1,1), mgp=c(1.35,.5,0))
largestMedianImprovement <- data.frame(matrix(NA,ncol=nSpecies,nrow=nrow(AUC)))
colMeans <- rep(NA,nSpecies)
bestConAUC <- data.frame(matrix(NA,nrow=100,ncol=nSpecies))
for(i in c(1:nSpecies)){
  percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  percentChange <- na.omit(percentChange)
  percentChange <- percentChange[!is.infinite(rowSums(percentChange)),]
  medPercentChange <- apply(FUN=median,X=percentChange,MARGIN=2)
  i.index <- which.max(medPercentChange)
  largestMedianImprovement[,i] <- 100*(AUC_c.list[[i]][,i.index]-AUC[,i])/AUC[,i]
  colMeans[i] <- max(medPercentChange)
  
  bestConAUC[,i] <- AUC_c.list[[i]][,i.index]
}
colMeans <- data.frame(cbind(mean=as.numeric(colMeans),id=1:nSpecies))
colMeans <- colMeans[order(colMeans$mean),]
clrs <- ifelse(colMeans$mean>.05,'gray','white')

bestConAUC_meds <- as.numeric(apply(X=bestConAUC,FUN=median,MARGIN=2,na.rm=TRUE))
AUC_meds <- apply(X=AUC,FUN=median,MARGIN=2,na.rm=TRUE)
plot_AUC_meds_data <- data.frame(con=bestConAUC_meds,marg=AUC_meds,percentDiff=(bestConAUC_meds-AUC_meds)/AUC_meds,speciesID=1:nSpecies)
plot_AUC_meds_data <- plot_AUC_meds_data[order(plot_AUC_meds_data$marg),]

plot(seq(1,nSpecies*1.5,by=1.5),plot_AUC_meds_data$con,
     #main='Marginal and best conditional median AUC',
     ylab='Median  AUC',xaxt = "n",xlab='',cex=.5,
     ylim=c(min(na.omit(plot_AUC_meds_data$marg)),max(na.omit(plot_AUC_meds_data$con))))
points(seq(1,nSpecies*1.5,by=1.5),plot_AUC_meds_data$marg,pch=16,cex=.5)

for(i in c(1:nSpecies)){
  if(plot_AUC_meds_data$percentDiff[i]>0){
    lines(rep(seq(1,nSpecies*1.5,by=1.5)[i],each=2),rbind(plot_AUC_meds_data$con[i]-.003,plot_AUC_meds_data$marg[i]))
  }else{
    lines(rep(seq(1,nSpecies*1.5,by=1.5)[i],each=2),rbind(plot_AUC_meds_data$con[i]+.003,plot_AUC_meds_data$marg[i]))
  }
}

text(x = seq(1,nSpecies*1.5,by=1.5)-2.2,
     y = par("usr")[3] - .03*(par("usr")[4]-par("usr")[3]),
     labels = plotNames[plot_AUC_meds_data$speciesID],
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 280,
     cex = 0.38,pos=4)
abline(h=0,lty=2)

#draw ticks
for(i in 1:nSpecies){
  x.i <- seq(1,nSpecies*1.5,by=1.5)[i]
  y.i <- lines(x=c(x.i,x.i),
               y=c(par("usr")[3], par("usr")[3] - .02*(par("usr")[4]-par("usr")[3])),col='black',
               xpd=NA)
}
abline(h=.9,lty=2)
abline(h=.8,lty=2)
abline(h=.7,lty=2)



#==========================================================================
# make table of biggest improvements ordered by all three: AUC, RMSE and K
#==========================================================================
best <- data.frame(matrix(NA,nrow=nSpecies,ncol=7))
K_percents <- data.frame(matrix(NA,nrow=nSpecies,ncol=13))
bestK_meds <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  #---------------------------------------------------------------
  # calculate improvements (actual, not percents) conditional in 
  # AUC, RMSE and K for reach set of conditional information
  #---------------------------------------------------------------
  AUC_percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  AUC_percentChange <- na.omit(AUC_percentChange)
  AUC_percentChange <- AUC_percentChange[!is.infinite(rowSums(AUC_percentChange)),]
  AUC_medPercentChange <- apply(FUN=median,X=AUC_percentChange,MARGIN=2)
  
  RMSE_percentChange <- (RMSE_c.list[[i]]-RMSE[,i])/RMSE[,i]
  RMSE_percentChange <- na.omit(RMSE_percentChange)
  RMSE_percentChange <- RMSE_percentChange[!is.infinite(rowSums(RMSE_percentChange)),]
  RMSE_medPercentChange <- apply(FUN=median,X=RMSE_percentChange,MARGIN=2)
  
  K_percentChange <- (K_c.list[[i]]-K[,i])/K[,i]
  K_percentChange <- na.omit(K_percentChange)
  K_percentChange <- K_percentChange[!is.infinite(rowSums(K_percentChange)),]
  K_medPercentChange <- apply(FUN=median,X=K_percentChange,MARGIN=2)

  #find conditional column index with largest median for all three
  medPercents <- data.frame(auc=AUC_medPercentChange,
                            rmse=(-1)*RMSE_medPercentChange,
                            k=K_medPercentChange)
  medPercents_strd <- apply(X=medPercents,FUN=scale,MARGIN=2)
  bestIndex <- which.max(rowMeans(medPercents_strd))
  
  best[i,] <- c(100*AUC_medPercentChange[bestIndex],
                median(na.omit(AUC_c.list[[i]][,bestIndex])),
                100*RMSE_medPercentChange[bestIndex],
                median(na.omit(RMSE_c.list[[i]][,bestIndex])),
                100*K_medPercentChange[bestIndex],
                median(na.omit(K_c.list[[i]][,bestIndex])),
                #median(na.omit(K[,i])),
                bestIndex)
  bestK_meds[i] <- median(na.omit(K_c.list[[i]][,bestIndex]))
}
names(best) <- c('auc_improv','auc','rmse_improv','rmse','k_improv','k','nConSpecies')
best <- data.frame(speciesName=colnames(ydat),best)
best <- best[order(best$auc_improv,decreasing=TRUE),]
best$nConSpecies <- ifelse(best$nConSpecies==11,20,
                           ifelse(best$nConSpecies==12,30,
                                  ifelse(best$nConSpecies==13,nSpecies-1,best$nConSpecies)))
best$rmse_improv <- -best$rmse_improv
bestRounded <- best
bestRounded[,2:7] <- round(best[,2:7],2)
bestRounded

#=====================================================
# check if any species changed categories of K or AUC
#=====================================================
#-----------------------------------------------------
best <- data.frame(matrix(NA,nrow=nSpecies,ncol=10))
bestK_meds <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  #---------------------------------------------------------------
  # calculate improvements (actual, not percents) conditional in 
  # AUC, RMSE and K for reach set of conditional information
  #---------------------------------------------------------------
  AUC_percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  AUC_percentChange <- na.omit(AUC_percentChange)
  AUC_percentChange <- AUC_percentChange[!is.infinite(rowSums(AUC_percentChange)),]
  AUC_medPercentChange <- apply(FUN=median,X=AUC_percentChange,MARGIN=2)
  
  RMSE_percentChange <- (RMSE_c.list[[i]]-RMSE[,i])/RMSE[,i]
  RMSE_percentChange <- na.omit(RMSE_percentChange)
  RMSE_percentChange <- RMSE_percentChange[!is.infinite(rowSums(RMSE_percentChange)),]
  RMSE_medPercentChange <- apply(FUN=median,X=RMSE_percentChange,MARGIN=2)
  
  K_percentChange <- (K_c.list[[i]]-K[,i])/K[,i]
  K_percentChange <- na.omit(K_percentChange)
  K_percentChange <- K_percentChange[!is.infinite(rowSums(K_percentChange)),]
  K_medPercentChange <- apply(FUN=median,X=K_percentChange,MARGIN=2)
  
  #find conditional column index with largest median for all three
  medPercents <- data.frame(auc=AUC_medPercentChange,
                            rmse=(-1)*RMSE_medPercentChange,
                            k=K_medPercentChange)
  medPercents_strd <- apply(X=medPercents,FUN=scale,MARGIN=2)
  bestIndex <- which.max(rowMeans(medPercents_strd))
  
  best[i,] <- c(100*AUC_medPercentChange[bestIndex],
                median(na.omit(AUC_c.list[[i]][,bestIndex])),
                median(na.omit(AUC[,i])),
                100*RMSE_medPercentChange[bestIndex],
                median(na.omit(RMSE_c.list[[i]][,bestIndex])),
                median(na.omit(RMSE[,i])),
                100*K_medPercentChange[bestIndex],
                median(na.omit(K_c.list[[i]][,bestIndex])),
                median(na.omit(K[,i])),
                #median(na.omit(K[,i])),
                bestIndex)
  bestK_meds[i] <- median(na.omit(K_c.list[[i]][,bestIndex]))
}
names(best) <- c('auc_improv','auc','auc_m','rmse_improv','rmse','rmse_m','k_improv','k','k_m','nConSpecies')
best <- data.frame(speciesName=colnames(ydat),best)
best <- best[order(best$auc_improv,decreasing=TRUE),]
best$nConSpecies <- ifelse(best$nConSpecies==11,20,
                           ifelse(best$nConSpecies==12,30,
                                  ifelse(best$nConSpecies==13,nSpecies-1,best$nConSpecies)))
best$rmse_improv <- -best$rmse_improv
bestRounded <- best
bestRounded[,2:11] <- round(best[,2:11],2)

#check for K category changes
best$k_cat_change <- NA
for(i in 1:nrow(best)){
  if(best$k[i]>=.81 & best$k_m[i]<=.8){
    best$k_cat_change[i] <- 1
  }else if(best$k[i]>=.61 & best$k_m[i]<=.6){
    best$k_cat_change[i] <- 1
  }else if(best$k[i]>=.41 & best$k_m[i]<=.4){
    best$k_cat_change[i] <- 1
  }else if(best$k[i]>=.21 & best$k_m[i]<=.2){
    best$k_cat_change[i] <- 1
  }else{
    best$k_cat_change[i] <- 0
  }
}
sum(best$k_cat_change)

#check for AUC category changes
best$auc_cat_change <- NA
for(i in 1:nrow(best)){
  if(best$auc[i]>.9 & best$auc_m[i]<=.9){
    best$auc_cat_change[i] <- 1
  }else if(best$auc[i]>.8 & best$auc_m[i]<=.8){
    best$auc_cat_change[i] <- 1
  }else if(best$auc[i]>.7 & best$auc_m[i]<=.7){
    best$auc_cat_change[i] <- 1
  }else{
    best$auc_cat_change[i] <- 0
  }
}
sum(best$auc_cat_change)
  
  
  

#==============================================================
# species most often in the top of residual correlated species
# for target species which observed improvements in AUC of at 
# least 5%
#==============================================================
head(allIDs_mat)
AUC_improv <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  AUC_percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  AUC_percentChange <- na.omit(AUC_percentChange)
  AUC_percentChange <- AUC_percentChange[!is.infinite(rowSums(AUC_percentChange)),]
  AUC_medPercentChange <- apply(FUN=median,X=AUC_percentChange,MARGIN=2)
  AUC_improv[i] <- max(na.omit(AUC_medPercentChange))
}
bestIndex <- which(AUC_improv>.05);length(bestIndex)
highestCorr_best <- corIndexes[bestIndex,1:5]
namesHighestCorr <- colnames(ydat)[highestCorr_best]
mostHighestCorr <- sort(table(namesHighestCorr))
mostHighestCorr

sort(match(names(mostHighestCorr)[(length(mostHighestCorr)-15):length(mostHighestCorr)],names(sort(colSums(ydat_pa)))))



#==========================================================================
# make table of biggest improvements just in AUC
#==========================================================================
best <- data.frame(matrix(NA,nrow=nSpecies,ncol=7))
K_percents <- data.frame(matrix(NA,nrow=nSpecies,ncol=13))
bestK_meds <- rep(NA,nSpecies)
for(i in c(1:nSpecies)){
  #---------------------------------------------------------------
  # calculate improvements (actual, not percents) conditional in 
  # AUC, RMSE and K for reach set of conditional information
  #---------------------------------------------------------------
  AUC_percentChange <- (AUC_c.list[[i]]-AUC[,i])/AUC[,i]
  AUC_percentChange <- na.omit(AUC_percentChange)
  AUC_percentChange <- AUC_percentChange[!is.infinite(rowSums(AUC_percentChange)),]
  AUC_medPercentChange <- apply(FUN=median,X=AUC_percentChange,MARGIN=2)
  
  RMSE_percentChange <- (RMSE_c.list[[i]]-RMSE[,i])/RMSE[,i]
  RMSE_percentChange <- na.omit(RMSE_percentChange)
  RMSE_percentChange <- RMSE_percentChange[!is.infinite(rowSums(RMSE_percentChange)),]
  RMSE_medPercentChange <- apply(FUN=median,X=RMSE_percentChange,MARGIN=2)
  
  K_percentChange <- (K_c.list[[i]]-K[,i])/K[,i]
  K_percentChange <- na.omit(K_percentChange)
  K_percentChange <- K_percentChange[!is.infinite(rowSums(K_percentChange)),]
  K_medPercentChange <- apply(FUN=median,X=K_percentChange,MARGIN=2)
  
  #find conditional column index with largest median for all three
  medPercents <- data.frame(auc=AUC_medPercentChange,
                            rmse=(-1)*RMSE_medPercentChange,
                            k=K_medPercentChange)
  medPercents_strd <- apply(X=medPercents,FUN=scale,MARGIN=2)
  bestIndex <- which.max(rowMeans(medPercents_strd))
  
  best[i,] <- c(100*AUC_medPercentChange[bestIndex],
                median(na.omit(AUC_c.list[[i]][,bestIndex])),
                100*RMSE_medPercentChange[bestIndex],
                median(na.omit(RMSE_c.list[[i]][,bestIndex])),
                100*K_medPercentChange[bestIndex],
                median(na.omit(K_c.list[[i]][,bestIndex])),
                #median(na.omit(K[,i])),
                bestIndex)
  bestK_meds[i] <- median(na.omit(K_c.list[[i]][,bestIndex]))
}
names(best) <- c('auc_improv','auc','rmse_improv','rmse','k_improv','k','nConSpecies')
best <- data.frame(speciesName=colnames(ydat),best)
#best <- best[order(best$auc_improv,decreasing=TRUE),]
best$nConSpecies <- ifelse(best$nConSpecies==11,20,
                           ifelse(best$nConSpecies==12,30,
                                  ifelse(best$nConSpecies==13,nSpecies-1,best$nConSpecies)))
best$rmse_improv <- -best$rmse_improv
bestRounded <- best
best$auc_m <- apply(X=AUC,FUN=median,MARGIN=2)
bestRounded[,2:7] <- round(best[,2:7],2)
bestRounded


#============================================================================
# make table 2 from supplmentary material with species names and # observers
#============================================================================
plotNames <- c('American eel','Banded darter','Blue ridge sculpin','Bluegill','Bluehead chub','Bluespotted sunfish','Bluntnose minnow','Brook trout','Brown bullhead','Brown trout','Central stoneroller','Chain pickerel','Comely shiner','Common shiner','Creek chub','Cutlips minnow','Eastern blacknose dace','Eastern creek chubsucker','Eastern mosquitofish','Eastern mudminnow','Fallfish','Fantail darter','Golden shiner','Green sunfish','Greenside darter','Johnny darter','Largemouth bass','Least brook lamprey','Longnose dace','Margined madtom','Mountain redbelly dace','Northern hog sucker','Pirate perch','Potomac sculpin','Pumpkinseed','Rainbow trout','Redbreast sunfish','Redfin pickerel','River chub','Rock bass','Rosefin shiner','Rosyface shiner','Rosyside dace','Satinfin shiner','Sea lamprey','Shield darter','Smallmouth bass','Spotfin shiner','Spottail shiner','Swallowtail shiner','Tessellated darter','Torrent sucker','White sucker','Yellow bullhead')
sNames <- c("Anguilla rostrata", "Etheostoma zonale", "Cottus caeruleomentum", "Lepomis macrochirus", "Nocomis leptocephalus", "Enneacanthus gloriosus", "Pimephales notatus", "Salvelinus fontinalis", "Ameiurus nebulosus", "Salmo trutta", "Campostoma anomalum", "Esox niger", "Notropis amoenus", "Luxilus cornutus", "Semotilus atromaculatus", "Exoglossum maxillingua", "Rhinichthys atratulus", "Erimyzon oblongus", "Gambusia holbrooki", "Umbra pygmaea", "Semotilus corporalis", "Etheostoma flabellare", "Notemigonus crysoleucas", "Lepomis cyanellus", "Etheostoma blennioides", "Etheostoma nigrum", "Micropterus salmoides", "Lethenteron appendix", "Rhinichthys cataractae", "Noturus insignis", "Chrosomus oreas", "Hypentelium nigricans", "Aphredoderus sayanus", "Cottus girardi", "Lepomis gibbosus", "Oncorhynchus mykiss", "Lepomis auritus", "Esox americanus americanus", "Nocomis micropogon", "Ambloplites rupestris", "Lythrurus ardens", "Notropis rubellus", "Clinostomus funduloides", "Notropis analostanus", "Petromyzon marinus", "Percina peltata", "Micropterus dolomieu", "Cyprinella spiloptera", "Notropis hudsonius", "Notropis procne", "Etheostoma olmstedi", "Thoburnia rhothoeca", "Catostomus commersonii", "Ameiurus natalis")
namesAndObs <- data.frame(comm_name=plotNames,
                          sci_name=sNames,
                          nObs=colSums(ydat_pa))
namesAndObs

#==================================================================
# analyze marginal predictions
#==================================================================
colnames(AUC) <- colnames(ydat)
colnames(K) <- colnames(ydat)
colnames(RMSE) <- colnames(ydat)
auc_m <- sort(apply(X=AUC,FUN=median,MARGIN=2,na.rm=T))
k_m <- sort(apply(X=K,FUN=median,MARGIN=2,na.rm=T))
rmse_m <- sort(apply(X=RMSE,FUN=median,MARGIN=2,na.rm=T,decreasing=T))

data.frame(names1 = names(auc_m),auc=auc_m,
      names2=names(k_m),k=k_m,
      names3=names(rmse_m),rmse=rmse_m)

#what percent fell into poor-good range for AUC (outstanding $>$ 0.9, excellent $>$ 0.8, acceptable $>$ 0.7, poor $\leq 0.7$\)
auc_cat <- rep(NA,nSpecies)
for(i in 1:nSpecies){
  auc.i <- auc_m[i]
  if(auc.i>=.9){
    auc_cat[i] <-'outstanding'
  }
  if(auc.i<.9 & auc.i>=.8){
    auc_cat[i] <- 'excelent'
  }
  if(auc.i<.8 & auc.i>=.7){
    auc_cat[i] <- 'acceptable'
  }
  if(auc.i<.7){
    auc_cat[i] <- 'poor'
  }
}
table(auc_cat)


#what percent fell into poor-good range for K almost perfect$>$0.81, substantial$>$.61, moderate$>$0.41, fair$>0.21$
k_cat <- rep(NA,nSpecies)
for(i in 1:nSpecies){
  k.i <- k_m[i]
  if(k.i>.8){
    k_cat[i] <-'almost perfect'
  }
  if(k.i<=.8 & k.i>.60){
    k_cat[i] <- 'substantial'
  }
  if(k.i<=.6 & k.i>.4){
    k_cat[i] <- 'moderate'
  }
  if(k.i<=.4 & k.i>.2){
    k_cat[i] <- 'fair'
  }
  if(k.i<.2){
    k_cat[i] <- 'no to slight agreement'
  }
}
table(k_cat)




#================================================================================
# Scatterplot of marginal AUC vs. mean of standardized standard deviations of the
# predictor variables each species was observed for.
#================================================================================
cvIndex <- c(4,5,7:10,12:15,17)
SDs <- data.frame(matrix(NA,nrow=nSpecies,ncol=length(cvIndex)))
for(i in c(1:nSpecies)){
  xdat.i <- orig_xdat[which(ydat_pa[,i]>0),]
  for(j in 1:length(cvIndex)){
    cn <- cvIndex[j]
    SDs[i,j] <-  as.numeric(sd(xdat.i[,cn]))
  }
}
SDs$speciesName <- colnames(ydat)
SDs$id <- 1:nSpecies
names(SDs) <- c(colnames(orig_xdat)[cvIndex],'speciesNames','id')


#mean of ranked IQR just color (save PDF: 5.91 x 3.86)
SDs_z <- SDs
for(i in 1:length(cvIndex)){
  for(j in 1:nSpecies){
    SDs_z[j,i] <- (SDs[j,i]-mean(SDs[,i]))/sd(SDs[,i])
  }
}
SDs_z$mean <- NA
negLengths <- rep(NA,nSpecies)
for(i in 1:nSpecies){
  SDs_z.i <- as.numeric(SDs_z[i,1:length(cvIndex)])
  SDs_z$mean[i] <- mean(SDs_z.i[which(SDs_z.i<0)])
  negLengths[i] <- length(SDs_z.i[which(SDs_z.i<0)])
}
#SDs_z$mean <- apply(X=SDs_z[,1:length(cvIndex)],FUN=mean,MARGIN=1) 
SDs_z$improvement <- best$auc_improv
SDs_z$auc_m <- best$auc_m
SDs_z$auc <- best$auc
cor.test(SDs_z$mean,SDs_z$auc_m)

ggplot(data=SDs_z,aes(mean,auc_m))+
  geom_point(aes(size=improvement,colour=improvement),alpha=.7)+
  geom_point(aes(size=improvement),shape=1)+
  scale_colour_gradient(low="white", high="black")+
  scale_size(range=c(0.01,11)) +
  guides(color=guide_legend(title='AUC \nimprovement'), size = guide_legend(title='AUC \nimprovement'))+
  theme_classic()+
  geom_text_repel(aes(mean,auc_m,label = speciesNames),
                  size = 1.3,max.overlaps = Inf)+
  labs(x='Mean of negative standardized standard devations',
       y='Marginal AUC')



#-----------------------------------------------------------------------
# Scatterplot of highest mean residual correlations vs AUC improvements
#-----------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dirname(getwd()))
load("fullRun_54Species/fullRun.RData")
out$parameters$corMu
cors_diagZero <- out$parameters$corMu
diag(cors_diagZero) <- NA
corMeans <- colMeans(abs(cors_diagZero),na.rm=T)
best$corMeans <- corMeans

plot(corMeans,best$auc_improv, pch=1,
     #col=ifelse(best$auc_m>.85,'grey','black'),pch=16,
     ylab='Percent improvement in AUC',
     xlab='Mean of absolute value of residual correlations',
     cex=match(best$auc_m,sort(best$auc_m))/20)
cor.test(corMeans,best$auc_improv)

ggplot(data=best,aes(corMeans,auc_improv))+
  geom_point(aes(size=auc_m,colour=auc_m),alpha=.7)+
  geom_point(aes(size=auc_m),shape=1)+
  scale_colour_gradient(low="white", high="black")+
  scale_size(range=c(0.01,5)) +
  guides(color=guide_legend(title='marginal \nAUC'), size = guide_legend(title='marginal \nAUC'))+
  theme_classic()+
  geom_text_repel(aes(corMeans,auc_improv,label = speciesName),
                  size = 1.3,max.overlaps = Inf)+
  labs(x='Mean of absoulute value of the residual correlations',
       y='Conditional AUC improvement')
cor.test(corMeans,best$auc_improv)







