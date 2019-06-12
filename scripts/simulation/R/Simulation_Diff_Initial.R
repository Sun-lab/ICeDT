
library(ICeDT)

dimX <- function(v){ if(is.null(dim(v))){r=length(v)} else{r=dim(v)}; r }

#-------------------------------------------------------------------#
# Call Simulation Function                                          #
#-------------------------------------------------------------------#

pctAb     = 20
nCT       = 5
nrep      = 5
subjPur   = 0.45
nStarts   = 1000

# simulate data for 135 individuals
set.seed(666)
data(mean_var_relation)
simData = simFunc(nS=135, nG=250, nrep=nrep, nCT=nCT, pctAb=pctAb,
                  meanVar_Rel=mean_var_relation)
names(simData)
lapply(simData, dimX)
simData$tumor_purity[1:5]
simData$PGene_Exp[1:2,1:7]

# choose the subjects with given tumor purity
chooseSubj = which.min(abs(simData$tumor_purity - subjPur))
simData$tumor_purity[chooseSubj]
simData$rho[,chooseSubj, drop=F]

simData$singleVec    = simData$MGene_Exp[,chooseSubj, drop=F]
simData$single_GnExp = matrix(simData$singleVec, nrow=nrow(simData$singleVec), 
                            ncol=nStarts)
simData$singleRho    = simData$rho[,chooseSubj, drop=F]


# simulate initia values from non-informative priors
startPoints = matrix(NA, nrow=nCT, ncol=nStarts)

for(i in 1:nStarts){
  startPoints[,i] = rdirichlet(n=1, alpha=rep(1,nCT))
}

simData$startPoints  = startPoints

#-------------------------------------------------------------------#
# estimate signature matrix
#-------------------------------------------------------------------#

X = simData$PGene_Exp
dim(X)
X[1:2,1:7]
X = X[,which(colnames(X) != "tumor")]

refE = refEstimate(X)
lapply(refE, dimX)

refE$refMat[1:2,]
refE$refVar[1:2,]
refE$ctMu[1:2,]
refE$ctVar[1:2,]

cor(c(refE$refVar), c(refE$ctVar))

#-------------------------------------------------------------------#
# estimation without giving initial values of rho
#-------------------------------------------------------------------#

Y = simData$single_GnExp
Z = refE$refMat
tumorPurity  = NULL
refVar       = NULL
rhoInit      = NULL
maxIter_prop = 100
maxIter_PP   = 100
rhoConverge  = 1e-3

startTime = proc.time()
Testing_Orig = ICeDT(Y = Y, Z = Z, tumorPurity = tumorPurity, 
                     refVar = refVar, rhoInit = rhoInit, 
                     maxIter_prop = maxIter_prop, 
                     maxIter_PP = maxIter_PP, 
                     rhoConverge = rhoConverge)

endTime = proc.time()
endTime - startTime

lapply(Testing_Orig, dimX)

# true value of rho's
simData$singleRho
# estimated values of rho's
summary(t(Testing_Orig$rho))

#-------------------------------------------------------------------#
# estimatino with initial values of rho, no purity
#-------------------------------------------------------------------#

rhoInit = simData$startPoints[-c(1),]
rownames(rhoInit) = colnames(Z)
dim(rhoInit)
rhoInit[,1:5]
summary(colSums(rhoInit))

startTime = proc.time()
Testing_w_rho = ICeDT(Y = Y, Z = Z, tumorPurity = tumorPurity, 
                     refVar = refVar, rhoInit = rhoInit, 
                     maxIter_prop = maxIter_prop, 
                     maxIter_PP = maxIter_PP, 
                     rhoConverge = rhoConverge)
endTime = proc.time()
endTime - startTime

lapply(Testing_w_rho, dimX)
# true value of rho's
simData$singleRho
# estimated values of rho's
summary(t(Testing_w_rho$rho))

#-------------------------------------------------------------------#
# estimatino with initial values of rho, with purity
#-------------------------------------------------------------------#

tumorPurity  = rep(simData$singleRho[1,1], ncol(Y))
rhoInit2 = t(t(rhoInit)/colSums(rhoInit) * (1 - tumorPurity))
dim(rhoInit2)
rhoInit2[,1:5]
summary(colSums(rhoInit2))

startTime = proc.time()
Testing_wr_pu = ICeDT(Y = Y, Z = Z, tumorPurity = tumorPurity, 
                     refVar = refVar, rhoInit = rhoInit2, 
                     maxIter_prop = maxIter_prop, 
                     maxIter_PP = maxIter_PP, 
                     rhoConverge = rhoConverge)
endTime = proc.time()
endTime - startTime

lapply(Testing_wr_pu, dimX)
# true value of rho's
simData$singleRho
# estimated values of rho's
summary(t(Testing_wr_pu$rho))

#-------------------------------------------------------------------#
# summarize the results
#-------------------------------------------------------------------#

png("../figures/simu_diff_initials.png", 
    width=6, height=6, units="in", res=400)
par(mfrow=c(2,2), mar=c(5,4,2,1), bty="n", cex=0.8)

rhoInit2t = t(rbind(tumorPurity,rhoInit2))
summary(colSums(rhoInit2))
rhoInit1t = t(rbind(1-colSums(rhoInit), rhoInit))
colnames(rhoInit1t)[1] = colnames(rhoInit2t)[1] = "purity"

boxplot(rhoInit2t, ylab="initial proportion", main="given tumor purity")
boxplot(rhoInit1t, ylab="initial proportion", main="without tumor purity")

plot(rep(simData$singleRho,each=1000),  c(t(Testing_wr_pu$rho)), 
     xlab="true proportion", ylab="estimated proportion", 
     main="given tumor purity")
abline(0,1)

plot(rep(simData$singleRho,each=1000),  c(t(Testing_w_rho$rho)), 
     xlab="true proportion", ylab="estimated proportion", 
     main="without tumor purity")
abline(0,1)
dev.off()

q(save="no")

