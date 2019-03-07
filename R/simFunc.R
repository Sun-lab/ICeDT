
#-------------------------------------------------------------------#
# utility functions                                                 #
#-------------------------------------------------------------------#

Identify_Closest <- function(x, pts2compare){
  return(which.min(abs(pts2compare-x)))
}

Id_Close_appVec <- function(x, pts2compare){
  out = apply(X=matrix(x,ncol=1), MARGIN = 1, FUN = Identify_Closest, 
              pts2compare=pts2compare)
  return(out)
}

Id_Close_appMat <- function(X, pts2compare){
  out = apply(X = X, MARGIN = 2, FUN = Id_Close_appVec, 
              pts2compare=pts2compare)
  return(out)
}

#-------------------------------------------------------------------#
# simulation function                                                 #
#-------------------------------------------------------------------#

simFunc <- function(nS, nG, nrep, nCT, pctAb, meanVar_Rel, abMech=1){
  #--------------------------------------------#
  # STEP 1  : Simulate Purified Reference      #
  #--------------------------------------------#
  ct.by.gene = rep(1:(nCT-1), each=floor(nG/(nCT-1)))
  ct.by.gene = c(ct.by.gene, rep((nCT-1), (nG-length(ct.by.gene))))
  ct.by.gene = ct.by.gene+1
  cellType   = rep(1:nCT, each=nrep)
  nSam       = length(cellType)
  
  muX = matrix(nrow=nG, ncol=nCT)
  
  # random number deciding Low, Medium, or High expression levels
  LMH = runif(n = nG, min = 0, max = 1)
  
  for(i in 1:nG){
    if(LMH[i]<0.33){
      muX[i,]   = runif(n = nCT, min = 2, max = 4)
      ct        = ct.by.gene[i]
      muX[i,ct] = runif(1, 3.5, 5)
    } else if(LMH[i]<0.66){
      muX[i,]   = runif(n = nCT, min = 4, max = 6)
      ct        = ct.by.gene[i]
      muX[i,ct] = runif(1, 5.5, 7)
    } else if(LMH[i]<1){
      muX[i,]   = runif(n = nCT, min = 6, max = 8)
      ct        = ct.by.gene[i]
      muX[i,ct] = runif(1, 7.5, 9)
    }
  }
  
  # muX0 saves the mu values before setting the values for tumor to be 0
  # and muX0 is used for senario 3. 
  muX0 = muX
  muX[,1] = rep(-10, nrow(muX))
  
  #-------------------------------------------------#
  # Mean-Variance Relationship                      #
  #-------------------------------------------------#
  
  sqrt_sigmaX_id = matrix(NA, nrow=nrow(muX), ncol=ncol(muX))
  sqrt_sigmaX_id = Id_Close_appMat(X=muX, pts2compare=meanVar_Rel$x)
  
  sqrt_sigmaX = matrix(NA, nrow=nrow(muX), ncol=ncol(muX))
  
  for(i in 1:ncol(sqrt_sigmaX_id)){
    sqrt_sigmaX[,i] = meanVar_Rel$y[sqrt_sigmaX_id[,i]]
  }
  
  sqrt_sigmaX = sqrt_sigmaX + rnorm(nG, 0, 0.085)
  sigmaX = (sqrt_sigmaX)^2
  
  #-------------------------------------------------#
  # Simulated Cell Type-specific Expression         #
  #-------------------------------------------------#  
  X = matrix(NA, nrow=nG, ncol=nSam)
  
  for(i in 1:nG){
    for(j in 1:nSam){
      ct = cellType[j]
      X[i,j] = exp(rnorm(1, mean=muX[i,ct], sd=sigmaX[i,ct]))
    }
  }
  
  #--------------------------------------------#
  # STEP 2  : Simulate Tumor Purity Levels     #
  #--------------------------------------------#
  tp   = rnorm(n = nS, mean = 0.60, sd = 0.15)
  
  if(any(tp<0.17)){
    wM     = which(tp<0.17)
    tp[wM] = 0.17
  }
  
  if(any(tp>0.95)){
    wM     = which(tp>0.95)
    tp[wM] = 0.95
  }
  
  rhos     = matrix(0, ncol=nS, nrow=nCT)
  rhos[1,] = tp
  
  if(nCT==3){
    aval = c(6, 4)
  } else if(nCT==4){
    aval = c(4, 2.5, 1.5)*10/8
  } else {
    aval = c(4, 2.5, 1.5, rep(2/(nCT-4), nCT-4))
  }
  
  temp = rdirichlet(n = nS, alpha = aval)
  rhos[-1,] = t(temp*(1-rhos[1,]))
  
  #--------------------------------------------#
  # STEP 3.A: Simulate Consistent Marker Genes #
  #--------------------------------------------#
  Y = matrix(0, nrow=nG, ncol=nS)
  
  for(k in 1:nS){
    X1   = matrix(NA, nrow=nG, ncol=nCT)
    rho1 = rhos[,k]
    
    for(i in 1:nG){
      for(ct in 1:nCT){
        X1[i,ct] = rnorm(1, mean=muX[i,ct], sd=sigmaX[i,ct])
      }
    }
    
    y = exp(X1) %*% rho1
    
    Y[,k] = y
  }  
  
  #--------------------------------------------#
  # STEP 3.B: Simulate Aberrant Marker Genes   #
  #--------------------------------------------#
  abInd = rbinom(n = nG, size = 1, prob = pctAb/100)
  wAb   = which(abInd==1)
  
  if(abMech == 1){
    for(j in wAb){
      wC = ct.by.gene[j]
      updntm = runif(n = 1, min = 0, max = 1)
      
      if(updntm <= 0.25){
        abInd[j]  = 1
        downVal   = runif(n = nCT, min = 0.25, max = 0.75)
        mnChg     = rep(0, nCT)
        mnChg[wC] = log(downVal[wC])
        tmp_muX   = muX[j,,drop=FALSE] + mnChg
        tmp_sigX  = sigmaX[j,,drop=FALSE]
      } else if(updntm <= 0.50){
        abInd[j]  = 2
        upVal     = runif(n = nCT, min = (4/3), max = 4)
        mnChg     = rep(0, nCT)
        mnChg[wC] = log(upVal[wC])
        tmp_muX   = muX[j,,drop=FALSE] + mnChg
        tmp_sigX  = sigmaX[j,,drop=FALSE]
      } else {
        abInd[j]  = 3
        tmp_muX   = muX0[j,,drop=FALSE]
        tmp_sigX  = sigmaX[j,,drop=FALSE] 
      }
      
      for(k in 1:nS){
        rho1 = rhos[,k]
        for(l in 1:nCT){
          X1[j,l] =  rnorm(1, mean=tmp_muX[l], sd=tmp_sigX[l])
        }
        
        Y[j,k] = exp(X1[j,])%*%rho1
      }
    }
  } else if(abMech==2){
    stop("abMech 2 has been deprecated!")
  } else {
    stop("abMech must equal 1 or 2!")
  }
  
  #--------------------------------------------#
  # STEP 4  : RETURN VALUES                    #
  #--------------------------------------------#
  
  colnames(X) = paste0("CT", cellType)
  colnames(X)[which(colnames(X) == "CT1")] = "tumor"
  
  outList = list()
  outList$MGene_Exp    = Y
  outList$PGene_Exp    = X
  outList$cellType     = cellType
  outList$aberrant     = drop(abInd)
  outList$rho          = rhos
  outList$tumor_purity = drop(rhos[1,])
  
  return(outList)
}
