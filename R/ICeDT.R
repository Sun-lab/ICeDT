
meanFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = mean)
  return(out)
}

varFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = var)
  return(out)
}

#-------------------------------------------------------------------#
# initial estimate of cell type composition
# yr = c(tumor purity, expression of nG genes)
# Z is the expression signature matrix
#-------------------------------------------------------------------#

lmInit<-function(yr, Z){
  rho_t = yr[1]
  y     = yr[-c(1)]
  
  initMod = lm(y~Z)
  
  rho  = coef(initMod)[-1]
  wneg = which(rho<0.01)
  rho[wneg] = 0.01
  rho = (rho/sum(rho))*(1-rho_t)
  
  return(rho)
}

#-------------------------------------------------------------------#
#  initial estimate of variance for consistent or abberant genes                                      #
#-------------------------------------------------------------------#

sigmaInit <- function(x, Z, nG){
  nCT    = ncol(Z)
  Y      = x[c(1:nG)]
  logY   = log(Y)
  rho    = x[c((nG+1):(nG+nCT))]
  eta_ij = Z %*% rho
  
  resid = Y - eta_ij
  Q3val = quantile(abs(resid), probs = c(0.75))
  g2use = which(abs(resid) > Q3val)
  
  init_val = c(0,0)
  
  # g2use is the initial set of abberant genes. It represents the first 
  # look at aberrant genes by selecting them from the the upper 
  # quartile of residuals from initial fit (worst fit 25%). 
  # The remaining 75% are considered "consistent marker genes".
  
  log_eta_ij = log(eta_ij)
  init_val[1] = sigma2_Update(logY = logY[-c(g2use)], 
                              log_eta_ij = log_eta_ij[-c(g2use)],
                              EM_wgt = rep(1,(nG-length(g2use))), 
                              AB_Up = FALSE)
  
  init_val[2] = sigma2_Update(logY = logY[c(g2use)], 
                              log_eta_ij = log_eta_ij[c(g2use)], 
                              EM_wgt = rep(0,length(g2use)), 
                              AB_Up = TRUE)
  return(init_val)
}



#------------------------ Support Functions ------------------------#

correctRho <- function(est, total){
  rho  = est
  wneg = which(rho<0.005)
  rho[wneg] = 0.005
  
  rho = total*rho/sum(rho)
  
  return(rho)
}

correctRho_v2 <- function(est){
  rho = est
  wneg = which(rho<0.005)
  rho[wneg] = 0.005
  
  return(rho)
}

extractFunction <- function(compList, element){
  out = mapply(compList, FUN = function(x){ get(element, x) })
  return(out)
}

#----------------------- Gradient Functions ------------------------#

gradFunc_noPurity <- function(x, logY, Z, sigma2C, sigma2A, EM_wgt){
  eta_ij = Z %*% x
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  out = t(Z) %*% matrix(c1,ncol=1)
  
  return(out)
}

gradFunc_givenPurity <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho = c(x, 1 - rho_i0 -sum(x))
  eta_ij = Z%*%rho
  
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  # the derivative was derived using all cell type composition of non-tumor 
  # cell types. Since the proportion of last cell type is 
  # rho_K = constant - sum(x), its derivative is the original derivative
  # multiples d [rho_K] / d [rho_k], for k = 1 to K-1, which is - 1 
  # so it is equivalent to replace Z by Z_star
  nCT    = ncol(Z)
  Z_star = Z %*% rbind(diag(rep(1,(nCT-1))), rep(-1,(nCT-1)))
  
  out = t(Z_star) %*% matrix(c1,ncol=1)
  
  return(out)
}

#-------------------------- Sigma Updates --------------------------#
sigma2_Update <- function(logY, log_eta_ij, EM_wgt, AB_Up=FALSE){
  
  varTheta_ij = logY - log_eta_ij
  
  if(AB_Up==FALSE){
    Sigma2_Up = 2*(sqrt((sum(EM_wgt*(varTheta_ij^2))/sum(EM_wgt))+1)-1)
  } else {
    Sigma2_Up = 2*(sqrt((sum((1-EM_wgt)*(varTheta_ij^2))/sum(1-EM_wgt))+1)-1)
  }
  
  return(Sigma2_Up)
}

#----------------- Likelihood (no purity) ----------------------#
# hin	and hin_jacob are needed for Augmented Lagrangian Minimization 
# Algorithm for optimization
# hin: a vector function specifying inequality constraints such that 
# hin[j] > 0 for all j
# hin.jac: Jacobian of hin

hin_func_noPurity <- function(x, ...){
  return(c(x-0.005, 1-sum(x)-0.005))
}

hin_jacob_noPurity <- function(x, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_noPurity <- function(x, logY, Z, sigma2C, sigma2A, EM_wgt){
  log_eta_ij = log(Z %*% x)
  mu_ijC = log_eta_ij - sigma2C/2
  mu_ijA = log_eta_ij - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#----------------- Likelihood (given purity) ----------------------#
hin_func_givenPurity <- function(x, rho_i0, ...){
  return(c(x-0.005, 1-rho_i0-sum(x)-0.005))
}

hin_jacob_givenPurity <-function(x, rho_i0, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_givenPurity <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho  = c(x, 1-rho_i0-sum(x))
  
  eta_ij = Z %*% rho
  mu_ijC = log(eta_ij) - sigma2C/2
  mu_ijA = log(eta_ij) - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#------------------- Actual Update Functions ------------------------#

updatePropn_Single <- function(x, Z, nCT, nG, maxIter_prop, givenPurity){
  #----------------------------------------#
  # Extract Info                           #
  #----------------------------------------#
  rho_i0 = x[1]
  urho_1 = x[c(2:(nCT+1))]
  logY   = x[c((nCT+2):(nCT+nG+1))]
  EM_wgt = x[c((nCT+nG+2):(nCT+2*nG+1))]
  
  #----------------------------------------#
  # Initialize Values                      #
  #----------------------------------------#
  eta_ij     = drop(Z %*% matrix(urho_1, ncol=1))
  log_eta_ij = log(eta_ij)
  
  sigma2C_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij,
                             EM_wgt = EM_wgt, AB_Up = FALSE)
  sigma2A_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij, 
                            EM_wgt = EM_wgt, AB_Up = TRUE)
  
  mu_ijC = log_eta_ij - sigma2C_1/2
  mu_ijA = log_eta_ij - sigma2A_1/2
  
  logLik_1 = sum(EM_wgt*dnorm(logY,mean=mu_ijC,sd=sqrt(sigma2C_1),log=TRUE))+
    sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(sigma2A_1),log=TRUE))
  
  for(k in 1:maxIter_prop){
    #----------------------------------#
    # Reset Current Values             #
    #----------------------------------#
    urho_0    = urho_1
    sigma2C_0 = sigma2C_1
    sigma2A_0 = sigma2A_1
    logLik_0  = logLik_1
    
    #----------------------------------#
    # Update Rho values                #
    #----------------------------------#
    if(givenPurity){
      auglagOut = auglag(par = urho_0[-c(nCT)], fn = logLik_givenPurity, 
                         gr = gradFunc_givenPurity, hin = hin_func_givenPurity, 
                         hin.jac = hin_jacob_givenPurity, logY = logY, 
                         rho_i0 = rho_i0, Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim = list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1[c(1:(nCT-1))] = auglagOut$par
      urho_1[nCT] = 1-rho_i0-sum(auglagOut$par)
      
      urho_1 = correctRho(est = urho_1, total = 1-rho_i0)
    } else {
      auglagOut = auglag(par = urho_0, fn = logLik_noPurity, 
                         gr = gradFunc_noPurity, hin = hin_func_noPurity, 
                         hin.jac = hin_jacob_noPurity, logY = logY, 
                         Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim=list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1 = auglagOut$par
      urho_1 = correctRho_v2(urho_1)
    }
    
    #----------------------------------#
    # Update Computational Values      #
    #----------------------------------#
    eta_ij     = drop(Z %*% matrix(urho_1, ncol=1))
    log_eta_ij = log(eta_ij)
    sigma2C_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij, 
                               EM_wgt = EM_wgt, AB_Up = FALSE)
    sigma2A_1  = sigma2_Update(logY = logY, log_eta_ij = log_eta_ij,
                               EM_wgt = EM_wgt, AB_Up = TRUE)
    
    mu_ijC    = log_eta_ij - sigma2C_1/2
    mu_ijA    = log_eta_ij - sigma2A_1/2
    logLik_1  = sum(EM_wgt*dnorm(logY,mean=mu_ijC,sd=sqrt(sigma2C_1),log=TRUE))+
      sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(sigma2A_1),log=TRUE))
    
    if(max(abs(urho_1-urho_0)) < 1e-3){ break }
  }
  return(list(rho_i=urho_1, iter=k, sigma2C_i=sigma2C_1, sigma2A_i=sigma2A_1))
}

#------------------- Update Multiple Subjects ----------------------#
updatePropn_All <- function(logY, rho_init, tumorPurity, Z, maxIter_prop, 
                            EM_wgt, givenPurity){
  
  PropnInfo = rbind(tumorPurity, rho_init, logY, EM_wgt)
  
  out = apply(X=PropnInfo, MARGIN = 2, FUN = updatePropn_Single, 
              Z=Z, givenPurity = givenPurity, maxIter_prop = maxIter_prop, 
              nCT=(nrow(rho_init)), nG = nrow(logY))
  
  #--- Reshape Output ---#
  rho_Mat   = extractFunction(compList = out, element = "rho_i")
  sig2C_Mat = extractFunction(compList = out, element = "sigma2C_i")
  sig2A_Mat = extractFunction(compList = out, element = "sigma2A_i")
  
  return(list(rho_Curr = rho_Mat, sig2C_Curr = sig2C_Mat, 
              sig2A_Curr = sig2A_Mat))
}

#-------------------------------------------------------------------#
#                            EM Weights                             #
#-------------------------------------------------------------------#

updateWgts <- function(logY, rho_init, sigma2C, sigma2A, Z, propC){
  
  EM_wgt = matrix(NA, nrow=nrow(logY), ncol=ncol(logY))
  
  for(i in 1:ncol(logY)){
    
    logY_i     = logY[,i]
    rho_init_i = rho_init[,i]
    
    #----------------------------------------------#
    # Cmu_ij for Consistent Marker Gene Probs      #
    # Amu_ij for Aberrant Marker Gene Probs        #
    #----------------------------------------------#
    
    eta_ij = Z %*% matrix(rho_init_i, ncol=1)
    Cmu_ij = log(eta_ij) - sigma2C[i]/2
    Amu_ij = log(eta_ij) - sigma2A[i]/2
    C_lLik = dnorm(logY_i, mean = Cmu_ij, sd = sqrt(sigma2C[i]), log = TRUE)
    A_lLik = dnorm(logY_i, mean = Amu_ij, sd = sqrt(sigma2A[i]), log = TRUE)
    
    #-----------------------------------#
    # Compiling Weights                 #
    #-----------------------------------#
    EM_wgt[,i] = 1/(1+((1-propC[i])/propC[i])*exp(A_lLik-C_lLik))
  }
  
  return(EM_wgt)
}

PropPlus_Update<- function(Y, rho_0, tumorPurity, givenPurity, 
                           sigma2C_0, sigma2A_0, Z, propC_0,
                           maxIter_PP, maxIter_prop, nG,
                           rhoConverge){
  logY = log(Y)
  
  rho_t1     = rho_0
  sigma2C_t1 = sigma2C_0
  sigma2A_t1 = sigma2A_0
  propC_t1   = propC_0
  
  for(j in 1:maxIter_PP){
    rho_t0     = rho_t1
    sigma2C_t0 = sigma2C_t1
    sigma2A_t0 = sigma2A_t1
    propC_t0   = propC_t1
    
    # Update EM Weights
    EM_wgt = updateWgts(logY = logY, rho_init = rho_t0, 
                        sigma2C = sigma2C_t0, sigma2A = sigma2A_t0, 
                        Z = Z, propC = propC_t0)
    
    # Update Proportions
    PropP_t1  = updatePropn_All(logY = logY, rho_init = rho_t0, 
                                tumorPurity = tumorPurity,
                                Z = Z, maxIter_prop = maxIter_prop, 
                                EM_wgt = EM_wgt, givenPurity = givenPurity)
    
    rho_t1     = PropP_t1$rho_Curr
    sigma2C_t1 = PropP_t1$sig2C_Curr
    sigma2A_t1 = PropP_t1$sig2A_Curr
    
    #---- Update Cons. MarkP ----#
    propC_t1 = colSums(EM_wgt)/nG
    
    #---- Check Convergence  ----#
    rho_diff = abs(rho_t0 - rho_t1)
    max_rho_diff = max(c(rho_diff))
    
    if(max_rho_diff < rhoConverge){ break }
    
    message("Current PropPlus Iter ",j,": Max Diff of ",max(max_rho_diff))
  }
  
  return(list(rho = rho_t1, sigma2C = sigma2C_t1, sigma2A = sigma2A_t1, 
              propC = propC_t1, Iter = j, max_rho_diff = max_rho_diff))
}

ICeDT <- function(Y, X, refMat=TRUE, tumorPurity=NULL, borrow4SD=TRUE, 
                  maxIter_prop=100, maxIter_PP=100, rhoConverge=1e-4){
  
  #-----------------------------------------------------#
  # Check Input                                         #
  #-----------------------------------------------------#
  
  cellType = colnames(X)
  
  if(any(cellType == "")){
    stop("colnames(X) are cell type labels and cannot be empty.")
  }
  
  if(nrow(Y)!=nrow(X)){
    stop("Expression Inputs do not have the same number of genes!")
  }

  if(any(is.na(Y))|any(is.na(X))){
     stop("Y and X must not contain any NA entries.")
  }
  
  if(any(Y < 0) | any(X < 0)){
    stop("Y and X should be non-negative.")
  }
  
  if(any(Y < 1e-4)){
    message("Adding 1e-5 to Y to ensure vallid log transformation.")
    Y = Y + 1e-5
  }
  
  if(any(X < 1e-4)){
    message("Adding 1e-5 to X to ensure vallid log transformation.")
    X = X + 1e-5
  }
  
  sortCT = sort(unique(cellType))
  
  nG  = nrow(X) # number of genes
  nP  = ncol(X) # number of purified samples
  nS  = ncol(Y) # number of bulk (mixed cell type) samples
  givenPurity = !(is.null(tumorPurity))

  if(!is.null(tumorPurity)){
    if(ncol(Y) != length(tumorPurity)){
      stop("ncol of Y does not equal to the length of tumorPurity.")
    } else if(!all(colnames(Y) == names(tumorPurity))){
      stop("colnames of Y do not match with the names of tumorPurity.")
    }
  }else{
    tumorPurity = rep(0, ncol(Y))
    names(tumorPurity) = colnames(Y)
  }
  
  if(refMat){
    if(any(duplicated(colnames(X)))){
      stop("Reference matirx X duplicated colnames.")
    }
    Z = X[,sortCT]
  }else{
    ctTable = table(cellType)
    
    if(min(ctTable)<2){
      stop("At least two pure samples are needed for each cell type.")
    }
    
    #-----------------------------------------------------#
    # estimate signature matrix                           #
    #-----------------------------------------------------#
    
    logX   = log(X) 
    
    CT_MU  = t(apply(X = logX, MARGIN = 1, FUN = meanFun, 
                     group = cellType))
    
    CT_var = t(apply(X = logX, MARGIN = 1, FUN = varFun,
                     group = cellType))
    
    if(any(colnames(CT_MU) != sortCT) || any(colnames(CT_var) != sortCT)){
      stop("Problems with tapply order!")
    }
    
    # Borrow information across cell types for SD estimation 
    if(borrow4SD){
      X1     = logX
      ntot_p = length(cellType)
      
      for(i in 1:length(sortCT)){
        wwi = which(cellType == sortCT[i])
        X1[,wwi] = logX[,wwi] - CT_MU[,i]
      }
      
      varXPop = apply(X1, 1, var)
      
      for(i in 1:length(sortCT)){
        wi = (sum(cellType == sortCT[i]))/ntot_p
        CT_var[,i] = wi*CT_var[,i] + (1-wi)*varXPop
      }
    }
    
    #-----------------------------------------------------#
    # Generate cell type-specific expressio matrix        #
    #-----------------------------------------------------#
    Z = exp(CT_MU + CT_var/2)
  }
  
  # lmInit function estimate cell type compositon by linear regression 
  rho_1 = apply(X = rbind(tumorPurity,Y), MARGIN=2, FUN=lmInit, Z=Z)

  # estimate residual variance for consistent/aberrant genes
  sigma_1 = apply(X = rbind(Y, rho_1), MARGIN=2, FUN=sigmaInit, Z=Z, nG=nG)
  
  sigma2C_1 = sigma_1[1,]
  sigma2A_1 = sigma_1[2,]
  
  # percent of consistent genes per sample 
  propC_1   = rep(0.5, nS)
  
  #-----------------------------------------------------#
  # EM algorithm                                        #
  #-----------------------------------------------------#
  propC_0   = propC_1
  rho_0     = rho_1
  sigma2C_0 = sigma2C_1
  sigma2A_0 = sigma2A_1
  
  PropPlus_Out = suppressWarnings(
    PropPlus_Update(Y = Y, rho_0 = rho_0, tumorPurity = tumorPurity, 
                    givenPurity = givenPurity, sigma2C_0 = sigma2C_0, 
                    sigma2A_0 = sigma2A_0, Z = Z, propC_0 = propC_0, 
                    maxIter_PP = maxIter_PP, maxIter_prop = maxIter_prop, 
                    nG = nG, rhoConverge = rhoConverge))
  
  rho_1     = PropPlus_Out$rho
  sigma2C_1 = PropPlus_Out$sigma2C
  sigma2A_1 = PropPlus_Out$sigma2A
  propC_1   = PropPlus_Out$propC
  
  EM_wgt =  updateWgts(logY = log(Y), rho_init = rho_1, 
                       sigma2C = sigma2C_1, sigma2A = sigma2A_1, 
                       Z = Z, propC = propC_1)
  
  #-----------------------------------------------------#
  # OUTPUT                                              #
  #-----------------------------------------------------#
  
  if(givenPurity){
    rho_final  = cbind(tumorPurity, t(rho_1))
  } else {
    tumorPurity_est = 1 - colSums(rho_1)
    rho_final  = cbind(tumorPurity_est, t(rho_1))
  }
  
  colnames(rho_final) = c("tumor", sortCT)
  rownames(rho_final) = colnames(Y)
  
  outList = list(rho            = rho_final,
                 fixedCT        = fixedCT,
                 sigma2C        = sigma2C_1,
                 sigma2A        = sigma2A_1,
                 Z              = Z,
                 CT_var         = CT_var,
                 P_Consistent   = propC_1,
                 PP_Consistent  = EM_wgt)
  
  return(outList)
}

