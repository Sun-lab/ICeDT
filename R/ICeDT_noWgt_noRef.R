
meanFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = mean)
  return(out)
}

varFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = var)
  return(out)
}

#-------------------------------------------------------------------#
# UTILITY FUNCTION 1: lmInit                                        #
#-------------------------------------------------------------------#
# INPUT:
#-------------------------------------------------------------------#
# Variable      | Description                                       #
#---------------|---------------------------------------------------#
#      yr       | Vector of length G+1 which contains the gene exp. #
#               | for a single subject across all genes and tumor   #
#               | purity in position 1.                             #
#---------------|---------------------------------------------------#
#      Zm       | Matrix of size (GxQ) which contains the reference #
#               | gene expression profiles for all cell types and   #
#               | genes.                                            #
#---------------|---------------------------------------------------#
#      Zt       | Vector of length G which contains the reference   #
#               | gene expression profile for the tumor tissues.    #
#-------------------------------------------------------------------#

# OUTPUT:
#-------------------------------------------------------------------#
# LABEL         | DESCRIPTION                                       #
#-------------------------------------------------------------------#
#     rho       | Estimate of subject's cell type abundance profile #
#-------------------------------------------------------------------#

lmInit<-function(yr, Zm, Zt){
  eta = yr[1]
  y   = yr[-c(1)]
  
  if(is.null(Zt)){
    offVal = rep(0, length(y))
  } else {
    offVal = Zt*eta
  }
  
  initMod = lm(y~Zm, offset=offVal)
  
  rho  = coef(initMod)[-1]
  wneg = which(rho<0.01)
  rho[wneg] = 0.01
  rho = (rho/sum(rho))*(1-eta)
  
  return(rho)
}

#-------------------------------------------------------------------#
# UTILITY FUNCTION 2: initial estimate of variance  for             #
# consistent or abberant genes                                      #
#-------------------------------------------------------------------#

sigmaInit <- function(x, Z, nG){
  # Qval is the length of rho in the folloiwng. 
  # rho was esimated by lmInit that ignores tumor cell type
  # hence its lenght is ncol(Z) - 1
  Qval   = ncol(Z) - 1
  Y      = x[c(1:nG)]
  logY   = log(Y)
  rho    = x[c((nG+1):(nG+Qval))]
  rho_i0 = x[(nG+Qval+1)]
  
  eta_ij = Z %*% c(0, rho)
  
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



#-------------------------------------------------------------------#
# Update Functions (HOMOSCEDASTIC)                                  # 
#                                                                   #
# DESCRIPTION:                                                      #
#   Contains all functions necessary for the update of the model    #
#   parameters, in the following order:                             #
#       (1) Subject Specific Proportions                            #
#       (2) Pooled/Scaled Aberrant Marker Profile                   #
#       (3) Cancer Profile                                          #
#       (4) EM Weights (Posterior Means)                            #
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
#                SECTION 1 - Update Proportions                     #
#-------------------------------------------------------------------#
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
# Z is still composed of a vector for "tumor" and everything else.
# These functions are used to compute gradients under the case where 
# the tumor purity is known (fix) vs unknown (nofix). 

gradFunc_noFix <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho = c(0, x)
  
  eta_ij = Z%*%rho
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  out = t(Z[,-c(1)]) %*% matrix(c1,ncol=1)
  
  return(out)
}

gradFunc_Fix <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho = c(x, 1 - rho_i0 -sum(x))
  rho = c(0, rho)
  
  eta_ij = Z%*%rho
  
  mu_ijC = log(eta_ij) - sigma2C/2
  d_ijC  = logY - mu_ijC
  
  mu_ijA = log(eta_ij) - sigma2A/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijC*EM_wgt/(sigma2C*eta_ij)) + c(d_ijA*(1-EM_wgt)/(sigma2A*eta_ij))
  
  # the derivative was derived using all cell type composition of non-tumor 
  # cell types. Here since the proportion of last cell type is 
  # rho_K = constant - sum(x), so its derivative is the original derivative
  # multiples d [rho_K] / d [rho_k], for k = 1 to K-1, which is - 1 
  # so it is equivalent to replace Z by Z_star
  Qval   = ncol(Z)-1
  Z_star = Z[,-c(1)] %*% rbind(diag(rep(1,(Qval-1))), rep(-1,(Qval-1)))
  
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

#----------------- Likelihood (No Fixed) ----------------------#
# hin	and hin_jacob are needed for Augmented Lagrangian Minimization 
# Algorithm for optimization
# hin: a vector function specifying inequality constraints such that 
# hin[j] > 0 for all j
# hin.jac: Jacobian of hin

hin_func_noFix <- function(x, rho_i0, ...){
  return(c(x-0.005, 1-sum(x)-0.005))
}

hin_jacob_noFix <- function(x, rho_i0, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_noFix <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho  = c(0, x)
  
  eta_ij = Z %*% rho
  # if(any(eta_ij < 0)){
  #   print(head(Z))
  #   print(rho)
  #   stop("0 in eta_ij\n")
  # }
  mu_ijC = log(eta_ij) - sigma2C/2
  mu_ijA = log(eta_ij) - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#----------------- Likelihood (Fixed) ----------------------#
hin_func_Fix <- function(x, rho_i0, ...){
  return(c(x-0.005, 1-rho_i0-sum(x)-0.005))
}

hin_jacob_Fix <-function(x, rho_i0, ...){
  return(rbind(diag(1,length(x)), rep(-1,length(x))))
}

logLik_Fix <- function(x, logY, rho_i0, Z, sigma2C, sigma2A, EM_wgt){
  rho  = c(0, x, 1-rho_i0-sum(x))
  
  eta_ij = Z %*% rho
  mu_ijC = log(eta_ij) - sigma2C/2
  mu_ijA = log(eta_ij) - sigma2A/2
  
  out = sum(EM_wgt*dnorm(logY, mean = mu_ijC, sd = sqrt(sigma2C), log = TRUE)) + 
    sum((1-EM_wgt)*dnorm(logY, mean = mu_ijA, sd = sqrt(sigma2A), log = TRUE))
  
  return(out)
}

#------------------- Actual Update Functions ------------------------#

updatePropn_Single <- function(x, Z, nCell, nG, maxIter_prop, useRho){
  #----------------------------------------#
  # Extract Info                           #
  #----------------------------------------#
  rho_i0 = x[1]
  urho_1 = x[c(2:(nCell+1))]
  logY   = x[c((nCell+2):(nCell+nG+1))]
  EM_wgt = x[c((nCell+nG+2):(nCell+2*nG+1))]
  
  #----------------------------------------#
  # Initialize Values                      #
  #----------------------------------------#
  
  eta_ij     = drop(Z %*% matrix(c(0,urho_1), ncol=1))
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
    if(useRho){
      auglagOut = auglag(par = urho_0[-c(nCell)], fn = logLik_Fix, 
                         gr = gradFunc_Fix, hin = hin_func_Fix, 
                         hin.jac = hin_jacob_Fix, logY = logY, 
                         rho_i0 = rho_i0, Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim = list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1[c(1:(nCell-1))] = auglagOut$par
      urho_1[nCell] = 1-rho_i0-sum(auglagOut$par)
      
      urho_1 = correctRho(est = urho_1, total = 1-rho_i0)
    } else {
      auglagOut = auglag(par = urho_0, fn = logLik_noFix, 
                         gr = gradFunc_noFix, hin = hin_func_noFix, 
                         hin.jac = hin_jacob_noFix, logY = logY, 
                         rho_i0 = rho_i0, Z=Z, sigma2C = sigma2C_0, 
                         sigma2A = sigma2A_0, EM_wgt = EM_wgt, 
                         control.optim=list(fnscale=-1), 
                         control.outer = list(trace=FALSE))
      
      urho_1 = auglagOut$par
      urho_1 = correctRho_v2(urho_1)
    }
    
    #----------------------------------#
    # Update Computational Values      #
    #----------------------------------#
    eta_ij     = drop(Z %*% matrix(c(0,urho_1), ncol=1))
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
updatePropn_All <- function(logY, rho_init, fixedCT_rho, Z, maxIter_prop, 
                            EM_wgt, useRho){
  
  PropnInfo = rbind(fixedCT_rho, rho_init, logY, EM_wgt)
  
  out = apply(X=PropnInfo, MARGIN = 2, FUN = updatePropn_Single, 
              Z=Z, useRho = useRho, maxIter_prop = maxIter_prop, 
              nCell=(nrow(rho_init)), nG = nrow(logY))
  
  #--- Reshape Output ---#
  rho_Mat   = extractFunction(compList = out, element = "rho_i")
  sig2M_Mat = extractFunction(compList = out, element = "sigma2C_i")
  sig2A_Mat = extractFunction(compList = out, element = "sigma2A_i")
  
  return(list(rho_Curr = rho_Mat, sig2M_Curr = sig2M_Mat, 
              sig2A_Curr = sig2A_Mat))
}

#-------------------------------------------------------------------#
#                SECTION 3 - EM Weights                             #
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
    
    eta_ij = Z %*% matrix(c(0, rho_init_i), ncol=1)
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

PropPlus_Update<- function(Y, rho_0, fixedCT_rho, useRho, 
                           sigma2C_0, sigma2A_0, Z, propC_0,
                           maxIter_PP, maxIter_prop, nG,
                           rhoConverge){
  logY = log(Y)
  
  rho_t1     = rho_0
  sigma2C_t1 = sigma2C_0
  sigma2A_t1 = sigma2A_0
  propC_t1   = propC_0
  
  for(j in 1:maxIter_PP){
    #----  Reset the Param   ----#
    rho_t0     = rho_t1
    sigma2C_t0 = sigma2C_t1
    sigma2A_t0 = sigma2A_t1
    propC_t0   = propC_t1
    
    #---- Update EM Weights  ----#
    EM_wgt = updateWgts(logY = logY, rho_init = rho_t0, 
                        sigma2C = sigma2C_t0, sigma2A = sigma2A_t0, 
                        Z = Z, propC = propC_t0)
    
    #---- Update Proportions ----#
    PropP_t1  = updatePropn_All(logY = logY, rho_init = rho_t0, 
                                fixedCT_rho = fixedCT_rho,
                                Z = Z, maxIter_prop = maxIter_prop, 
                                EM_wgt = EM_wgt, useRho = useRho)
    
    rho_t1     = PropP_t1$rho_Curr
    sigma2C_t1 = PropP_t1$sig2M_Curr
    sigma2A_t1 = PropP_t1$sig2A_Curr
    
    #---- Update Cons. MarkP ----#
    propC_t1 = colSums(EM_wgt)/nG
    
    #---- Check Convergence  ----#
    rho_diff = abs(rho_t0 - rho_t1)
    max_rho_diff = max(c(rho_diff))
    
    if(max_rho_diff < rhoConverge){ break }
    
    message("Current PropPlus Iter ",j,": Max Diff of ",max(max_rho_diff))
  }
  
  return(list(Rho = rho_t1, sigma2C = sigma2C_t1, sigma2A = sigma2A_t1, 
              propC = propC_t1, Iter = j, max_rho_diff = max_rho_diff))
}

# no weights means no gene-specific weight, 
# no reference means cell type-specific expression will be estimated
ICeDT_noWgt_noRef <- function(Y, X, cellType, fixedCT = NULL, 
                              fixedCT_rho = NULL, useRho = FALSE, 
                              borrow4SD = TRUE, maxIter_prop = 100, 
                              maxIter_PP = 100, rhoConverge = 1e-4){
  
  #-----------------------------------------------------#
  # Check Input                                         #
  #-----------------------------------------------------#
  
  if(nrow(Y)!=nrow(X)){
    stop("Expression Inputs do not have the same number of genes!")
  }

  if(any(is.na(Y))|any(is.na(X))){
     stop("Y and X must not contain any NA entries. Please remove
           any rows with NA and try again.")
  }
  
  if(any(Y<0) | any(X<0)){
    stop("As normalized expression values, Y and X should be
          non-negative. Please correct entries less than 0
          and try again.")
  }
  
  if(any(Y<1e-4) | any(X<1e-4)){
    message("X or Y contain expression values smaller than 1e-4. Adding
             small correction (1e-5) to ensure log transformation
             is viable.")
    Y = Y+1e-5
    X = X+1e-5
  }
  
  if(is.null(fixedCT)){
    message("No fixed cell type set, All proportions will be estimated.")
    fixedCT_rho = rep(0,length(ncol(Y)))
    if(!is.null(fixedCT_rho)){
      message("since fixeCT is null, values of fixedCT_rho is discarded.")
    }
  } else {
    if(!(fixedCT %in% unique(cellType))){
      stop("The value of fixedCT is not present in cellType vector!")
    }
    if(ncol(Y) != length(fixedCT_rho)){
      stop("ncol of Y does not equal to length of fixedCT_rho.")
    } else if(!all(colnames(Y) == names(fixedCT_rho))){
      stop("colnames of Y does not match with names of fixedCT_rho.")
    }
  }
  
  if(ncol(X) != length(cellType)){
    stop("Number of cell type labels does not match the number of 
         samples in pure references!")
  }
  
  ctTable = table(cellType)
  
  if(min(ctTable)<2){
    stop("At least two pure samples are needed for each cell type.")
  }
  
  #-----------------------------------------------------#
  # EXTRACTING SAMPLE INFO                              #
  #-----------------------------------------------------#

  nnct = tapply(cellType, cellType, length)
  ntot = sum(nnct)
  
  sortCT = sort(unique(cellType))
  
  nG  = nrow(X) # number of genes
  nP  = ncol(X) # number of purified samples
  nS  = ncol(Y) # number of bulk (mixed cell type) samples
  nCT = length(sortCT) # number of cell types
  
  #-----------------------------------------------------#
  # Pure Sample Estimation                              #
  #-----------------------------------------------------#
  
  logX = log(X) 
  
  CT_MU  = t(apply(X = logX, MARGIN = 1, FUN = meanFun, 
                   group = cellType))
  
  CT_var = t(apply(X = logX, MARGIN = 1, FUN = varFun,
                   group = cellType))
  
  #------------------ Quick Check ----------------------#
  
  if(any(colnames(CT_MU) != sortCT)){
    message("Problems with tapply order!")
  }
  
  if(any(colnames(CT_var) != sortCT)){
    message("Problems with tapply order!")
  }
  
  if(any(names(nnct) != sortCT)){
    message("Problems with tapply order!")
  }
  
  #-- Borrow information across cell types for SD estimation --#
  
  if(borrow4SD){
    X1     = logX
    ntot_p = length(cellType) - length(which(cellType==fixedCT))
    
    for(i in 1:nCT){
      if(sortCT[i] == fixedCT){ next }
      wwi = which(cellType == sortCT[i])
      X1[,wwi] = logX[,wwi] - CT_MU[,i]
    }
    
    varXPop = apply(X1[,-which(cellType==fixedCT)], 1, var)
    
    for(i in 1:nCT){
      if(sortCT[i] == fixedCT){ next }
      wi = (nnct[i])/(ntot_p)
      CT_var[,i] = wi*CT_var[,i] + (1-wi)*varXPop
    }
  }
  
  #----------------- Correct Order ---------------------#
  # Rearrange order so that fixed cell type (tumor) is the first one
  fIdx   = which(colnames(CT_MU)==fixedCT)
  CT_MU  = cbind(CT_MU[,fIdx],  CT_MU[,-c(fIdx)])
  CT_var = cbind(CT_var[,fIdx], CT_var[,-c(fIdx)])
  colnames(CT_MU)[1]  = fixedCT
  colnames(CT_var)[1] = fixedCT
  
  #-----------------------------------------------------#
  # INITIALIZATION                                      #
  #-----------------------------------------------------#
  # Z are cell type-specific gene expression in original 
  # scale, it was denoted by gamma in the manuscript
  # Edit to ensure tumor contribution is 0
  Z     = exp(CT_MU + CT_var/2)
  Z[,1] = rep(0, nG)
  
  #----------- Mixture Sample Proportions --------------#
  # lmInit function estimate cell type compositon after 
  # excluding tumor cell type, and thus length(rho_1) = ncol(Z_) -1
  rho_1 = apply(X = rbind(fixedCT_rho, Y), MARGIN=2, FUN=lmInit,
                Zm=Z[,-c(1)], Zt=Z[,1])

  #----------- Aberrant Profile Initial ----------------#
  sigma_1 = apply(X = rbind(Y, rho_1, fixedCT_rho), MARGIN=2, 
                  FUN = sigmaInit, Z = Z, nG = nG)
  
  sigma2C_1 = sigma_1[1,]
  sigma2A_1 = sigma_1[2,]
  
  #---- percent of consistent genes per sample -------#
  propC_1   = rep(0.5, nS)
  
  #-----------------------------------------------------#
  # IMPLEMENT ALGORITHMIC FIT                           #
  #-----------------------------------------------------#
  propC_0   = propC_1
  rho_0     = rho_1
  sigma2C_0 = sigma2C_1
  sigma2A_0 = sigma2A_1
  
  #------ Update Proportions + -----#
  PropPlus_Out = suppressWarnings(
    PropPlus_Update(Y = Y, rho_0 = rho_0, fixedCT_rho = fixedCT_rho, 
                    useRho = useRho, sigma2C_0 = sigma2C_0, 
                    sigma2A_0 = sigma2A_0, Z = Z, propC_0 = propC_0, 
                    maxIter_PP = maxIter_PP, maxIter_prop = maxIter_prop, 
                    nG = nG, rhoConverge = rhoConverge))
  
  rho_1     = PropPlus_Out$Rho
  sigma2C_1 = PropPlus_Out$sigma2C
  sigma2A_1 = PropPlus_Out$sigma2A
  propC_1   = PropPlus_Out$propC
  
  EM_wgt =  updateWgts(logY = log(Y), rho_init = rho_1, 
                       sigma2C = sigma2C_1, sigma2A = sigma2A_1, 
                       Z = Z, propC = propC_1)
  
  #-----------------------------------------------------#
  # OUTPUT                                              #
  #-----------------------------------------------------#
  
  if(useRho){
    fidx = which(sortCT==fixedCT)
    
    rho_final  = cbind(fixedCT_rho,t(rho_1))
    colnames(rho_final) = c(fixedCT,sortCT[-fidx])
    rownames(rho_final) = colnames(Y)
  } else {
    fidx = which(sortCT==fixedCT)
    
    fixedCT_rho_est = 1-colSums(rho_1)
    
    rho_final  = cbind(fixedCT_rho_est,t(rho_1))
    colnames(rho_final) = c(fixedCT,sortCT[-fidx])
    rownames(rho_final) = colnames(Y)
  }
  
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

