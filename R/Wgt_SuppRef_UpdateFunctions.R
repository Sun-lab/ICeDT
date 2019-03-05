#####################################################################
# PROJECT 3: Update Functions (HOMOSCEDASTIC)                       #
#####################################################################
# PROGRAM NAME:                                                     #
#   P3_UpdateFunctions.R                                            #
# PROGRAMMER:                                                       #
#   Douglas Roy Wilson, Jr.                                         #
# DATE CREATED:                                                     #
#   5/22/2017                                                       #
# LAST EDIT:                                                        #
#   5/22/2017                                                       #
# VERSION:                                                          #
#   R-3.3.1                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   Contains all functions necessary for the update of the model    #
#   parameters, in the following order:                             #
#       (1) Subject Specific Proportions                            #
#       (2) Pooled/Scaled Aberrant Marker Profile                   #
#       (3) Cancer Profile                                          #
#       (4) EM Weights (Posterior Means)                            #
#####################################################################

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 1 - Update Proportions                     #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
#------------------------ Support Functions ------------------------#
# correctRho <- function(est,total){
#   rho  = est
#   wneg = which(rho<0.005)
#   rho[wneg] = 0.005
#   
#   rho = total*rho/sum(rho)
#   
#   return(rho)
# }
# 
# correctRho_v2 <- function(est){
#   rho = est
#   wneg = which(rho<0.005)
#   rho[wneg] = 0.005
#   
#   return(rho)
# }
# 
# extractFunction<-function(compList,element){
#   out = mapply(compList,FUN = function(x){get(element,x)})
#   return(out)
# }

#------------------------ Hessian Functions ------------------------#
# <Matrix Descriptions for multiple Functions>
# rho    = (rho_{i,1},...,rho_{i,Q})
#
# logY   = (log[y_i1],...,log[y_iG])
#
# Z      = [Gamma_{01},  Gamma_{11},   ....   ....  .... Gamma_{Q1}]
#          [Gamma_{02},  Gamma_{12},   ....   ....  .... Gamma_{Q2}]
#          [    .                       .                          ]
#          [    .                              .                   ]
#          [    .                                    .             ]
#          [Gamma_{0G},  Gamma_{1G},   ....   ....  .... Gamma_{QG}]
#
# Zm     = [Gamma_{11},   ....   ....  .... Gamma_{Q1}]
#          [Gamma_{12},   ....   ....  .... Gamma_{Q2}]
#          [    .          .                          ]
#          [    .                 .                   ]
#          [    .                       .             ]
#          [Gamma_{1G},   ....   ....  .... Gamma_{QG}]
#
# rho_i0 = Tumor Purity
#
# EM_wgt = (1_{T1=1},...,1_{TG=1})

# HS_HessFunc<-function(rho,logY,rho_i0,Z,Zm,eta_ij,Sigma2,Sigma2A,EM_wgt,debugInd=FALSE){
#   rho  = c(rho_i0,rho)
#   
#   if(debugInd==TRUE){
#     eta_ij = NULL
#     eta_ij = Z%*%rho
#   }
#   mu_ij  = log(eta_ij)-Sigma2/2
#   mu_ijA = log(eta_ij)-Sigma2A/2
#   dij    = logY-mu_ij
#   dij_A  = logY-mu_ijA
#   
#   c1  = drop(EM_wgt/((eta_ij^2)*Sigma2)) 
#   c2  = drop((1-EM_wgt)/(Sigma2A*(eta_ij^2)))
#   
#   #---- Aberrant Mixed In ----#
#   out = -t(Zm)%*%(c1*(diag(drop(1+dij)))+c2*diag(drop(1+dij_A)))%*%Zm
#   
#   return(out)
# }

#----------------------- Gradient Functions ------------------------#
# Z is still composed of a vector for "tumor" and everything else.
HS_GradFunc_noFix_Wgt_suppRef <- function(x,logY,rho_i0,Z,Z_star,Sigma2,Sigma2A,
                                          var_wgt,EM_wgt){
  rho = c(0,x)
  
  eta_ij = Z%*%rho
  
  Sigma2_wgt  = var_wgt*Sigma2
  Sigma2A_wgt = var_wgt*Sigma2A
  
  mu_ijM = log(eta_ij)-Sigma2_wgt/2
  d_ijM  = logY-mu_ijM
  
  mu_ijA = log(eta_ij)-Sigma2A_wgt/2
  d_ijA  = logY-mu_ijA
  
  c1 = c(d_ijM*EM_wgt/(Sigma2_wgt*eta_ij))+c(d_ijA*(1-EM_wgt)/(Sigma2A_wgt*eta_ij))
  
  out = t(Z[,-c(1)])%*%matrix(c1,ncol=1)
  
  return(out)
}

HS_GradFunc_Fix_Wgt_suppRef <- function(x,logY,rho_i0,Z,Z_star,Sigma2,Sigma2A,
                                        var_wgt,EM_wgt){
  rho = c(x,1-rho_i0-sum(x))
  rho = c(0,rho)
  
  eta_ij = Z%*%rho
  
  Sigma2_wgt  = var_wgt*Sigma2
  Sigma2A_wgt = var_wgt*Sigma2A
  
  mu_ijM = log(eta_ij)-Sigma2_wgt/2
  d_ijM  = logY - mu_ijM
  
  mu_ijA = log(eta_ij)-Sigma2A_wgt/2
  d_ijA  = logY - mu_ijA
  
  c1 = c(d_ijM*EM_wgt/(Sigma2_wgt*eta_ij))+c(d_ijA*(1-EM_wgt)/(Sigma2A_wgt*eta_ij))
    
  out = t(Z_star)%*%matrix(c1,ncol=1)
    
  return(out)
}

#-------------------------- Sigma Updates --------------------------#
Het_Sigma2_Func_Wgt_suppRef<-function(x,logY,eta_ij,var_wgt,EM_wgt,AB_Up=FALSE){
  Sigma2_wgt  = var_wgt*x
  
  mu_ij = log(eta_ij)-Sigma2_wgt/2
  
  if(AB_Up==FALSE){
    out = sum(EM_wgt*dnorm(x = logY,mean = mu_ij,sd = sqrt(Sigma2_wgt),log = TRUE))
  } else {
    out = sum((1-EM_wgt)*dnorm(x = logY,mean = mu_ij,sd = sqrt(Sigma2_wgt),log = TRUE))
  }
  return(out)
}

HS_Sigma2_Update_Wgt_suppRef<-function(logY,eta_ij,EM_wgt,var_wgt = var_wgt,AB_Up=FALSE){
  
  optimizeOut = optimize(f = Het_Sigma2_Func_Wgt_suppRef,logY = logY,
                         eta_ij=eta_ij,var_wgt=var_wgt,AB_Up=AB_Up,
                         EM_wgt=EM_wgt,lower = 0,upper=10,maximum = TRUE)
  
  Sigma2_Up = optimizeOut$maximum
  
  return(Sigma2_Up)
}

#----------------- Likelihood (No Fixed) ----------------------#
hin_func_noFix_Wgt_suppRef<-function(x,rho_i0,...){
  return(c(x-0.005,1-sum(x)-0.005))
}

hin_jacob_noFix_Wgt_suppRef <-function(x,rho_i0,...){
  return(rbind(diag(1,length(x)),rep(-1,length(x))))
}

HS_ComputeLik_noFix_Wgt_suppRef<-function(x,logY,rho_i0,Z,Z_star,Sigma2,Sigma2A,
                                          var_wgt,EM_wgt){
  rho  = c(0,x)
  
  eta_ij = Z%*%rho
  
  Sigma2_wgt  = var_wgt*Sigma2
  Sigma2A_wgt = var_wgt*Sigma2A
  
  mu_ijM = log(eta_ij)-Sigma2_wgt/2
  mu_ijA = log(eta_ij)-Sigma2A_wgt/2
  
  out = sum(EM_wgt*dnorm(x = logY,mean = mu_ijM,sd = sqrt(Sigma2_wgt),log = TRUE))+
    sum((1-EM_wgt)*dnorm(x = logY,mean = mu_ijA,sd = sqrt(Sigma2A_wgt),log = TRUE))
  return(out)
}

#----------------- Likelihood (No Fixed) ----------------------#
hin_func_Fix_Wgt_suppRef<-function(x,rho_i0,...){
  return(c(x-0.005,1-rho_i0-sum(x)-0.005))
}

hin_jacob_Fix_Wgt_suppRef <-function(x,rho_i0,...){
  return(rbind(diag(1,length(x)),rep(-1,length(x))))
}

HS_ComputeLik_Fix_Wgt_suppRef<-function(x,logY,rho_i0,Z,Z_star,Sigma2,Sigma2A,
                                        var_wgt,EM_wgt){
  rho  = c(0,x,1-rho_i0-sum(x))
  
  Sigma2_wgt  = var_wgt*Sigma2
  Sigma2A_wgt = var_wgt*Sigma2A
  
  eta_ij = Z%*%rho
  mu_ijM = log(eta_ij)-Sigma2_wgt/2
  mu_ijA = log(eta_ij)-Sigma2A_wgt/2
  
  out = sum(EM_wgt*dnorm(x = logY,mean = mu_ijM,sd = sqrt(Sigma2_wgt),log = TRUE))+
        sum((1-EM_wgt)*dnorm(x = logY,mean = mu_ijA,sd = sqrt(Sigma2A_wgt),log = TRUE))
  return(out)
}

#------------------- Actual Update Functions ------------------------#
HS_UpdatePropn_Single_Wgt_suppRef <- function(x,Z,var_wgt,nCell,nG,maxIter_prop,useRho){
  #----------------------------------------#
  # Extract Info                           #
  #----------------------------------------#
  rho_i0 = x[1]
  urho_1 = x[c(2:(nCell+1))]
  logY   = x[c((nCell+2):(nCell+nG+1))]
  EM_wgt = x[c((nCell+nG+2):(nCell+2*nG+1))]
  
  Qval = ncol(Z)-1
  
  #message("Log Y: ",exp(logY[1]))
  
  #----------------------------------------#
  # Initialize Values                      #
  #----------------------------------------#
  Zm = Z[,-c(1)]
  Z_star = Zm%*%rbind(diag(rep(1,(Qval-1))),rep(-1,(Qval-1)))
  
  eta_ij   = drop(Z%*%matrix(c(0,urho_1),ncol=1))
  theta_ij = log(eta_ij)
  Sigma2M_1 = HS_Sigma2_Update_Wgt_suppRef(logY    = logY,
                                           eta_ij  = eta_ij,
                                           var_wgt = var_wgt,
                                           EM_wgt  = EM_wgt,AB_Up = FALSE)
  Sigma2A_1 = HS_Sigma2_Update_Wgt_suppRef(logY    = logY,
                                           eta_ij  = eta_ij,
                                           var_wgt = var_wgt,
                                           EM_wgt  = EM_wgt,AB_Up = TRUE)
  
  Sigma2M_1_wgt = var_wgt*Sigma2M_1
  Sigma2A_1_wgt = var_wgt*Sigma2A_1
  
  mu_ijM    = theta_ij - Sigma2M_1_wgt/2
  mu_ijA    = theta_ij - Sigma2A_1_wgt/2
  logLik_1  = sum(EM_wgt*dnorm(logY,mean=mu_ijM,sd=sqrt(Sigma2M_1_wgt),log=TRUE))+
              sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(Sigma2A_1_wgt),log=TRUE))
  
  for(k in 1:maxIter_prop){
    #----------------------------------#
    # Reset Current Values             #
    #----------------------------------#
    urho_0    = urho_1
    Sigma2M_0 = Sigma2M_1
    Sigma2A_0 = Sigma2A_1
    logLik_0  = logLik_1
    
    #----------------------------------#
    # Update Rho values                #
    #----------------------------------#
    if(useRho){
      auglagOut = auglag(par = urho_0[-c(nCell)],fn = HS_ComputeLik_Fix_Wgt_suppRef,
                         gr = HS_GradFunc_Fix_Wgt_suppRef,hin = hin_func_Fix_Wgt_suppRef,
                         hin.jac = hin_jacob_Fix_Wgt_suppRef,
                         logY = logY,rho_i0=rho_i0,Z=Z,Z_star=Z_star,Sigma2 = Sigma2M_0,
                         Sigma2A = Sigma2A_0,var_wgt=var_wgt,EM_wgt=EM_wgt,
                         control.optim=list(fnscale=-1),control.outer = list(trace=FALSE))
      
      urho_1[c(1:(nCell-1))] = auglagOut$par
      urho_1[nCell] = 1-rho_i0-sum(auglagOut$par)
      
      urho_1 = correctRho(est = urho_1,total = 1-rho_i0)
    } else {
      auglagOut = auglag(par = urho_0,fn = HS_ComputeLik_noFix_Wgt_suppRef,
                         gr = HS_GradFunc_noFix_Wgt_suppRef,hin = hin_func_noFix_Wgt_suppRef,
                         hin.jac = hin_jacob_noFix_Wgt_suppRef,
                         logY = logY,rho_i0=rho_i0,Z=Z,Z_star=Z_star,Sigma2 = Sigma2M_0,
                         Sigma2A = Sigma2A_0,var_wgt=var_wgt,EM_wgt=EM_wgt,
                         control.optim=list(fnscale=-1),control.outer = list(trace=FALSE))
      
      urho_1 = auglagOut$par
      urho_1 = correctRho_v2(urho_1)
    }

    #----------------------------------#
    # Update Computational Values      #
    #----------------------------------#
    eta_ij    = drop(Z%*%matrix(c(0,urho_1),ncol=1))
    theta_ij  = log(eta_ij)
    Sigma2M_1 = HS_Sigma2_Update_Wgt_suppRef(logY    = logY,
                                             eta_ij  = eta_ij,
                                             var_wgt = var_wgt,
                                             EM_wgt  = EM_wgt,AB_Up = FALSE)
    Sigma2A_1 = HS_Sigma2_Update_Wgt_suppRef(logY    = logY,
                                             eta_ij  = eta_ij,
                                             var_wgt = var_wgt,
                                             EM_wgt  = EM_wgt,AB_Up = TRUE)
    
    Sigma2M_1_wgt = var_wgt*Sigma2M_1
    Sigma2A_1_wgt = var_wgt*Sigma2A_1
    
    mu_ijM    = theta_ij - Sigma2M_1_wgt/2
    mu_ijA    = theta_ij - Sigma2A_1_wgt/2
    logLik_1  = sum(EM_wgt*dnorm(logY,mean=mu_ijM,sd=sqrt(Sigma2M_1_wgt),log=TRUE))+
                sum((1-EM_wgt)*dnorm(logY,mean=mu_ijA,sd=sqrt(Sigma2A_1_wgt),log=TRUE))
    
    if(max(abs(urho_1-urho_0)) < 1e-3){ break }
  }
  return(list(rho_i = urho_1,iter=k,Sigma2M_i=Sigma2M_1,Sigma2A_i=Sigma2A_1))
}

#------------------- Update Multiple Subjects ----------------------#
HS_UpdatePropn_All_Wgt_suppRef<-function(logY,Rho_init,fixed_rho,Z,maxIter_prop,
                                         var_wgt,EM_wgt,useRho){
  PropnInfo = rbind(fixed_rho,Rho_init,logY,EM_wgt)
  out = apply(X=PropnInfo,MARGIN = 2,FUN = HS_UpdatePropn_Single_Wgt_suppRef,Z=Z,var_wgt=var_wgt,
              useRho = useRho,maxIter_prop = maxIter_prop,nCell=(nrow(Rho_init)),
              nG = nrow(logY))
  
  #--- Reshape Output ---#
  propMat   = extractFunction(compList = out,element = "rho_i")
  sig2M_Mat = extractFunction(compList = out,element = "Sigma2M_i")
  sig2A_Mat = extractFunction(compList = out,element = "Sigma2A_i")
  
  return(list(propCurr = propMat,sig2M_Curr = sig2M_Mat,sig2A_Curr = sig2A_Mat))
}

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 3 - EM Weights                             #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
HS2_UpdateWgts_Single_Wgt_suppRef <-function(x,Z,var_wgt,nG,Qval){
  #-----------------------------------#
  # Extract Information               #
  #-----------------------------------#
  logY      = x[c(1:nG)]
  Rho_init  = x[c((nG+1):(nG+Qval))]
  fixed_rho = x[(nG+Qval+1)]
  Sigma2M   = x[(nG+Qval+2)]
  Sigma2A   = x[(nG+Qval+3)]
  p_m       = x[(nG+Qval+4)]
  
  #-----------------------------------#
  # Consistent Marker Gene Probs      #
  #-----------------------------------#
  eta_ij  = Z%*%matrix(c(0,Rho_init),ncol=1)
  CMmu_ij = log(eta_ij)-var_wgt*Sigma2M/2
  
  CM_LLik = dnorm(x = logY,mean = CMmu_ij,sd = sqrt(Sigma2M*var_wgt),log = TRUE)
  
  #-----------------------------------#
  # Aberrant Marker Gene Probs        #
  #-----------------------------------#
  AMmu_ij = log(eta_ij)-var_wgt*Sigma2A/2
  
  PM_LLik = dnorm(x = logY,mean = AMmu_ij,sd = sqrt(Sigma2A*var_wgt),log = TRUE)
  
  #-----------------------------------#
  # Compiling Weights                 #
  #-----------------------------------#
  EM_wgt = 1/(1+((1-p_m)/p_m)*exp(PM_LLik-CM_LLik))
  
  return(EM_wgt)
}

HS2_UpdateWgts_All_Wgt_suppRef <- function(logY,Rho_init,fixed_rho,Sigma2M,Sigma2A,
                               var_wgt,Z,p_m){
  Wgt_Info = rbind(logY,Rho_init,fixed_rho,Sigma2M,Sigma2A,p_m)
  
  EM_wgt = apply(X = Wgt_Info,MARGIN = 2,FUN = HS2_UpdateWgts_Single_Wgt_suppRef,
                 nG = nrow(logY),Qval = nrow(Rho_init),Z=Z,var_wgt=var_wgt)
  
  return(EM_wgt)
}

RP_UpdatePM_Wgt_suppRef <- function(nG,EM_wgt){
  return(colSums(EM_wgt)/nG)
}