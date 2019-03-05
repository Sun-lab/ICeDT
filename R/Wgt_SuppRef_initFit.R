#####################################################################
# INITIALIZATION FUNCTIONS                                          #
#####################################################################
# PROGRAM NAME:                                                     #
#   P3_InitFit.R                                                    #
# PROGRAMMER:                                                       #
#   Douglas Roy Wilson, Jr.                                         #
# DATE CREATED:                                                     #
#   05/22/2017                                                      #
# LAST EDIT:                                                        #
#   05/22/2017                                                      #
# VERSION:                                                          #
#   R-3.3.1                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   Contains the code necessary to implement the initial estimates  #
#   for several important parameters in the ASICeD model fit.       #
#####################################################################

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
lmInit_Wgt_suppRef<-function(yr,Zm,Zt){
  eta = yr[1]
  y   = yr[-c(1)]
  
  if(is.null(Zt)){
    offVal = rep(0,length(y))
  } else {
    offVal = Zt*eta
  }
  
  initMod = lm(y~Zm,offset=offVal)
  
  rho  = coef(initMod)[-1]
  wneg = which(rho<0.01)
  rho[wneg] = 0.01
  rho = (rho/sum(rho))*(1-eta)
  
  return(rho)
}

#-------------------------------------------------------------------#
# UTILITY FUNCTION 2: Aberrant Profile Init                         #
#-------------------------------------------------------------------#
AbProf_Init_Wgt_suppRef <- function(x,Z,var_wgt,nG,Qval){
  Y      = x[c(1:nG)]
  logY   = log(Y)
  rho    = x[c((nG+1):(nG+Qval))]
  rho_i0 = x[(nG+Qval+1)]
  
  #eta_ij = Z%*%c(rho_i0,rho)
  eta_ij = Z%*%c(0,rho)
  
  resid = Y-eta_ij
  Q3val = quantile(abs(resid),probs = c(0.75))
  g2use = which(abs(resid)>Q3val)
  
  init_val = c(0,0)
  
  init_val[1] = HS_Sigma2_Update_Wgt_suppRef(logY = logY[-c(g2use)],eta_ij = eta_ij[-c(g2use)],
                                             EM_wgt = rep(1,(nG-length(g2use))),AB_Up = FALSE,var_wgt = var_wgt)
  init_val[2] = HS_Sigma2_Update_Wgt_suppRef(logY = logY[c(g2use)],eta_ij = eta_ij[c(g2use)],
                                             EM_wgt = rep(0,length(g2use)),AB_Up = TRUE,var_wgt = var_wgt)
  
  
  return(init_val)
}
