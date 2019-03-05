##########################################################################
# PROJECT 3: Fitting Algorithm                                           #
##########################################################################
# PROGRAM NAME:                                                          #
#   P3_FittingAlgorithm.R                                                #
# PROGRAMMER:                                                            #
#   Douglas Roy Wilson, Jr.                                              #
# DATE CREATED:                                                          #
#   05/22/2017                                                           #
# LAST EDIT:                                                             #
#   05/22/2017                                                           #
# VERSION:                                                               #
#   R-3.3.1                                                              #
#------------------------------------------------------------------------#
# DESCRIPTION:                                                           #  
#   The following functions use the modular fit pieces to describe       #
#   and conduct the proposed fit algorithm for project 3 - examine       #
#   immune cell type abundance from within a bulk tumor tissue.          #
##########################################################################

#------------------------------------------------------------------------#
# SUPPORTING FUNCTIONS                                                   #
#------------------------------------------------------------------------#
# For testing purposes only.
IQRFun <- function(x,group){
  out = tapply(X = x,INDEX = group,FUN = IQR)
  return(out)
}

medianFun <- function(x,group){
  out = tapply(X = x,INDEX = group,FUN = median)
  return(out)
}


#------------------------------------------------------------------------#
# ASICeD_fit_ABpool_CP()                                                 #
#------------------------------------------------------------------------# 
# INPUT:
#========================================================================#
# VARIABLE       |   RANGE    | DESCRIPTION                              #
#========================================================================#
#       Y        |  [0,Inf]   | TReC and Genomic Feature length          #
#                |            | corrected expression for each bulk       #
#                |            | tumor sample. Matrix of size nG x nS.    #
#----------------|------------|------------------------------------------#
#       X        |  [0,Inf]   | TReC and Genomic Feature length          #
#                |            | corrected expression for each pure       #
#                |            | reference sample. Matrix of size         #
#                |            | nG x nP                                  #
#----------------|------------|------------------------------------------#
#   cellType     | {1,2,..,K} | Vector of length nP describing cell-     #
#                |    -or-    | type of each column of matrix X. Order   #
#                |{CT1,..,CTK}| must match that found in the matrix X.   #
#----------------|------------|------------------------------------------#
#    fixedCT     | {1,2,..,K} | Single value matching one of the labels  #
#                |    -or-    | in cellType that identifies the tumor    #
#                |{CT1,..,CTK}| samples used for "references." If NULL,  #
#                |            | all cell type proportions are estimated  #
#                |            | and it is assumed all cell types have    #
#                |            | purified references.                     #
#----------------|------------|------------------------------------------#
#  fixedCT_rho   | [0.0,1.0]  | Vector of length nS recording the tumor  #
#                |            | purity of each mixture sample. Order     #
#                |            | must match column order in Y.            #
#----------------|------------|------------------------------------------#
#     PurCO      | [0.0,1.0]  | In order to initialize a pure sample     #
#                |            | reference for the fixed cell type, we    #
#                |            | will incorporate all pure samples of this#
#                |            | type as well as all mixtures with a      #
#                |            | fixedCT_rho>PurCO. Recommended >= .90    #
#----------------|------------|------------------------------------------#
#   borrow4SD    | TRUE/FALSE | Indicator of whether or not to borrow    #
#                |            | information across cell types in order   #
#                |            | to estimate Std Dev in pure references.  #
#----------------|------------|------------------------------------------#
#    maxIter     |  [1,500]   | Maximum number of iteration used for     #
#                |            | updating model parameters.               #
#----------------|------------|------------------------------------------#
#    FitTech     | {"PH","R"} | Fit Technique used to update proportions #
#                |            | (PH - Post Hoc, R - constrOptim).        #
#========================================================================#
#
# OUTPUT:
#========================================================================#
# VARIABLE       | DESCRIPTION                                           #
#========================================================================#
# $IC_Abundance  | Estimated immune cell type abundances (also contains  #
#                | tumor purity estimates provided by user). Matrix of   #
#                | size nS x (K+1).                                      #
#----------------|-------------------------------------------------------#
# $LogExp_Mean   | Estimated mean expression profile for each cell type  #
#                | on the log scale. This includes the tumor/fixed cell  #
#                | type expression profile. Matrix of size nG x (K+1).   #
#----------------|-------------------------------------------------------#
# $LogExp_StdDev | Estimated standard deviation of expression for each   #
#                | cell type on the log scale. This includes the tumor/  #
#                | fixed cell type. Matrix of size nG x (K+1).           #
#----------------|-------------------------------------------------------#
# $Fixed_CellType| The label of the tumor/fixed cell type specified by   #
#                | the user. Equals "None" if the user did not provide   #
#                | a fixed cell type.                                    #
#----------------|-------------------------------------------------------#
# $Convg_Ind     | Indicator of convergence status of the fit procedure: #
#                |     0 - Convergence                                   #
#                |     1 -                                               #
#                |     2 -                                               #
#========================================================================#

ICeDT_Wgt_suppRef<-function(Y,fixedCT_rho,useRho=FALSE,RefMat,RefVar,varLog=FALSE,
                                limitWgt = TRUE,limitQuant=c(0.10,0.90),userWgt=NULL,
                                useIQR=FALSE,maxIter_prop = 100,
                                maxIter_PP=100,RhoConv_CO = 1e-4, Subj_CO){
  #-----------------------------------------------------#
  # Check Input                                         #
  #-----------------------------------------------------#
  if(nrow(Y)!=nrow(RefMat)){
    stop("Mixture Expresion and Reference Matrix have different
         numbers of genes!")
  }
  
  if(!identical(rownames(Y),rownames(RefMat))){
    stop("Mixture Expressions and Reference Matrix genes
         not aligned. Check Row Order!")
  }
  
  if(ncol(Y)!=length(fixedCT_rho)){
    stop("Mixture Expressions and Tumor Purity
         mismatch - Different Number of Subjects!")
  } else if(!all(colnames(Y)==names(fixedCT_rho))){
    stop("Mixture Expression and Tumor Purity 
         mismatch! -- Different order of subjects or labels incorrect.")
  }
  
  if(!identical(dim(RefMat),dim(RefVar))){
    stop("Reference Expression and Reference Variance mismatch!")
  }
  
  if(!identical(rownames(RefMat),rownames(RefVar))||
     !identical(colnames(RefMat),colnames(RefVar))){
    stop("Reference Expression and Variance matrices are 
          not matched!")
  }
  
  #-----------------------------------------------------#
  # EXTRACTING SAMPLE INFO                              #
  #-----------------------------------------------------#
  sortCT = colnames(RefMat)
  
  if(any(duplicated(sortCT))){
    stop("Reference Matrix has more than one column with
         the same cell type label!")
  }
  
  nG  = nrow(RefMat)
  nS  = ncol(Y)
  Qp1 = length(sortCT)+1
  
  #-----------------------------------------------------#
  # ORDER REFERENCE MATRIX                              #
  #-----------------------------------------------------#
  # Rearrange order so that fixed cell type (tumor) is
  # first.
  RefMat_1 = cbind(rep(0,nG),RefMat[,sortCT])
  colnames(RefMat_1) = c("Dummy_Tumor",sortCT)
  
  #-----------------------------------------------------#
  # INITIALIZATION                                      #
  #-----------------------------------------------------#
  Z_1      = RefMat_1
  
  #----------- Estimating Variance Weight --------------#
  if(useIQR%in%c(0,1,2,3,4)&varLog==FALSE){
    stop("Recompute Variance on Log-Scale!")
  } else if(useIQR==5&varLog==TRUE){
    stop("Recompute Variance on normal scale!")
  }
  
  if(useIQR==0){
    maxVar = apply(X = RefVar,MARGIN = 1,FUN = max)
    var_wgt = maxVar/median(maxVar)
  } else if(useIQR==6){
    var_wgt = userWgt
  }
  
  if(limitWgt){
    quants = quantile(var_wgt,prob=limitQuant)
    w10 = which(var_wgt<quants[1])
    w90 = which(var_wgt>quants[2])
    var_wgt[w10] = quants[1]
    var_wgt[w90] = quants[2]
  }
  
  #----------- Mixture Sample Proportions --------------#
  P_Est = apply(X = rbind(fixedCT_rho,Y),MARGIN = 2,FUN = lmInit_Wgt_suppRef,
                Zm=Z_1[,-c(1)],Zt=Z_1[,1])
  Rho_1 = P_Est
  
  #----------- Aberrant Profile Initial ----------------#
  SigmaInit = apply(X = rbind(Y,Rho_1,fixedCT_rho),MARGIN=2,
                    FUN=AbProf_Init_Wgt_suppRef,Z = Z_1,var_wgt=var_wgt,
                    nG = nG,Qval = (ncol(Z_1)-1))
  
  Sigma2M_1 = SigmaInit[1,]
  Sigma2A_1 = SigmaInit[2,]
  
  #---- Additional Values ----#
  Pm_1  = rep(0.5,nS)
  
  #-----------------------------------------------------#
  # IMPLEMENT ALGORITHMIC FIT                           #
  #-----------------------------------------------------#
  Rho_0     = Rho_1
  Sigma2M_0 = Sigma2M_1
  Sigma2A_0 = Sigma2A_1
  
  Z_0       = Z_1
  
  Pm_0  = Pm_1

  #------ Update Proportions + -----#
  PropPlus_Out   = PropPlus_Update_Wgt_suppRef(Y = Y,Rho_0 = Rho_0,fixedCT_rho = fixedCT_rho,useRho=useRho,
                                               Sigma2M_0 = Sigma2M_0,Sigma2A_0 = Sigma2A_0,var_wgt = var_wgt,
                                               Z_0 = Z_0,Pm_0 = Pm_0,maxIter_PP = maxIter_PP,
                                               maxIter_prop = maxIter_prop,nG = nG,
                                               RhoConv_CO = RhoConv_CO,Subj_CO = Subj_CO)
  
  Rho_1     = PropPlus_Out$Rho
  Sigma2M_1 = PropPlus_Out$Sigma2M
  Sigma2A_1 = PropPlus_Out$Sigma2A
  Pm_1      = PropPlus_Out$Pm
  
  EM_wgt =  HS2_UpdateWgts_All_Wgt_suppRef(logY = log(Y),Rho_init = Rho_1,fixed_rho = fixedCT_rho,
                                           Sigma2M = Sigma2M_1,Sigma2A = Sigma2A_1,var_wgt = var_wgt,
                                           Z = Z_1,p_m = Pm_1)
  
  #-----------------------------------------------------#
  # OUTPUT                                              #
  #-----------------------------------------------------#
  if(useRho){
    IC_Abundance  = cbind(fixedCT_rho,t(Rho_1))
    colnames(IC_Abundance) = c("Dummy_Tumor",sortCT)
    rownames(IC_Abundance) = colnames(Y)
    
    Fixed_CellType = "Dummy_Tumor"
  } else {
    fixedCT_rho_est = 1-colSums(Rho_1)
    
    IC_Abundance  = cbind(fixedCT_rho_est,t(Rho_1))
    colnames(IC_Abundance) = c("Dummy_Tumor",sortCT)
    rownames(IC_Abundance) = colnames(Y)
    
    Fixed_CellType = "Dummy_Tumor"
  }
  
  outList = list(IC_Abundance   = IC_Abundance,
                 Fixed_CellType = Fixed_CellType,
                 Sigma2M        = Sigma2M_1,
                 Sigma2A        = Sigma2A_1,
                 var_wgt        = var_wgt,
                 Z              = Z_1,
                 var_wgt        = var_wgt,
                 P_Consistent   = Pm_1,
                 PP_Consistent  = EM_wgt)
  
  return(outList)
}

PropPlus_Update_Wgt_suppRef<- function(Y,Rho_0,fixedCT_rho,useRho,
                                       Sigma2M_0,Sigma2A_0,var_wgt,
                                       Z_0,Pm_0,
                                       maxIter_PP,maxIter_prop,nG,
                                       RhoConv_CO,Subj_CO){
  logY = log(Y)
  
  Rho_t1     = Rho_0
  Sigma2M_t1 = Sigma2M_0
  Sigma2A_t1 = Sigma2A_0
  Pm_t1      = Pm_0
  
  for(j in 1:maxIter_PP){
    #----  Reset the Param   ----#
    Rho_t0     = Rho_t1
    Sigma2M_t0 = Sigma2M_t1
    Sigma2A_t0 = Sigma2A_t1
    Pm_t0      = Pm_t1
    
    #---- Update EM Weights  ----#
    EM_wgt = HS2_UpdateWgts_All_Wgt_suppRef(logY = logY,Rho_init = Rho_t0,fixed_rho = fixedCT_rho,
                                            Sigma2M = Sigma2M_t0,Sigma2A = Sigma2A_t0,var_wgt = var_wgt,
                                            Z = Z_0,p_m = Pm_t0)
    
    #---- Update Proportions ----#
    PropP_t1  = HS_UpdatePropn_All_Wgt_suppRef(logY = logY,Rho_init = Rho_t0,fixed_rho = fixedCT_rho,
                                               Z = Z_0,maxIter_prop = maxIter_prop,EM_wgt = EM_wgt,
                                               useRho = useRho,var_wgt = var_wgt)
    
    Rho_t1     = PropP_t1$propCurr
    Sigma2M_t1 = PropP_t1$sig2M_Curr
    Sigma2A_t1 = PropP_t1$sig2A_Curr
    
    #---- Update Cons. MarkP ----#
    Pm_t1 = RP_UpdatePM_Wgt_suppRef(nG = nG,EM_wgt = EM_wgt)
    
    #---- Check Convergence  ----#
    Rho_Diff = abs(Rho_t0-Rho_t1)
    max_Diff = apply(X = Rho_Diff,MARGIN = 2,FUN = max)
    
    if(sum(max_Diff<RhoConv_CO)>=Subj_CO){
      break
    }
    
    if(j%%15){
      message("Update Iteration ",j,"Complete")
    }
    
  }
  
  return(list(Rho = Rho_t1,Sigma2M = Sigma2M_t1,Sigma2A = Sigma2A_t1,Pm = Pm_t1,Iter = j))
}