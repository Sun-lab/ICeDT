#####################################################################
# CIBERSORT FLOW CYTOMETRY DATA: ICeDT Fit with LM22 reference      #
#####################################################################
# PROGRAM NAME:                                                     #
#   FlowCytometryFit_LM22Matrix.R                                   #
# PROGRAMMER:                                                       #
#   Douglas Roy Wilson, Jr.                                         #
# DATE CREATED:                                                     #
#   2/22/2018                                                       #
# LAST EDIT:                                                        #
#   2/22/2018                                                       #
# VERSION:                                                          #
#   R-3.3.1                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   After altering the ICeD-T programs to allow model fit with a    #
#   predefined reference expression matrix, we are refitting to the #
#   CIBERSORT Flow Cytometry data using their LM22 matrix.          #
#####################################################################

library(alabama)

useWgts   = TRUE
useTarget = TRUE
wgtOpt    = 0

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 0 - LIBRARY and NECESSARY CODE             #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
library(preprocessCore)

meanFun <- function(x,group){
  out = tapply(X = x,INDEX = group,FUN = mean)
  return(out)
}

varFun <- function(x,group){
  out = tapply(X = x,INDEX = group,FUN = var)
  return(out)
}

IQRFun <- function(x,group){
  out = tapply(X = x,INDEX = group,FUN = IQR)
  return(out)
}

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 1 - Fit CIBERSORT                          #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
setwd("D://DougData/Documents/Dissertation/Paper 3 - Immune Cell Deconv/Sun_Requested_Analysis/DrSun_GrantFigures/data/FlowCytometry/")

#load ground truth
flow = read.table("PBMCs-Fig3a-Flow-Cytometry.txt", sep="\t", row.names=1, header=T)
dim(flow)
flow[1:2,]

# Load Web - App Cibersort:
CSORT_FlowCyto = read.table(file = "../../programs/FlowCytometry/Fits/CIBERSORT.Output_FlowCyto.csv",
                            header=TRUE,sep = ",",row.names = 1)
results = CSORT_FlowCyto[,1:22]

#clean results for comparison
results_clean = cbind(results[,c(1:2,4:7,10)],apply(results[,11:12],1,sum),apply(results[,13:16],1,sum))

#normalize to 1
results_clean = 100*results_clean / apply(results_clean,1,sum)

colnames(results_clean)[8]="NK cells"
colnames(results_clean)[9]="Monocytes"

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 2 - PROCESS DATA                           #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
exprDat = read.table(file = "PBMCs-Fig3a-HumanHT-12-V4.txt",header = TRUE,sep = "\t")
exprDat_v2 = exprDat[,-c(1)]
rownames(exprDat_v2) = exprDat[,1]

### CIBERSORT Load
lm22 = read.table("LM22-ref-sample_logged.txt", , sep = "\t",
                  header = TRUE, as.is = TRUE)
dim(lm22)
lm22[1:2,1:2]
lm22 = data.matrix(lm22)

sam = read.table("LM22-ref-sample_info.txt", , sep = "\t",
                 header = TRUE, as.is = TRUE)
dim(sam)
sam[1:2,]

table(sam$label)

lm22_new = exp(lm22)

g2use = intersect(rownames(lm22_new),rownames(exprDat_v2))
lm22_new   = lm22_new[g2use,]
exprDat_v2 = exprDat_v2[g2use,]

### Option 1: Quantile normalize the reference data, match mixture to that
### Option 2: Quantile normalize within major cell groups, then quantile normalize
###           mixture to that. 

### Let's go with option 1 so that we can see what is happening.
lm22_qnorm = normalize.quantiles(lm22_new)
rownames(lm22_qnorm) = rownames(lm22_new)
colnames(lm22_qnorm) = colnames(lm22_new)
if(useTarget){
  target = normalize.quantiles.determine.target(x=lm22_qnorm)
  exprDat_qnorm = normalize.quantiles.use.target(x = as.matrix(exprDat_v2),target = target)
  rownames(exprDat_qnorm) = rownames(exprDat_v2)
  colnames(exprDat_qnorm) = colnames(exprDat_v2)
} else {
  exprDat_qnorm = normalize.quantiles(as.matrix(exprDat_v2))
  rownames(exprDat_qnorm) = rownames(exprDat_v2)
  colnames(exprDat_qnorm) = colnames(exprDat_v2)
}

# restrict to LM22 genes:
lm22_ref = read.table(file = "LM22.txt",header = TRUE,sep = "\t",row.names = 1)
g2use = rownames(lm22_ref)

g2use_fin = intersect(g2use,rownames(lm22_qnorm))
g2use_fin = intersect(g2use_fin,rownames(exprDat_qnorm))

lm22_ref2  = lm22_ref[g2use_fin,]
MixDat_Use = exprDat_qnorm[g2use_fin,]
PurDat_Use = lm22_qnorm[g2use_fin,]
cellType   = sam$label

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 3 - FIT with CIBERSORT MATRIX              #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
RefMat = lm22_ref2

CT_VAR = t(apply(X = log(PurDat_Use),MARGIN = 1,FUN = varFun,
                 group = cellType))

CT_VAR = CT_VAR[,c("B naive","B memory","Plasma","CD8+ T","CD4+ T","CD4+ T memory resting",
                   "CD4+ T memory activated","Tfh","Tregs","T gd","NK resting","NK activated","Monocytes",
                   "Macrophages M0","Macrophages M1","Macrophages M2","Dendritic resting","Dendritic activated",
                   "Mast resting","Mast activated","Eosinophils","Neutrophils")]

colnames(CT_VAR) = colnames(lm22_ref2)

if(useWgts){
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_initFit.R")
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_FittingAlgorithm.R")
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_UpdateFunctions.R")
  
  fit = ASICeD_fit(Y=MixDat_Use,fixedCT_rho=rep(0,ncol(MixDat_Use)),useRho = FALSE,RefMat = as.matrix(RefMat),RefVar = CT_VAR,varLog = TRUE,
                   limitWgt = TRUE,limitQuant = c(0.15,0.85),userWgt = NULL,nu = 0.0001,useIQR = wgtOpt,maxIter = 500,maxIter_prop = 500,
                   maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(MixDat_Use));
  
  setwd("D://DougData/Documents/Dissertation/Paper 3 - Immune Cell Deconv/Sun_Requested_Analysis/DrSun_GrantFigures/programs/FlowCytometry/Fits/")
  fnm = sprintf("./CIBERSORT_useWgts_%s/ICeDT%s_target_%s_fit.RData",
                useWgts,wgtOpt,useTarget)
  
  save(fit,file=fnm)
} else {
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_initFit.R")
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_FittingAlgorithm.R")
  source("../../programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_UpdateFunctions.R")
  
  fit = ASICeD_fit(Y=MixDat_Use,fixedCT_rho=rep(0,ncol(MixDat_Use)),useRho = FALSE,
                   RefMat = as.matrix(RefMat),maxIter = 500,maxIter_prop = 500,
                   maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(MixDat_Use));
  
  setwd("D://DougData/Documents/Dissertation/Paper 3 - Immune Cell Deconv/Sun_Requested_Analysis/DrSun_GrantFigures/programs/FlowCytometry/Fits/")
  fnm = sprintf("./CIBERSORT_useWgts_%s/ICeDT%s_target_%s_fit.RData",
                useWgts,wgtOpt,useTarget)
  
  save(fit,file=fnm)
}

setwd("D://DougData/Documents/Dissertation/Paper 3 - Immune Cell Deconv/Sun_Requested_Analysis/DrSun_GrantFigures/programs/FlowCytometry/Fits/")

fnm = sprintf("./CIBERSORT_useWgts_%s/ICeDT%s_target_%s_fit.RData",
              useWgts,wgtOpt,useTarget)

load(fnm)

resDj = fit$IC_Abundance[,-c(1)]

results_clean_dj = cbind(resDj[,c(1:2,4:7,10)],apply(resDj[,11:12],1,sum),apply(resDj[,13:16],1,sum))
#normalize to 1
results_clean_dj = 100*results_clean_dj / apply(results_clean_dj,1,sum)

colnames(results_clean_dj)[8]="NK cells"
colnames(results_clean_dj)[9]="Monocytes"

### Correlation Pre CT Correction
sum((results_clean_dj-flow)^2)

fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_ICeDT%s_target_%s.pdf",
              useWgts,wgtOpt,useTarget)
pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_dj)[2]){
  plot(flow[,i],results_clean_dj[,i], main=colnames(results_clean_dj)[i],xlab="Flow cytometry (%)", ylab="ICeD-T (%)", las=1)	
  abline(lm(results_clean_dj[,i]~flow[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow[,i],results_clean_dj[,i]),2)),paste("P=", round(cor.test(flow[,i],results_clean_dj[,i])$p.value,4))))
  i = i + 1
}

dev.off()

### Correlation Post CT COrrection
ctSize = c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.42,1.4)

results_clean_dj_ctCorr = t(t(results_clean_dj)/ctSize)
results_clean_dj_ctCorr = 100*(results_clean_dj_ctCorr/rowSums(results_clean_dj_ctCorr))

results_clean_ctCorr = t(t(results_clean)/ctSize)
results_clean_ctCorr = 100*(results_clean_ctCorr/rowSums(results_clean_ctCorr))

fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_ICeDT%s_target_%s_ctCorr.pdf",
              useWgts,wgtOpt,useTarget)
pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_dj_ctCorr)[2]){
  plot(flow[,i],results_clean_dj_ctCorr[,i], main=colnames(results_clean_dj_ctCorr)[i],xlab="Flow cytometry (%)", ylab="ICeD-T (%)", las=1)	
  abline(lm(results_clean_dj_ctCorr[,i]~flow[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow[,i],results_clean_dj_ctCorr[,i]),2)),paste("P=", round(cor.test(flow[,i],results_clean_dj_ctCorr[,i])$p.value,4))))
  i = i + 1
}

dev.off()

sum((results_clean_dj_ctCorr-flow)^2)
sum((results_clean_ctCorr-flow)^2)

cor(c(results_clean_dj_ctCorr),c(as.matrix(flow)))
cor(c(results_clean_ctCorr),c(as.matrix(flow)))

### Plot the Web CIBERSORT (just in case)
fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_WebCibersort.pdf",
              useWgts)

pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_ctCorr)[2]){
  plot(flow[,i],results_clean_ctCorr[,i], main=colnames(results_clean_ctCorr)[i],xlab="Flow cytometry (%)", ylab="CIBERSORT (%)", las=1)	
  abline(lm(results_clean_ctCorr[,i]~flow[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow[,i],results_clean_ctCorr[,i]),2)),paste("P=", round(cor.test(flow[,i],results_clean_ctCorr[,i])$p.value,4))))
  i = i + 1
}

dev.off()


#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 4 - EPIC FIT                               #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
### EPIC GENES
FList = c("BANK1","CD79A", "CD79B", "FCER2", "FCRL2", "FCRL5", "MS4A1", "PAX5", "POU2AF1", "STAP1", "TCL1A",
          "ADAM33", "CLDN11", "COL1A1", "COL3A1", "COL14A1", "CRISPLD2", "CXCL14", "DPT", "F3", "FBLN1", "ISLR", "LUM", "MEG3", "MFAP5", "PRELP", "PTGIS", "SFRP2", "SFRP4", "SYNPO2", "TMEM119",
          "ANKRD55", "DGKA", "FOXP3", "GCNT4", "IL2RA", "MDS2", "RCAN3", "TBC1D4", "TRAT1",
          "CD8B", "HAUS3", "JAKMIP1","NAA16", "TSPYL1",
          "CDH5", "CLDN5", "CLEC14A", "CXorf36", "ECSCR", "F2RL3", "FLT1", "FLT4", "GPR4", "GPR182", "KDR", "MMRN1", "MMRN2", "MYCT1", "PTPRB", "RHOJ", "SLCO2A1", "SOX18", "STAB2", "VWF",
          "APOC1", "C1QC", "CD14", "CD163", "CD300C", "CD300E", "CSF1R", "F13A1", "FPR3", "HAMP", "IL1B", "LILRB4", "MS4A6A", "MSR1", "SIGLEC1", "VSIG4",
          "CD33", "CD300C", "CD300E", "CECR1", "CLEC6A", "CPVL", "EGR2", "EREG", "MS4A6A", "NAGA", "SLC37A2",
          "CEACAM3", "CNTNAP3", "CXCR1", "CYP4F3", "FFAR2", "HIST1H2BC", "HIST1H3D", "KY", "MMP25", "PGLYRP1", "SLC12A1", "TAS2R40",
          "CD160", "CLIC3", "FGFBP2", "GNLY", "GNPTAB", "KLRF1", "NCR1", "NMUR1", "S1PR5", "SH2D1B",
          "BCL11B", "CD5", "CD28", "IL7R", "ITK", "THEMIS", "UBASH3A")

table(FList%in%rownames(lm22_ref2))

### FITTING EPIC
library(EPIC)

refProfiles.1 = t(apply(X = PurDat_Use,MARGIN = 1,FUN = meanFun,
                        group = cellType))
refProfiles.2 = lm22_ref2

refProfiles.var1 = t(apply(X = PurDat_Use,MARGIN = 1,FUN = varFun,
                           group = cellType))
refProfiles.var1.reOrd = refProfiles.var1[,c("B naive","B memory","Plasma","CD8+ T","CD4+ T","CD4+ T memory resting",
                                             "CD4+ T memory activated","Tfh","Tregs","T gd","NK resting","NK activated","Monocytes",
                                             "Macrophages M0","Macrophages M1","Macrophages M2","Dendritic resting","Dendritic activated",
                                             "Mast resting","Mast activated","Eosinophils","Neutrophils")]
colnames(refProfiles.var1.reOrd) = colnames(refProfiles.2)

refProfiles.var2 = t(apply(X = PurDat_Use,MARGIN = 1,FUN = IQRFun,
                           group = cellType))
refProfiles.var2.reOrd = refProfiles.var2[,c("B naive","B memory","Plasma","CD8+ T","CD4+ T","CD4+ T memory resting",
                                             "CD4+ T memory activated","Tfh","Tregs","T gd","NK resting","NK activated","Monocytes",
                                             "Macrophages M0","Macrophages M1","Macrophages M2","Dendritic resting","Dendritic activated",
                                             "Mast resting","Mast activated","Eosinophils","Neutrophils")]
colnames(refProfiles.var2.reOrd) = colnames(refProfiles.2)

sigGenes = FList
sigGenes2 = intersect(rownames(lm22_ref2),rownames(refProfiles.2))

### Use QNORM - Profile 1 
EPIC_fit_p1v1 = EPIC(bulk = MixDat_Use,reference = list(refProfiles=refProfiles.1,sigGenes=sigGenes,
                                                        refProfiles.var=refProfiles.var1),scaleExprs = FALSE)
EPIC_fit_p1v2 = EPIC(bulk = MixDat_Use,reference = list(refProfiles=refProfiles.1,sigGenes=sigGenes,
                                                        refProfiles.var=refProfiles.var2),scaleExprs = FALSE)

### Use QNORM - Profile 2 
EPIC_fit_p2v1 = EPIC(bulk = MixDat_Use,list(refProfiles=as.matrix(refProfiles.2),sigGenes=sigGenes2,
                                            refProfiles.var=as.matrix(refProfiles.var1.reOrd[sigGenes2,])),scaleExprs = FALSE)
EPIC_fit_p2v2 = EPIC(bulk = MixDat_Use,reference = list(refProfiles=as.matrix(refProfiles.2),sigGenes=sigGenes2,
                                                        refProfiles.var=as.matrix(refProfiles.var2.reOrd[sigGenes2,])),scaleExprs = FALSE)

### Correct
EPIC_p1v1 = EPIC_fit_p1v1$mRNAProportions[,c("B naive","B memory","Plasma","CD8+ T","CD4+ T","CD4+ T memory resting",
                                             "CD4+ T memory activated","Tfh","Tregs","T gd","NK resting","NK activated","Monocytes",
                                             "Macrophages M0","Macrophages M1","Macrophages M2","Dendritic resting","Dendritic activated",
                                             "Mast resting","Mast activated","Eosinophils","Neutrophils")]
rclean_EPICp1v1 = cbind(EPIC_p1v1[,c(1:2,4:7,10)],apply(EPIC_p1v1[,11:12],1,sum),apply(EPIC_p1v1[,13:16],1,sum))
rclean_EPICp1v1 = 100*rclean_EPICp1v1 / apply(rclean_EPICp1v1,1,sum)

EPIC_p1v2 = EPIC_fit_p1v2$mRNAProportions[,c("B naive","B memory","Plasma","CD8+ T","CD4+ T","CD4+ T memory resting",
                                             "CD4+ T memory activated","Tfh","Tregs","T gd","NK resting","NK activated","Monocytes",
                                             "Macrophages M0","Macrophages M1","Macrophages M2","Dendritic resting","Dendritic activated",
                                             "Mast resting","Mast activated","Eosinophils","Neutrophils")]
rclean_EPICp1v2 = cbind(EPIC_p1v2[,c(1:2,4:7,10)],apply(EPIC_p1v2[,11:12],1,sum),apply(EPIC_p1v2[,13:16],1,sum))
rclean_EPICp1v2 = 100*rclean_EPICp1v2 / apply(rclean_EPICp1v2,1,sum)

EPIC_p2v1 = EPIC_fit_p2v1$mRNAProportions
rclean_EPICp2v1 = cbind(EPIC_p2v1[,c(1:2,4:7,10)],apply(EPIC_p2v1[,11:12],1,sum),apply(EPIC_p2v1[,13:16],1,sum))
rclean_EPICp2v1 = 100*rclean_EPICp2v1 / apply(rclean_EPICp2v1,1,sum)

EPIC_p2v2 = EPIC_fit_p2v2$mRNAProportions
rclean_EPICp2v2 = cbind(EPIC_p2v2[,c(1:2,4:7,10)],apply(EPIC_p2v2[,11:12],1,sum),apply(EPIC_p2v2[,13:16],1,sum))
rclean_EPICp2v2 = 100*rclean_EPICp2v2 / apply(rclean_EPICp2v2,1,sum)

colnames(rclean_EPICp1v1)[8] = colnames(rclean_EPICp1v2)[8] = colnames(rclean_EPICp2v1)[8] = colnames(rclean_EPICp2v2)[8] ="NK cells"
colnames(rclean_EPICp1v1)[9] = colnames(rclean_EPICp1v2)[9] = colnames(rclean_EPICp2v1)[9] = colnames(rclean_EPICp2v2)[9] ="Monocytes"

sum((rclean_EPICp1v1-flow)^2)
sum((rclean_EPICp1v2-flow)^2)
sum((rclean_EPICp2v1-flow)^2)
sum((rclean_EPICp2v2-flow)^2)

ctSize = c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.42,1.4)

rclean_Ep1v1_ct = t(t(rclean_EPICp1v1)/ctSize)
rclean_Ep1v1_ct = 100*(rclean_Ep1v1_ct/rowSums(rclean_Ep1v1_ct))

rclean_Ep1v2_ct = t(t(rclean_EPICp1v2)/ctSize)
rclean_Ep1v2_ct = 100*(rclean_Ep1v2_ct/rowSums(rclean_Ep1v2_ct))

rclean_Ep2v1_ct = t(t(rclean_EPICp2v1)/ctSize)
rclean_Ep2v1_ct = 100*(rclean_Ep2v1_ct/rowSums(rclean_Ep2v1_ct))

rclean_Ep2v2_ct = t(t(rclean_EPICp2v2)/ctSize)
rclean_Ep2v2_ct = 100*(rclean_Ep2v2_ct/rowSums(rclean_Ep2v2_ct))

sum((rclean_Ep1v1_ct-flow)^2)
sum((rclean_Ep1v2_ct-flow)^2)
sum((rclean_Ep2v1_ct-flow)^2)
sum((rclean_Ep2v2_ct-flow)^2)

cor(c(rclean_Ep1v1_ct),c(as.matrix(flow)))
cor(c(rclean_Ep1v2_ct),c(as.matrix(flow)))
cor(c(rclean_Ep2v1_ct),c(as.matrix(flow)))
cor(c(rclean_Ep2v2_ct),c(as.matrix(flow)))

### EPIC Plots
prof = c("p1v1","p1v2","p2v1","p2v2")
for(j in prof){
  pdf(file = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_EPIC_%s.pdf",useWgts,j), height=9, width=9)
  
  par( mfrow = c( 3, 3 ) )
  
  eval(parse(text=sprintf("results_clean_dj = rclean_E%s_ct",j)))
  
  i = 1
  while(i <= dim(results_clean_dj)[2]){
    plot(flow[,i],results_clean_dj[,i], main=colnames(results_clean_dj)[i],xlab="Flow cytometry (%)", ylab="EPIC (%)", las=1)	
    abline(lm(results_clean_dj[,i]~flow[,i]),col="red")
    legend("topleft", legend=c(paste("R=",round(cor(flow[,i],results_clean_dj[,i]),2)),paste("P=", round(cor.test(flow[,i],results_clean_dj[,i])$p.value,4))))
    i = i + 1
  }
  
  dev.off()
  
}

### USes EPIC_p2v2 because it has the highest correlation and (of those with
### the highest correlation), it has the lowest error.

#-------------------------------------------------------------------#
# Grouping Similar Cell Types                                       #
#-------------------------------------------------------------------#
flow_ed = cbind(rowSums(flow[,1:2]),flow[,3,drop=FALSE],rowSums(flow[,4:6]),flow[,7:9])
colnames(flow_ed)[c(1,3)] = c("B-cells","CD4")

results_clean_dj_cted = cbind(rowSums(results_clean_dj_ctCorr[,1:2]),results_clean_dj_ctCorr[,3,drop=FALSE],rowSums(results_clean_dj_ctCorr[,4:6]),results_clean_dj_ctCorr[,7:9])
colnames(results_clean_dj_cted)[c(1,3)] = c("B-cells","CD4")

results_clean_cted = cbind(rowSums(results_clean_ctCorr[,1:2]),results_clean_ctCorr[,3,drop=FALSE],rowSums(results_clean_ctCorr[,4:6]),results_clean_ctCorr[,7:9])
colnames(results_clean_cted)[c(1,3)] = c("B-cells","CD4")

results_clean_epic_cted = cbind(rowSums(rclean_Ep2v2_ct[,1:2]),rclean_Ep2v2_ct[,3,drop=FALSE],rowSums(rclean_Ep2v2_ct[,4:6]),rclean_Ep2v2_ct[,7:9])
colnames(results_clean_epic_cted)[c(1,3)] = c("B-cells","CD4")

cor(c(as.matrix(flow_ed)),c(results_clean_dj_cted))
cor(c(as.matrix(flow_ed)),c(results_clean_cted))
cor(c(as.matrix(flow_ed)),c(results_clean_epic_cted))

sum((results_clean_cted-flow_ed)^2)
sum((results_clean_dj_cted-flow_ed)^2)
sum((results_clean_epic_cted-flow_ed)^2)

#-------------------------------------------------------------------#
# Revised Plots                                                     #
#-------------------------------------------------------------------#
fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_ICeDT%s_target_%s_ctCorr_ed.pdf",
              useWgts,wgtOpt,useTarget)
pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_dj_cted)[2]){
  plot(flow_ed[,i],results_clean_dj_cted[,i], main=colnames(results_clean_dj_cted)[i],xlab="Flow cytometry (%)", ylab="ICeD-T (%)", las=1)	
  abline(lm(results_clean_dj_cted[,i]~flow_ed[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow_ed[,i],results_clean_dj_cted[,i]),2)),paste("P=", round(cor.test(flow_ed[,i],results_clean_dj_cted[,i])$p.value,4))))
  i = i + 1
}

dev.off()

fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_CIBERSORT_web_ctCorr_ed.pdf",
              useWgts,wgtOpt,useTarget)
pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_cted)[2]){
  plot(flow_ed[,i],results_clean_cted[,i], main=colnames(results_clean_cted)[i],xlab="Flow cytometry (%)", ylab="CIBERSORT (%)", las=1)	
  abline(lm(results_clean_cted[,i]~flow_ed[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow_ed[,i],results_clean_cted[,i]),2)),paste("P=", round(cor.test(flow_ed[,i],results_clean_cted[,i])$p.value,4))))
  i = i + 1
}

dev.off()

fnm = sprintf("./CIBERSORT_useWgts_%s/Fig3a-plots_EPIC_ctCorr_ed.pdf",
              useWgts,wgtOpt,useTarget)
pdf(file = fnm, height=9, width=9)

par( mfrow = c( 3, 3 ) )

i = 1
while(i <= dim(results_clean_epic_cted)[2]){
  plot(flow_ed[,i],results_clean_epic_cted[,i], main=colnames(results_clean_epic_cted)[i],xlab="Flow cytometry (%)", ylab="ICeD-T (%)", las=1)	
  abline(lm(results_clean_epic_cted[,i]~flow_ed[,i]),col="red")
  legend("topleft", legend=c(paste("R=",round(cor(flow_ed[,i],results_clean_epic_cted[,i]),2)),paste("P=", round(cor.test(flow_ed[,i],results_clean_epic_cted[,i])$p.value,4))))
  i = i + 1
}

dev.off()

