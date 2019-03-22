#####################################################################
# EPIC TITRATION DATA: TESTING THEIR BEST                           #
#####################################################################
# PROGRAM NAME:                                                     #
#    EPIC_titrationData.R                                           #
# PROGRAMMER:                                                       #
#    Douglas Roy Wilson, Jr.                                        #
# DATE CREATED:                                                     #
#    02/28/2018                                                     #
# LAST EDIT:                                                        #
#    02/28/2018                                                     #
# VERSION:                                                          #
#    R-3.3.1                                                        #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#    Testing against EPIC's validation data. First Melanoma data    #
#####################################################################

#-------------------------------------------------------------------#
# LIBRARY and CODE LOAD                                             #
#-------------------------------------------------------------------#
setwd("D://DougData/Documents/Dissertation/Paper 3 - Immune Cell Deconv/")
setwd("./Sun_Requested_Analysis/DrSun_GrantFigures/")

library(EPIC)
library(alabama)

source("./programs/EPIC_Extract.R")

#-------------------------------------------------------------------#
# LOAD Melanoma Data                                                #
#-------------------------------------------------------------------#
load("./programs/EPIC-master/data/melanoma_data.rda")

EPICmdat = EPIC(bulk = melanoma_data$counts,reference = "TRef")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RESULTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Recast EPIC PREDICTION
EPIC.predr = EPICmdat$cellFractions[,c("Bcells","CD4_Tcells","CD8_Tcells","NKcells","otherCells")]
EPIC.predr[,5] = rowSums(EPICmdat$cellFractions[,c("CAFs","Endothelial","Macrophages","otherCells")])

EPIC.GT = melanoma_data$cellFractions.obs[,-c(6)]
EPIC.GT[,5] = rowSums(melanoma_data$cellFractions.obs[,c(5:6)])

### EPIC RECAST PREDICTION
#           Bcells   CD4_Tcells  CD8_Tcells    NKcells    otherCells
# LAU125  0.01016774 0.03041567 0.009302975 9.058697e-09  0.9501136
# LAU1255 0.04136029 0.05849291 0.129023216 4.706589e-07  0.7711231
# LAU1314 0.67562140 0.08198931 0.015668658 4.095899e-03  0.2226247
# LAU355  0.45493081 0.26813661 0.016429308 4.601164e-08  0.2605032

### GROUND TRUTH
#         Bcells   CD4_Tcells CD8_Tcells  NKcells  Cancer_cells  other_cells
# LAU125  0.1812     0.0082     0.0035    0.0050       0.6803      0.1218
# LAU1255 0.0579     0.0276     0.0376    0.0017       0.3756      0.4997
# LAU1314 0.4667     0.1815     0.0454    0.0025       0.0007      0.3031
# LAU355  0.3248     0.2315     0.0582    0.0017       0.0006      0.3832

### EPIC PREDICTION
#           Bcells       CAFs     CD4_Tcells  CD8_Tcells  Endothelial  Macrophages   NKcells    otherCells
# LAU125  0.01016774 2.359222e-04 0.03041567 0.009302975 2.554933e-02  0.01193535 9.058697e-09  0.9123930
# LAU1255 0.04136029 5.590262e-04 0.05849291 0.129023216 4.649865e-07  0.01961274 4.706589e-07  0.7509509
# LAU1314 0.67562140 1.163239e-08 0.08198931 0.015668658 2.799781e-07  0.00144306 4.095899e-03  0.2211814
# LAU355  0.45493081 7.457349e-09 0.26813661 0.016429308 4.095153e-09  0.00861752 4.601164e-08  0.2518857

#-------------------------------------------------------------------#
# HEADER.B                                                          #
#-------------------------------------------------------------------#
### Load Variance Data
TRefVar  = load("./data/TRef_Data/TRef_Var_weights.RData")
TRefPdat = load("./data/TRef_Data/TRef_purData.RData")

### Run EPIC Extract
EPIC_res_ICeDT = EPIC_Extract(bulk = melanoma_data$counts,reference = "TRef")

bulk = EPIC_res_ICeDT$bulk
refProfiles = EPIC_res_ICeDT$ref
refVar      = TRef_vprof_r[rownames(bulk),]

colnames(refVar) = colnames(refProfiles)

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 5 - LM22 FIT                               #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
setwd("./data/")

sig_matrix = './FlowCytometry/LM22.txt'
X = read.table(sig_matrix, header=TRUE, sep="\t", row.names=1, 
               check.names=FALSE)

X = data.matrix(X)

dim(X)
X[1:2,1:5]

colnames(X)
summary(X)

### Process Ref
ref  = X
cgenes = intersect(rownames(ref),rownames(melanoma_data$counts))

ref = ref[cgenes,]

### Process Bulk
bulk   = melanoma_data$counts
bulk   = bulk[cgenes,]

bulkstd = matrix(0,nrow=nrow(bulk),ncol=ncol(bulk))
rownames(bulkstd) = rownames(bulk)
colnames(bulkstd) = colnames(bulk)
for(i in 1:ncol(bulk)){
  bulkstd[,i] = (bulk[,i]-mean(bulk[,i]))/sd(bulk[,i])
}

### Put Bulk and Reference on Same playing field
mu_X = mean(ref)
sd_X = sd(as.vector(ref))

bulkREV = bulkstd*sd_X+mu_X

all(rownames(bulkREV)==rownames(ref))
dim(bulkREV)
dim(ref)

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 5 - Run Data                               #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
setwd("../")

source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_UpdateFunctions.R")

fitnw_nT = ASICeD_fit(Y=(bulkREV+0.000001),fixedCT_rho=rep(0,ncol(bulkREV)),useRho = FALSE,RefMat = ref,
                      maxIter = 500,maxIter_prop = 500,
                      maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulkREV));

source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_UpdateFunctions.R")

fitw0_nT = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=rep(0,ncol(bulk)),useRho = FALSE,RefMat = as.matrix(refProfiles+0.00001),useIQR = 0,
                      RefVar = refVar,varLog = TRUE,limitWgt = TRUE,limitQuant = c(0.20,0.80), maxIter = 500,maxIter_prop = 500,
                      maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

fitw2_nT = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=rep(0,ncol(bulk)),useRho = FALSE,RefMat = as.matrix(refProfiles+0.00001),useIQR = 2,
                      RefVar = refVar,varLog = TRUE,limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter = 500,maxIter_prop = 500,
                      maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

#-------------------------------------------------------------------#
# CONTROL for TUMOR                                                 #
#-------------------------------------------------------------------#
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_UpdateFunctions.R")

fitnw_T = ASICeD_fit(Y=(bulk+0.000001),fixedCT_rho=EPIC.GT[,5],useRho = TRUE,RefMat = as.matrix(refProfiles+0.000001),
                     maxIter = 500,maxIter_prop = 500,
                     maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_UpdateFunctions.R")

fitw0_T = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=EPIC.GT[,5],useRho = TRUE,RefMat = as.matrix(refProfiles+0.00001),useIQR = 0,
                     RefVar = refVar,varLog = TRUE,limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter = 500,maxIter_prop = 500,
                     maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

fitw2_T = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=EPIC.GT[,5],useRho = TRUE,RefMat = as.matrix(refProfiles+0.00001),useIQR = 2,
                     RefVar = refVar,varLog = TRUE,limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter = 500,maxIter_prop = 500,
                     maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

#-------------------------------------------------------------------#
# CONTROL for TUMOR (VERSION 2)                                     #
#-------------------------------------------------------------------#
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 4 - Supp Ref/P3_M2_UpdateFunctions.R")

fitnw_T2 = ASICeD_fit(Y=(bulk+0.000001),fixedCT_rho=melanoma_data$cellFractions.obs[,"Cancer_cells"],
                      useRho = TRUE,RefMat = as.matrix(refProfiles+0.000001),maxIter = 500,
                      maxIter_prop = 500,maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_initFit.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_FittingAlgorithm.R")
source("./programs/Full Model - Versions 4+5 with Diff EM Weights/Version 5 - Supp Ref/P3_M2_UpdateFunctions.R")

fitw0_T2 = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=melanoma_data$cellFractions.obs[,"Cancer_cells"],useRho = TRUE,
                      RefMat = as.matrix(refProfiles+0.00001),useIQR = 0,RefVar = refVar,varLog = TRUE,
                      limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter = 500,maxIter_prop = 500,
                      maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

fitw0_T2 = ASICeD_fit(Y=(bulk+0.00001),fixedCT_rho=melanoma_data$cellFractions.obs[,"Cancer_cells"],useRho = TRUE,
                      RefMat = as.matrix(refProfiles+0.00001),useIQR = 0,RefVar = refVar,varLog = TRUE,
                      limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter = 500,maxIter_prop = 500,
                      maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));


save(file="./data/EPIC_titrationData_fitsv2.RData",fitnw_T,fitnw_nT,fitw0_T,fitw0_nT,fitw2_T,fitw2_nT)
load(file="./data/EPIC_titrationData_fitsv2.RData")


fitTypes = c("nw_nT","nw_T","nw_T2","w0_nT","w0_T","w0_T2","w2_nT","w2_T","w2_T2")

for(j in fitTypes){
  eval(parse(text=sprintf("renorm4 = fit%s$IC_Abundance",j)))
  renorm4 = t(t(renorm4)/c(0.40,0.40,0.40,0.40,0.40,0.40,1.42,0.43))
  renorm4 = renorm4/rowSums(renorm4)
  
  renorm4.recast = renorm4[,c("Bcells","CD4_Tcells","CD8_Tcells","NKcells","Dummy_Tumor")]
  renorm4.recast[,5] = rowSums(renorm4[,c("Dummy_Tumor","CAFs","Endothelial",'Macrophages')])
  
  ### Raw Comparison
  if(j%in%c("nw_nT","w0_nT","w2_nT")){
    sum((EPIC.GT - EPIC.predr)^2)
    sum((EPIC.GT - renorm4.recast)^2)
    
    cor(c(EPIC.GT),c(EPIC.predr))
    cor(c(EPIC.GT),c(renorm4.recast))
  } else {
    sum((EPIC.GT[,-c(5)] - EPIC.predr[,-c(5)])^2)
    sum((EPIC.GT[,-c(5)] - renorm4.recast[,-c(5)])^2)
    
    cor(c(EPIC.GT[,-c(5)]),c(EPIC.predr[,-c(5)]))
    cor(c(EPIC.GT[,-c(5)]),c(renorm4.recast[,-c(5)]))
  }
  
  
  ### Relative Propn Immune Cells
  EPIC.predr.rel = EPIC.predr[,-c(5)]/rowSums(EPIC.predr[,-c(5)])
  EPIC.GT.rel    = EPIC.GT[,-c(5)]/rowSums(EPIC.GT[,-c(5)])
  renorm4.rel    = renorm4.recast[,-c(5)]/rowSums(renorm4.recast[,-c(5)])
  
  sum((EPIC.GT.rel-EPIC.predr.rel)^2)
  sum((EPIC.GT.rel-renorm4.rel)^2)
  
  cor(c(EPIC.GT.rel),c(EPIC.predr.rel))
  cor(c(EPIC.GT.rel),c(renorm4.rel))
}

# Extras:
EPIC.GT.rel
EPIC.predr.rel
renorm4.rel

rowSums((EPIC.GT.rel-renorm4.rel)^2)
sum((EPIC.GT.rel[-c(1),]-renorm4.rel[-c(1),])^2)
cor(c(EPIC.GT.rel[-c(1),]),c(renorm4.rel[-c(1),]))

rowSums((EPIC.GT.rel-EPIC.predr.rel)^2)
sum((EPIC.GT.rel[-c(1),]-EPIC.predr.rel[-c(1),])^2)
cor(c(EPIC.GT.rel[-c(1),]),c(EPIC.predr.rel[-c(1),]))


#-------------------------------------------------------------------#
# EXTRAS                                                            #
#-------------------------------------------------------------------#
renorm4 = fitnw_nT$IC_Abundance[,-c(1)]
renorm4 = cbind(rowSums(renorm4[,c(1:2)]),rowSums(renorm4[,c(5:7)]),
                         renorm4[,c(4)],rowSums(renorm4[,c(13:16)]),
                         rowSums(renorm4[,c(11:12)]))
colnames(renorm4) = c("Bcells","CD4_Tcells","CD8_Tcells","Macrophages","NKcells")

renorm4ct  = t(t(renorm4)/c(0.4,0.4,0.4,1.42,0.43))
renorm4ct  = renorm4ct/rowSums(renorm4ct)

renorm4rel = renorm4ct[,-c(4)]/rowSums(renorm4ct[,-c(4)])

rowSums((EPIC.GT.rel-CSORT_Full_LM22rel)^2)
sum((EPIC.GT.rel[-c(1),]-CSORT_Full_LM22rel[-c(1),])^2)
cor(c(EPIC.GT.rel[-c(1),]),c(CSORT_Full_LM22rel[-c(1),]))

### Relative Propn Immune Cells
EPIC.predr.rel = EPIC.predr[,-c(5)]/rowSums(EPIC.predr[,-c(5)])
EPIC.GT.rel    = EPIC.GT[,-c(5)]/rowSums(EPIC.GT[,-c(5)])
renorm4.rel    = renorm4.recast[,-c(5)]/rowSums(renorm4.recast[,-c(5)])

rowSums((EPIC.predr.rel-EPIC.GT.rel)^2)
rowSums((renorm4.rel-EPIC.GT.rel)^2)

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#              SECTION 2 - WRITE OUT CIBERSORT DATA                 #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
### Set directory
setwd("./data/")

### Save Counts
CSORT_Exp = melanoma_data$counts
CSORT_Out = data.frame(Gene_title=rownames(CSORT_Exp),CSORT_Exp)

write.table(file="Melanoma_Exp_Full.txt",x = CSORT_Out,quote = FALSE,row.names=FALSE,col.names = TRUE,sep="\t")

### Save TRef
RefProf = refProfiles
TREF_out = data.frame(Gene_title=rownames(RefProf),RefProf)

write.table(file="TRef_Output.txt",x=TREF_out,quote = FALSE,row.names=FALSE,col.names = TRUE,sep="\t")

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 3 - CIBERSORT Read IN                      #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
### Full, LM22 Ref
CSORT_Full_LM22 = read.table(file = "./data/CIBERSORT.Output_FULL_lm22.csv",
                             sep=",",header=TRUE,row.names = 1)

CSORT_Full_LM22r = cbind(rowSums(CSORT_Full_LM22[,c(1:2)]),rowSums(CSORT_Full_LM22[,c(5:7)]),
                         CSORT_Full_LM22[,c(4)],rowSums(CSORT_Full_LM22[,c(13:16)]),
                         rowSums(CSORT_Full_LM22[,c(11:12)]))
colnames(CSORT_Full_LM22r) = c("Bcells","CD4_Tcells","CD8_Tcells","Macrophages","NKcells")

CSORT_Full_LM22ct  = t(t(CSORT_Full_LM22r)/c(0.4,0.4,0.4,1.42,0.43))
CSORT_Full_LM22ct  = CSORT_Full_LM22ct/rowSums(CSORT_Full_LM22ct)

CSORT_Full_LM22rel = CSORT_Full_LM22ct[,-c(4)]/rowSums(CSORT_Full_LM22ct[,-c(4)])

rowSums((EPIC.GT.rel-CSORT_Full_LM22rel)^2)
sum((EPIC.GT.rel[-c(1),]-CSORT_Full_LM22rel[-c(1),])^2)
cor(c(EPIC.GT.rel[-c(1),]),c(CSORT_Full_LM22rel[-c(1),]))

### Rescale, TRef
CSORT_Resc_TRef = read.table(file = "./data/CIBERSORT.Output_Rescale_TRef.csv",
                             sep=",",header=TRUE,row.names=1)

CSORT_Resc_TRefr = t(t(CSORT_Resc_TRef[,-c(8:10)])/c(0.40,0.40,0.40,0.40,0.40,1.42,0.43))
CSORT_Resc_TRefr = CSORT_Resc_TRefr/rowSums(CSORT_Resc_TRefr)

CSORT_Resc_TRefrel = CSORT_Resc_TRefr[,-c(2,5,6)]/rowSums(CSORT_Resc_TRefr[,-c(2,5,6)])

rowSums((EPIC.GT.rel-CSORT_Resc_TRefrel)^2)
sum((EPIC.GT.rel[-c(1),]-CSORT_Resc_TRefrel[-c(1),])^2)
cor(c(EPIC.GT.rel[-c(1),]),c(CSORT_Resc_TRefrel[-c(1),]))

#-------------------------------------------------------------------#
# Plots                                                             #
#-------------------------------------------------------------------#
## Run the above for fit type ("w0_T")
colnames(renorm4.recast)[5] = "otherCells"
ICeDT_Dat <- melt(renorm4.recast)
colnames(ICeDT_Dat) = c("Subject","CellType","ICeD_T")

EPIC_Dat <- melt(EPIC.predr)
colnames(EPIC_Dat) = c("Subject","CellType","EPIC")

colnames(EPIC.GT)[5] = "otherCells"
GT_Dat <- melt(EPIC.GT)
colnames(GT_Dat) = c("Subject","CellType","Prop")

ICeDT_Plot_Dat <- merge(x = ICeDT_Dat,y = GT_Dat,by=c("Subject","CellType"))
ICeDT_Plot_Dat$CellGroup = ifelse(test = ICeDT_Plot_Dat$CellType=="otherCells",yes = "Non-Immune",no="Immune")
EPIC_Plot_Dat <-  merge(x = EPIC_Dat,y = GT_Dat,by=c("Subject","CellType"))
EPIC_Plot_Dat$CellGroup = ifelse(test = EPIC_Plot_Dat$CellType=="otherCells",yes = "Non-Immune",no="Immune")


EPIC_Plot <- ggplot(data=EPIC_Plot_Dat,aes(x=Prop,y=EPIC,colour=CellGroup))+
  geom_point()

ICeDT_Plot <-  ggplot(data=ICeDT_Plot_Dat,aes(x=Prop,y=ICeD_T,colour=CellGroup))+
  geom_point()

pdf(file=sprintf("./programs/EPIC_TitrationData/ICeDT_wwgt_wTP_plot.pdf"),height = 4,width = 4)
print(ICeDT_Plot)
dev.off()

pdf(file=sprintf("./programs/EPIC_TitrationData/EPIC_plot.pdf"),height = 4,width = 4)
print(EPIC_Plot)
dev.off()