### LOAD Edited EPIC functions and Reference Data
source("EPIC_Extract.R")
TRefPdat = load("../data/TRef_Data/TRef_purData.RData")

### Use the EPIC gene set but rescale pure data for reweighting
load("../programs/EPIC-master/data/TRef.rda")

Egenes = TRef$sigGenes

## Let MixDat_tot be your mixture expression matrix across all subjects
commonGenes = intersect(rownames(Bdat_r),rownames(MixDat_tot))
Egenes = Egenes[which(Egenes%in%commonGenes)]

## Rescale the reference data using EPIC technique (renormalize TPM at common
## Genes to 10e6)

## XXXdat_r: Contains an nG by nP structure of TPM counts used by Racle et al to 
##           typify various immune populations. 
##   - B   : B cells
##   - CAF : Cancer Associated Fibroblasts
##   - CD4 : CD4 + T-cells
##   - CD8 : CD8 + T-cells
##   - E   : Endothelials
##   - M   : Macrophages
##   - NK  : Natural Killers

Bscl   = scaleCounts(counts = Bdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
CAFscl = scaleCounts(counts = CAFdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
CD4scl = scaleCounts(counts = CD4dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
CD8scl = scaleCounts(counts = CD8dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
Escl   = scaleCounts(counts = Edat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
Mscl   = scaleCounts(counts = Mdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
NKscl  = scaleCounts(counts = NKdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts

Brvar    = apply(X = log(Bscl+0.00001),MARGIN = 1,FUN = var)
CAFrvar  = apply(X = log(CAFscl+0.00001),MARGIN = 1,FUN = var)
CD4rvar  = apply(X = log(CD4scl+0.00001),MARGIN = 1,FUN = var)
CD8rvar  = apply(X = log(CD8scl+0.00001),MARGIN = 1,FUN = var)
Ervar    = apply(X = log(Escl+0.00001),MARGIN = 1,FUN = var)
Mrvar    = apply(X = log(Mscl+0.00001),MARGIN = 1,FUN = var)
NKrvar   = apply(X = log(NKscl+0.00001),MARGIN = 1,FUN = var)

EPIC_res_ICeDT = EPIC_Extract(bulk = MixDat_tot,scaleExprs = TRUE,reference = "TRef")

g2useFin = intersect(names(Brvar),commonGenes)
g2useFin = intersect(g2useFin,Egenes)

bulk = EPIC_res_ICeDT$bulk[g2useFin,]
refProfiles = EPIC_res_ICeDT$ref[g2useFin,]
refVar      = cbind(Brvar,CAFrvar,CD4rvar,CD8rvar,Ervar,Mrvar,NKrvar)
refVar      = refVar[g2useFin,]
colnames(refVar) = c("Bcells","CAFs","CD4_Tcells","CD8_Tcells","Endothelial","Macrophages","NKcells")


#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 2 - FIT with ICED-T                        #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
### ICeD-T with no weights and with a supplied reference
fitnw = ICeDT_fit_noWgt_suppRef(Y=bulk,fixedCT_rho=rep(0,ncol(bulk)),useRho = FALSE,
                              RefMat = as.matrix(refProfiles+0.000001),
                              maxIter_prop = 500,
                              maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));

### ICeD-T with weights and a supplied reference:
fitww = ICeDT_fit_Wgt_suppRef(Y=bulk,fixedCT_rho=rep(0,ncol(bulk)),useRho = FALSE,RefMat = as.matrix(refProfiles+0.00001),useIQR = 0,
                              RefVar = refVar,varLog = TRUE,limitWgt = TRUE,limitQuant = c(0.15,0.85), maxIter_prop = 500,
                              maxIter_PP = 250,RhoConv_CO = 1e-3,Subj_CO = ncol(bulk));
