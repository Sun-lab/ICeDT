#####################################################################
# STEP 3 - Reprocess Data                                           #
#####################################################################
# PROGRAM NAME:                                                     #
#   running_ICeDT_on_Hugo_data.R                                    #
# PROGRAMMER:                                                       #
#   Douglas Roy Wilson, Jr.                                         #
#   Chong Jin                                                       #
# DATE CREATED:                                                     #
#   02/19/2018                                                      #
# LAST EDIT:                                                        #
#   03/21/2018                                                      #
# VERSION:                                                          #
#   R-3.5.1                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   Reads in data prepared by Chong Jin and Wei Sun, runs the       #
#   ICeD-T deconvolution procedure, and saves output for further    #
#   processing.                                                     #
#                                                                   #
#   Outputs from ICeD-T:                                            #
#     fitnw -- no weight                                            #
#     fitw0 -- use weight                                           #
#####################################################################

#-------------------------------------------------------------------#
#                SECTION 0 - RUN PARAMETERS                         #
#-------------------------------------------------------------------#
#                                                                   #
# Geneset: choose one of {"Revised", "Original"}                    #
#     - "Revised"  = EPIC Genes, LM22 Genes, MCP-Counter genes      #
#     - "Original" = EPIC utilized genes                            #
#-------------------------------------------------------------------#

Geneset = "Original"

#-------------------------------------------------------------------#
#                SECTION 1 - LIBRARIES and CODE                     #
#-------------------------------------------------------------------#

library(nnls)
library(quantreg)
library(hqreg)
library(gplots)
library(org.Hs.eg.db)
library(alabama)
library(EPIC)
library(clinfun)
library(ICeDT)

#-------------------------------------------------------------------#
#                SECTION 2 - READ-IN PROVIDED DATA                  #
#-------------------------------------------------------------------#

# -----------------------------------------------------------------
# OBTAIN LM22 Genes
# -----------------------------------------------------------------
setwd("..")
sig_matrix = './data/FlowCytometry/LM22.txt'
X = read.table(sig_matrix, header=TRUE, sep="\t", row.names=1, 
               check.names=FALSE)

X = data.matrix(X)

dim(X)
X[1:2,1:5]

colnames(X)
summary(X)

# ------------------------------------------------------------
# LOAD gene Annotation Data
# ------------------------------------------------------------
geneInfo = read.table("./data/gencode.v27.genes.txt", sep="\t", as.is=TRUE, 
                      header=TRUE,quote="")

geneInfo$geneId2 = unlist(lapply(strsplit(x = geneInfo$geneId,split = "\\."),
                                 '[[',1))

dim(geneInfo)
geneInfo[1:2,]
length(unique(geneInfo$geneId))

# ------------------------------------------------------------
# EPIC Genes
# ------------------------------------------------------------

FList2 = c("BANK1","CD79A", "CD79B", "FCER2", "FCRL2", "FCRL5", "MS4A1", 
           "PAX5", "POU2AF1", "STAP1", "TCL1A", "ADAM33", "CLDN11", "COL1A1", 
           "COL3A1", "COL14A1", "CRISPLD2", "CXCL14", "DPT", "F3", "FBLN1", 
           "ISLR", "LUM", "MEG3", "MFAP5", "PRELP", "PTGIS", "SFRP2", "SFRP4", 
           "SYNPO2", "TMEM119", "ANKRD55", "DGKA", "FOXP3", "GCNT4", "IL2RA", 
           "MDS2", "RCAN3", "TBC1D4", "TRAT1", "CD8B", "HAUS3", "JAKMIP1",
           "NAA16", "TSPYL1", "CDH5", "CLDN5", "CLEC14A", "CXorf36", "ECSCR", 
           "F2RL3", "FLT1", "FLT4", "GPR4", "GPR182", "KDR", "MMRN1", "MMRN2", 
           "MYCT1", "PTPRB", "RHOJ", "SLCO2A1", "SOX18", "STAB2", "VWF",
           "APOC1", "C1QC", "CD14", "CD163", "CD300C", "CD300E", "CSF1R", 
           "F13A1", "FPR3", "HAMP", "IL1B", "LILRB4", "MS4A6A", "MSR1", 
           "SIGLEC1", "VSIG4", "CD33", "CD300C", "CD300E", "CECR1", "CLEC6A", 
           "CPVL", "EGR2", "EREG", "MS4A6A", "NAGA", "SLC37A2",
           "CEACAM3", "CNTNAP3", "CXCR1", "CYP4F3", "FFAR2", "HIST1H2BC", 
           "HIST1H3D", "KY", "MMP25", "PGLYRP1", "SLC12A1", "TAS2R40",
           "CD160", "CLIC3", "FGFBP2", "GNLY", "GNPTAB", "KLRF1", "NCR1", 
           "NMUR1", "S1PR5", "SH2D1B", "BCL11B", "CD5", "CD28", "IL7R", 
           "ITK", "THEMIS", "UBASH3A")

FList_tot = FList2

# ------------------------------------------------------------
# load gene expression data parsed from Hugo et al. 2016
# ------------------------------------------------------------

load("./data/FPKM_and_counts_filtered.RData")
ls()

dim(FPKM)
FPKM[1:2,1:5]

dim(counts)
counts[1:2,1:5]
table(rownames(counts) == geneInfo$geneId)
table(rownames(X) %in% geneInfo$hgnc_symbol)

### Match FLIst to ensembl_gene_id
mat_FList = match(FList_tot,geneInfo$hgnc_symbol)
table(FList_tot==geneInfo$hgnc_symbol[mat_FList])
FList_ens = geneInfo$ensembl_gene_id[mat_FList]
FList_ens = FList_ens[!is.na(FList_ens)]

# --------------------------------------------------------------------
# try to identify those genes that we do not have matched gene symbols
# --------------------------------------------------------------------

w2check = which(! rownames(X) %in% geneInfo$hgnc_symbol)
w2check
rownames(X)[w2check]
keys = rownames(X)[w2check]

columns(org.Hs.eg.db)
einfo = select(org.Hs.eg.db, keys=keys, columns = c("SYMBOL","ENSEMBL"), 
               keytype="ALIAS")

dim(einfo)
einfo

mat1 = match(einfo$SYMBOL, geneInfo$hgnc_symbol)
mat2 = match(einfo$ENSEMBL, geneInfo$ensembl_gene_id)

table(mat1 == mat2, useNA="ifany")
table(is.na(mat1), is.na(mat2))

einfo[which(is.na(mat2)),]

mat3 = na.omit(mat2)
cbind(einfo[which(!is.na(mat2)),], geneInfo[mat3,1:7])

einfo[which(is.na(mat2) & (! is.na(mat1))),]
einfo[which(! is.na(mat2) & (is.na(mat1))),]

geneSym = keys
newSym  = einfo$SYMBOL[match(keys, einfo$ALIAS)]
cbind(geneSym, newSym)

w2update = which(!is.na(newSym))
geneSym[w2update] = newSym[w2update]
cbind(geneSym, newSym)

geneInfo$hgnc_symbol[which(geneInfo$hgnc_symbol == "KIRREL1")] = "KIRREL"
geneInfo$hgnc_symbol[which(geneInfo$hgnc_symbol == "MACF1")] = "KIAA0754"

rownames(X)[w2check] = geneSym
table(rownames(X) %in% geneInfo$hgnc_symbol)

rownames(X)[!rownames(X) %in% geneInfo$hgnc_symbol]

#-------------------------------------------------------------------#
#                SECTION 3 - MERGING DATASETS                       #
#-------------------------------------------------------------------#

### Remove .XX from ENSG00000000003.XX to make consisent with Linsley
counts_rn = unlist(lapply(strsplit(x = rownames(counts),split="\\."),
                          '[[',1))
counts_rn[1:5]
#rownames(counts)

### Handling Duplicates from Paralogous Y
id2chk = which(duplicated(counts_rn))
rownames(counts)[id2chk]

rows2chk = rownames(counts)[id2chk]
max.r = apply(X = counts[rows2chk,],1,max)

counts2 = counts[-which(rownames(counts)%in%rows2chk),]
counts  = counts2

counts_rn = unlist(lapply(strsplit(x = rownames(counts),split="\\."),
                          '[[',1))
rownames(counts) = counts_rn

geneInfo2 = geneInfo[-c(which(geneInfo$geneId%in%rows2chk)),]
geneInfo  = geneInfo2

### Normalize Data TPM
countsS = counts

Exp2use = countsS

r25  = apply(Exp2use, 1, quantile, probs = 0.25)
w2kp = which((r25 > 10)|(rownames(Exp2use)%in%FList_ens))

Exp2use = Exp2use[w2kp,]
tot = colSums(Exp2use)
s75 = apply(Exp2use, 2, quantile, prob = 0.75)

cor(tot, s75)

GeneLengths_Mat = readRDS("./data/gene_lengths_v27.rds")
dim(GeneLengths_Mat)
GeneLengths_Mat[1:2,]

w2rm = which(GeneLengths_Mat$Gencode_ID%in%rows2chk)
GeneLengths_Mat = GeneLengths_Mat[-c(w2rm),]

mat1 = match(rownames(Exp2use), 
             unlist(lapply(X = strsplit(x = GeneLengths_Mat$Gencode_ID,split ="\\."),'[[',1)))
L = GeneLengths_Mat[mat1, "Exonic"]
TPM = (Exp2use+1)/L
TPM = (1e6)*t(t(TPM)/colSums(TPM))

geneInfo = geneInfo[which(geneInfo$geneId2%in%rownames(TPM)),]

all(geneInfo$geneId2==rownames(TPM))
table(is.na(geneInfo$hgnc_symbol))
length(unique(geneInfo$hgnc_symbol))

#-------------------------------------------------------------------#
#                SECTION 4 - Reducing Unlabeled Genes               #
#-------------------------------------------------------------------#
t1 = table(geneInfo$hgnc_symbol)
sort(t1, decreasing=TRUE)[1:5]

ww1 = which(geneInfo$hgnc_symbol == "COG8")
ww1
apply(TPM[ww1,], 1, summary)

ww1 = which(geneInfo$hgnc_symbol == "POLR2J4")
ww1
apply(TPM[ww1,], 1, summary)

w2rm = union(which(geneInfo$hgnc_symbol == ""), c(13745,13904))

TPM = TPM[-c(w2rm),]
geneInfo = geneInfo[-c(w2rm),]

dim(TPM)
dim(geneInfo)
all(rownames(TPM)==geneInfo$ensembl_gene_id)

#-------------------------------------------------------------------#
#                SECTION 5 - GeneSet and Weights                    #
#-------------------------------------------------------------------#
### Run Epic and use EPIC_Extract.R
MixDat_tot = TPM
rownames(MixDat_tot) = geneInfo$hgnc_symbol

EPIC_Bref = EPIC(bulk = MixDat_tot,scaleExprs = TRUE,reference = "BRef")
EPIC_Tref = EPIC(bulk = MixDat_tot,scaleExprs = TRUE,reference = "TRef")

### Prepare data for ICeD-T run
source("./programs/EPIC_Extract.R")
# Expression and variance summarized by cell type
# "TRef_leprof_r" "TRef_eprof_r"  "TRef_vprof_r" 
TRefVar  = load("./data/TRef_Data/TRef_Var_weights.RData")
# Raw expression by cell -- one object per cell type
# "Bdat_r"   "CAFdat_r" "CD4dat_r" "CD8dat_r" "Edat_r"   "Mdat_r"   "NKdat_r" 
TRefPdat = load("./data/TRef_Data/TRef_purData.RData")

if(Geneset=="Original"){
  ### Use the EPIC gene set but rescale pure data for reweighting
  load("./programs/EPIC-master/data/TRef.rda")
  
  Egenes = TRef$sigGenes
  
  commonGenes = intersect(rownames(Bdat_r),rownames(MixDat_tot))
  Egenes = Egenes[which(Egenes%in%commonGenes)]
  
  Bscl   = scaleCounts(counts = Bdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CAFscl = scaleCounts(counts = CAFdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CD4scl = scaleCounts(counts = CD4dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CD8scl = scaleCounts(counts = CD8dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  Escl   = scaleCounts(counts = Edat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  Mscl   = scaleCounts(counts = Mdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  NKscl  = scaleCounts(counts = NKdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  
  Brvar    = apply(X = log(Bscl+1e-5),MARGIN = 1,FUN = var)
  CAFrvar  = apply(X = log(CAFscl+1e-5),MARGIN = 1,FUN = var)
  CD4rvar  = apply(X = log(CD4scl+1e-5),MARGIN = 1,FUN = var)
  CD8rvar  = apply(X = log(CD8scl+1e-5),MARGIN = 1,FUN = var)
  Ervar    = apply(X = log(Escl+1e-5),MARGIN = 1,FUN = var)
  Mrvar    = apply(X = log(Mscl+1e-5),MARGIN = 1,FUN = var)
  NKrvar   = apply(X = log(NKscl+1e-5),MARGIN = 1,FUN = var)
  
  EPIC_res_ICeDT = EPIC_Extract(bulk = MixDat_tot,scaleExprs = TRUE,reference = "TRef")
  
  g2useFin = intersect(names(Brvar),commonGenes)
  g2useFin = intersect(g2useFin,Egenes)
  
  bulk = EPIC_res_ICeDT$bulk[g2useFin,]
  refProfiles = EPIC_res_ICeDT$ref[g2useFin,]
  refVar      = cbind(Brvar,CAFrvar,CD4rvar,CD8rvar,Ervar,Mrvar,NKrvar)
  refVar      = refVar[g2useFin,]
  colnames(refVar) = c("Bcells","CAFs","CD4_Tcells","CD8_Tcells","Endothelial","Macrophages","NKcells")
  
} else if(Geneset=="Revised"){
  load("./programs/EPIC-master/data/TRef.rda")
  
  Egenes = union(FList_tot,rownames(X))
  
  commonGenes = intersect(rownames(Bdat_r),rownames(MixDat_tot))
  Egenes = Egenes[which(Egenes%in%commonGenes)]
  
  MIXscl = scaleCounts(counts = MixDat_tot,sigGenes = Egenes,renormGenes = commonGenes)$counts
  Bscl   = scaleCounts(counts = Bdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CAFscl = scaleCounts(counts = CAFdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CD4scl = scaleCounts(counts = CD4dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  CD8scl = scaleCounts(counts = CD8dat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  Escl   = scaleCounts(counts = Edat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  Mscl   = scaleCounts(counts = Mdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  NKscl  = scaleCounts(counts = NKdat_r,sigGenes = Egenes,renormGenes = commonGenes)$counts
  
  Bexp     = apply(X = Bscl,MARGIN=1,FUN=mean)
  CAFexp   = apply(X = CAFscl,MARGIN=1,FUN=mean)
  CD4exp   = apply(X = CD4scl,MARGIN=1,FUN=mean)
  CD8exp   = apply(X = CD8scl,MARGIN=1,FUN=mean)
  Eexp     = apply(X = Escl,MARGIN=1,FUN=mean)
  Mexp     = apply(X = Mscl,MARGIN=1,FUN=mean)
  NKexp    = apply(X = NKscl,MARGIN=1,FUN=mean)
  
  Brvar    = apply(X = log(Bscl+1e-5),MARGIN = 1,FUN = var)
  CAFrvar  = apply(X = log(CAFscl+1e-5),MARGIN = 1,FUN = var)
  CD4rvar  = apply(X = log(CD4scl+1e-5),MARGIN = 1,FUN = var)
  CD8rvar  = apply(X = log(CD8scl+1e-5),MARGIN = 1,FUN = var)
  Ervar    = apply(X = log(Escl+1e-5),MARGIN = 1,FUN = var)
  Mrvar    = apply(X = log(Mscl+1e-5),MARGIN = 1,FUN = var)
  NKrvar   = apply(X = log(NKscl+1e-5),MARGIN = 1,FUN = var)
  
  bulk = MIXscl[Egenes,]
  
  refProfiles = cbind(Bexp,CAFexp,CD4exp,CD8exp,Eexp,Mexp,NKexp)
  refProfiles = refProfiles[Egenes,]
  colnames(refProfiles) = c("Bcells","CAFs","CD4_Tcells","CD8_Tcells",
                            "Endothelial","Macrophages","NKcells")
  
  refVar      = cbind(Brvar,CAFrvar,CD4rvar,CD8rvar,Ervar,Mrvar,NKrvar)
  refVar      = refVar[Egenes,]
  colnames(refVar) = c("Bcells","CAFs","CD4_Tcells","CD8_Tcells",
                       "Endothelial","Macrophages","NKcells")
  
}

#-------------------------------------------------------------------#
#                SECTION 6 - Run Data                               #
#-------------------------------------------------------------------#

# Using ICeDT package
fitnw = ICeDT::ICeDT(Y=bulk, Z=refProfiles, tumorPurity=NULL, refVar=NULL,
                     rhoInit=NULL, maxIter_prop = 500, maxIter_PP = 250, 
                     rhoConverge = 1e-3)

fitw0 = ICeDT::ICeDT(Y=bulk, Z=refProfiles, tumorPurity=NULL, refVar=refVar,
                     rhoInit=NULL, maxIter_prop = 500, maxIter_PP = 250, 
                     rhoConverge = 1e-3)

save(fitnw,fitw0,file=sprintf("./data/ICeDT_ExpandedFits_GeneSet%s.RData",Geneset))
save(EPIC_Bref, EPIC_Tref, file="./data/EPIC_Fits.RData")

load(sprintf("./data/ICeDT_ExpandedFits_GeneSet%s.RData",Geneset))

#-------------------------------------------------------------------#
#                SECTION 7 - PROCESS OUTPUT                         #
#-------------------------------------------------------------------#
pinfo = read.table("./data/patient_info.txt", sep="\t", as.is=TRUE, 
                   header=TRUE)
dim(pinfo)
pinfo[1:2,]

pinfo[,c("irRECIST", "SRA.tumor.RNA")]

pinfo$SRA.tumor.RNA
table(rownames(EPIC_Tref$cellFractions) %in% pinfo$SRA.tumor.RNA)

pinfo = pinfo[match(rownames(EPIC_Tref$cellFractions), pinfo$SRA.tumor.RNA),]
table(rownames(EPIC_Tref$cellFractions) == pinfo$SRA.tumor.RNA)
table(pinfo$irRECIST)

#-------------------------------------------------------------------#
#                SECTION 8 - PLOTTING RESULTS                       #
#-------------------------------------------------------------------#
figures_dir = sprintf("./figures/GeneSet%s/",Geneset)
dir.create(figures_dir, recursive = TRUE)

for(j in c("nw","w0")){
  eval(parse(text=sprintf("renorm4 = t(fit%s$rho)",j)))
  # normalize using cell size factors
  renorm4 = t(t(renorm4)/c(0.40,0.40,0.40,0.40,0.40,0.40,1.42,0.43))
  renorm4 = renorm4/rowSums(renorm4)
  
  renorm4 = renorm4[,-c(1)]
  renorm4 = renorm4/rowSums(renorm4)
  
  if(j=="nw"){renorm_ICeDT_nw = renorm4}
  if(j=="w0"){renorm_ICeDT_w0 = renorm4}
  
  toplot = renorm4
  toplot[which(toplot > 50)] = 50
  
  palette.gr.marray = colorRampPalette(c("blue", "white", "red"))
  subtype_cols = rep("darkgreen", nrow(pinfo))
  subtype_cols[which(pinfo$irRECIST == "Progressive Disease")] = "orange"
  table(pinfo$irRECIST, subtype_cols)
  
  pdf(file.path(figures_dir, sprintf("cell_types_ICeDT_%s.pdf",j)), width=8, height=8)
  
  heatmap.2(as.matrix(toplot), trace = "none", col = palette.gr.marray, 
            Rowv = TRUE, Colv = TRUE, dendrogram = "both", key = TRUE, 
            sepwidth = c(0.01, 0.01), offsetCol=0.1, srtCol=45, 
            key.title=NA, key.ylab=NA, margin=c(10,8), 
            RowSideColors = subtype_cols)
  
  dev.off()
  
  pdf(file.path(figures_dir, sprintf("cell_types_ICeDT_%s_boxplot.pdf",j)), width=8, height=6)
  par(mar=c(12,4,1,1), bty="n")
  boxplot(renorm4, las=2)
  dev.off()
  
  pdf(file.path(figures_dir, sprintf("ICeDT_%s_CD8_vs_response.pdf",j)), width=4, height=4)
  par(mar=c(10,4,1,1), bty="n", las=2)
  boxplot(renorm4[,4] ~ pinfo$irRECIST, ylab="CD8+ T cell proportion", outline=FALSE, 
          ylim=c(0, max(renorm4[,4])))
  set.seed(1234)
  stripchart(renorm4[,4] ~ pinfo$irRECIST, method="jitter", jitter = 0.2, 
             vertical=TRUE, add=TRUE, pch=20, col="darkred")
  dev.off()
}

### Tumor References
renorm4 = EPIC_Tref$cellFractions[,-c(8)]
renorm4 = renorm4/rowSums(renorm4)

renorm_EPIC = renorm4

toplot = renorm4
toplot[which(toplot > 50)] = 50

palette.gr.marray = colorRampPalette(c("blue", "white", "red"))
subtype_cols = rep("darkgreen", nrow(pinfo))
subtype_cols[which(pinfo$irRECIST == "Progressive Disease")] = "orange"
table(pinfo$irRECIST, subtype_cols)

pdf(file.path(figures_dir, "cell_types_EPIC_Tref.pdf"), width=8, height=8)

heatmap.2(as.matrix(toplot), trace = "none", col = palette.gr.marray, 
          Rowv = TRUE, Colv = TRUE, dendrogram = "both", key = TRUE, 
          sepwidth = c(0.01, 0.01), offsetCol=0.1, srtCol=45, 
          key.title=NA, key.ylab=NA, margin=c(10,8), 
          RowSideColors = subtype_cols)

dev.off()

pdf(file.path(figures_dir, "cell_types_EPIC_TREF_boxplot.pdf"), width=8, height=6)
par(mar=c(12,4,1,1), bty="n")
boxplot(renorm4, las=2)
dev.off()

pdf(file.path(figures_dir, "EPICTRef_CD8_vs_response.pdf"), width=4, height=4)
par(mar=c(10,4,1,1), bty="n", las=2)
boxplot(renorm4[,4] ~ pinfo$irRECIST, ylab="CD8+ T cell proportion", outline=FALSE, 
        ylim=c(0, max(renorm4[,4])))
set.seed(1234)
stripchart(renorm4[,4] ~ pinfo$irRECIST, method="jitter", jitter = 0.2, 
           vertical=TRUE, add=TRUE, pch=20, col="darkred")
dev.off()

#-------------------------------------------------------------------#
# CIBERSORT Read In                                                 #
#-------------------------------------------------------------------#
### Full, LM22 Ref
CSORT_Full_LM22 = read.table(file = "./data/CIBERSORT.Output_SunMelExt_LM22NQNorm.csv",
                             sep=",",header=TRUE,row.names = 1)

CSORT_CTSize = rep(0.4,22)

BLoc    = c(1,2)
PlasLoc = c(3)
TLoc    = c(4,5,6,7,8,9,10)
NKLoc   = c(11,12)
MonoLoc = c(13,14,15,16,17,18)
GranLoc = c(19,20,21,22)

CSORT_CTSize[NKLoc]   = 0.43
CSORT_CTSize[MonoLoc] = 1.42
CSORT_CTSize[GranLoc] = 0.15

renorm_CSORT = CSORT_Full_LM22[,1:22]
renorm_CSORT = t(t(renorm_CSORT)/CSORT_CTSize)
renorm_CSORT = renorm_CSORT/rowSums(renorm_CSORT)

par(mar=c(10,4,1,1), bty="n", las=2)

pdf("./figures/CIBERSORT_CD8_vs_response.pdf", width=4, height=4)
par(mar=c(10,4,1,1), bty="n", las=2)
boxplot(renorm_CSORT[,4] ~ pinfo$irRECIST, ylab="CD8+ T cell proportion", outline=FALSE, 
        ylim=c(0, max(renorm_CSORT[,4])))
set.seed(1234)
stripchart(renorm_CSORT[,4] ~ pinfo$irRECIST, method="jitter", jitter = 0.2, 
           vertical=TRUE, add=TRUE, pch=20, col="darkred")
dev.off()

### Using the TRef Data
CSORT_Full_TRef = read.table(file = "./data/CIBERSORT.Output_SunMelExt_TRefNQNorm.csv",
                             sep=",",header=TRUE,row.names = 1)

renorm_CSORT_TRef = CSORT_Full_TRef[,1:7]

renorm_CSORT_TRef = t(t(renorm_CSORT_TRef)/c(0.40,0.40,0.40,0.40,0.40,1.42,0.43))
renorm_CSORT_TRef = renorm_CSORT_TRef/rowSums(renorm_CSORT_TRef)

pdf("./figures/CIBERSORT_TRef_CD8_vs_response.pdf", width=4, height=4)
par(mar=c(10,4,1,1), bty="n", las=2)
boxplot(renorm_CSORT_TRef[,4] ~ pinfo$irRECIST, ylab="CD8+ T cell proportion", outline=FALSE, 
        ylim=c(0, max(renorm_CSORT_TRef[,4])))
set.seed(1234)
stripchart(renorm_CSORT_TRef[,4] ~ pinfo$irRECIST, method="jitter", jitter = 0.2, 
           vertical=TRUE, add=TRUE, pch=20, col="darkred")
dev.off()


#-------------------------------------------------------------------#
# JONCKHEERE-TERPSTRA test                                          #
#-------------------------------------------------------------------#

## Set category labels in numerical order (Complete = 1 , Partial = 2 , Progressive = 3)
disCat_num = ifelse(test = pinfo$irRECIST=="Complete Response",yes = 1,
                    no = ifelse(test = pinfo$irRECIST=="Partial Response",yes = 2,
                                no = ifelse(test = pinfo$irRECIST=="Progressive Disease",yes = 3,no=NA)))

table(disCat_num,pinfo$irRECIST)

### 
jonckheere.test(x=renorm_ICeDT_nw[,4],g = disCat_num,alternative = "decreasing")
jonckheere.test(x=renorm_ICeDT_w0[,4],g = disCat_num,alternative = "decreasing")
jonckheere.test(x=renorm_EPIC[,4],g = disCat_num,alternative = "decreasing")
jonckheere.test(x=renorm_CSORT[,4],g = disCat_num,alternative = "decreasing")

#-------------------------------------------------------------------#
#                SECTION 9 - CHECK ABBERANT PROBABILITY             #
#-------------------------------------------------------------------#

p0 = fitnw$cProb
dim(p0)

p1 = fitw0$cProb
dim(p1)


p0[1:2,1:5]
p1[1:2,1:5]

p0 = data.matrix(p0)
p1 = data.matrix(p1)

q90 <- function(v){
  qs = quantile(v, probs=c(0.10, 0.90))
  qs[2] - qs[1]
}

pdf(sprintf("./figures/probConsistent_GeneSet%s.pdf",Geneset), 
    width=9, height=3)
par(mar=c(5,4,1,1), bty="n", mfrow=c(1,3), cex=0.8)
plot(pmax(density(c(p0))$y, density(c(p1))$y), main="", xlim=c(0,1),
     xlab="probability consistent", ylab="density", type="n")
lines(density(c(p0)), lty=1, col="black")
lines(density(c(p1)), lty=2, col="red")
legend("topright", c("no weight", "w/ weight"), lty=c(1,2), 
       col=c("black", "red"), bty="n")

plot(apply(p0, 1, median), apply(p0, 1, q90), 
     xlab="median prob. consistent", ylab="90 percentile - 10 percentile")
plot(apply(p1, 1, median), apply(p1, 1, q90), 
     xlab="median prob. consistent", ylab="90 percentile - 10 percentile")

dev.off()

# Scatterplot of predicted vs. observed gene expression
dim(fitw0$rho[-1,])
dim(refProfiles)
dim(bulk)
predicted_bulk_nw = refProfiles %*% fitnw$rho[-1,]
predicted_bulk_w0 = refProfiles %*% fitw0$rho[-1,]
p0_cutoffs = quantile(p0, c(0.333,0.666))
p1_cutoffs = quantile(p1, c(0.333,0.666))

cat(sprintf("Consistent probability cutoffs for model w/o weight: %.3f, %.3f \n", 
            p0_cutoffs[1], p0_cutoffs[2]))
cat(sprintf("Consistent probability cutoffs for model w/ weight: %.3f, %.3f \n", 
            p1_cutoffs[1], p1_cutoffs[2]))

# Geneset=="Revised":
# Consistent probability cutoffs for model w/o weight: 0.558, 0.680 
# Consistent probability cutoffs for model w/ weight: 0.946, 0.983
# Geneset=="Original":
# Consistent probability cutoffs for model w/o weight: 0.509, 0.549 
# Consistent probability cutoffs for model w/ weight: 0.512, 0.538 

pdf(sprintf("./figures/ExpectedVsObservedExpr_GeneSet%s.pdf", 
            Geneset), width=9, height=6)
par(mar=c(5,4,1,1), bty="n", mfrow=c(2,3), cex=0.8)

# plot observed vs. expected expression, stratified by 3-quantiles
plot_log1p = function(x, y, ...) {
  smoothScatter(log(x+1e-5), log(y+1e-5), xlim=c(-5, 10), ylim=c(-5, 10), ...)
  legend("bottomright", bty="n",
         legend=sprintf("Pearson correlation = %.2f", cor(log(x+1e-5), log(y+1e-5))))
}

plot_log1p(c(predicted_bulk_nw)[p0 < p0_cutoffs[1]], c(bulk)[p0 < p0_cutoffs[1]],
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/o weight", main="low prob of being consistent")
plot_log1p(c(predicted_bulk_nw)[p0 >= p0_cutoffs[1] & p0 <= p0_cutoffs[2]], 
           c(bulk)[p0 >= p0_cutoffs[1] & p0 <= p0_cutoffs[2]], 
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/o weight", main="med prob of being consistent")
plot_log1p(c(predicted_bulk_nw)[p0 > p0_cutoffs[2]], c(bulk)[p0 > p0_cutoffs[2]], 
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/o weight", main="high prob of being consistent")

plot_log1p(c(predicted_bulk_w0)[p1 < p1_cutoffs[1]], c(bulk)[p1 < p1_cutoffs[1]], 
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/ weight", main="low prob of being consistent")
plot_log1p(c(predicted_bulk_w0)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]], 
           c(bulk)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]], 
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/ weight", main="med prob of being consistent")
plot_log1p(c(predicted_bulk_w0)[p1 > p1_cutoffs[2]], c(bulk)[p1 > p1_cutoffs[2]], 
     xlab="Predicted gene expression", ylab="Observed gene expression",
     sub="model w/ weight", main="high prob of being consistent")

dev.off()


setwd("./programs")
sessionInfo()

q(save="no")

