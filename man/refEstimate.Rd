\name{refEstimate}
\alias{refEstimate}
\title{
 estimate reference matrix of cell type-specific expression. 
}
\description{
estimate reference matrix of cell type-specific expression using observed gene expression from purified samples of a set of pre-selected marker genes. 
}
\usage{
refEstimate(X, borrow4SD = TRUE)
}
\arguments{
  \item{X}{
A matrix of observed cell type-specific gene expression with each row for a gene and each column for a purified sample.  There should be at least two samples per cell type.
}
  \item{borrow4SD}{
Indicator (TRUE/FALSE) of whether or not to borrow information across cell types in order to estimate the standard deviations in pure reference samples.
}
}
\details{
 a reference matrix is a matrix of gene expression where each row corresponds to a gene and each column corresponds to a cell type. The genes should be pre-selected so that they have cell type-specific expression. Cell type-specific mean expression of each gene is estimated assuming the expression of one gene in each cell type follows a log-normal distribution. 
}
\value{
A list
\item{refMat }{ expected gene expression (in linear scale) for each gene (row) and each cell type (column). }
\item{refVar }{ variance of gene expression (in log scale) for each gene (row) and each cell type (column). }
\item{ctMu }{ expected gene expression (in log scale) for each gene (row) and each cell type (column). }
\item{ctVar }{ variance of gene expression (in log scale) for each gene (row) and each cell type (column), after borrowing information across cell types if \code{borrow4SD} is TRUE. }
}
\examples{
set.seed(666)
data(mean_var_relation)
simData = simFunc(nS=135, nG=250, nrep=5, nCT=5, pctAb=20,
                  meanVar_Rel=mean_var_relation)
X = simData$PGene_Exp
dim(X)
X[1:2,1:7]
X = X[,which(colnames(X) != "tumor")]

refE = refEstimate(X)
names(refE)

refE$refMat[1:2,]
refE$refVar[1:2,]
refE$ctMu[1:2,]
refE$ctVar[1:2,]
}
\author{
Wei Sun, Douglas Wilson
}
\keyword{ utilities }
