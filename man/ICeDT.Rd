\name{ICeDT}
\alias{ICeDT}
\title{
ICeDT Model Fit
}
\description{
Estimate cell type composition of tumor samples using gene expression data of bulk tumor samples as well as cell type-specific expression from a set of signature genes. 
}
\usage{
ICeDT(Y, Z, tumorPurity = NULL, refVar = NULL, rhoInit = NULL, 
maxIter_prop = 100, maxIter_PP = 100, rhoConverge = 1e-03)
}
\arguments{
  \item{Y}{
A matrix of normalized gene expression (e.g. FPKM, TPM) from bulk tumor samples. It is of dimensions nG by nS, where nG is the number of genes and nS is the number of bulk tumor samples. 
}
  \item{Z}{
A matrix of estimated cell type-specific gene expression with each row for a gene and each column for a cell type. The colnames should be cell types. 
}
  \item{tumorPurity}{
A vector of length nS which contains the tumor purity for each sample in \code{Y}. Order must match that given in the columns of \code{Y}.
}
  \item{refVar}{
A matrix of the same dimension as \code{Z}, which gives cell type-specific variance for each gene (a row) and each sample (a column). 
}
  \item{rhoInit}{
A matrix of initial values or rho. Each column is the cell type composition for each sample in \code{Y}, and each row is one cell type. The rows should match to the columns of \code{Z}, and the columns should match to the columns of \code{Y}. 
}
  \item{maxIter_prop}{
Maximum number of iterations used for updating proportion parameters within a single EM iteration.
}
  \item{maxIter_PP}{
Maximum number of iterations used for updating posterior probabilities. 
}
  \item{rhoConverge}{
The maximum difference in estimated cell type proportions for defining convergence.
}
}
\details{
Estimate cell type composition using a mixture of log-normal distribution, where one mixture component includes consistent genes that behave consistently in bulk tumor samples and refence expression matrix, and the other mixture component includes aberrant genes that behave differently between bulk tumor samples and reference expression matrix. 
}
\value{
A list item containing:

- rho: A matrix of estimated immune cell type abundances. One row per subject, one column per cell type.

- sigma2C: The estimated variance parameter for consistent marker genes. One parameter for each subject in the order specified in the rows of rho. 

- sigma2A: The estimated variance parameter for Aberrant genes. One parameter for each subject in the order specified in the rows of rho. 

- cProp: The estimated proportion of consistent genes, one value for each subject in the order established by rho. 

- cProb: A matrix of probabilities that a gene is consistent given the model estimates (one row per gene, one column per subject). 

}
\references{
Wilson, Douglas R; Jin, Chong; Ibrahim, Joseph G.; Sun, Wei "ICeD-T Provides Accurate Estimates of Immune Cell Abundance in Tumor Samples by Allowing for Aberrant Gene Expression Patterns" <<In Review>>, 2019.
}
\author{
Douglas Roy Wilson, Jr., Wei Sun
}
\keyword{ methods }
