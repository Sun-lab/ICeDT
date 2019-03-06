\name{ICeDT}
\alias{ICeDT}
\title{
ICeDT Model Fit
}
\description{
This function uses purified, immune reference sample data and mixed expression profiles to estimate the cell type composition of tumor samples. This function can utilize information regarding tumor purity if available.
}
\usage{
ICeDT(Y, X, refMat=TRUE, tumorPurity=NULL, borrow4SD = TRUE,
maxIter_prop = 100, maxIter_PP = 100, rhoConverge = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A matrix of normalized mixture expression profiles (e.g. FPKM, TPM) of dimensions nG by nS. This matrix contains one row per gene and one column per mixture expression profile.
}
  \item{X}{
A matrix of observed or estimated cell type-specific gene expression with each row for a gene. If it is the observed gene expression, each column is a purified sample, and there should be at least two samples per cell type. If it is the latter, each column is one cell type. The colnames should be cell types. 
}
  \item{refMat}{
Indicate whether X is estimated reference matrix (TRUE) or observed gene expression (FALSE). 
}
  \item{tumorPurity}{
A vector of length nS which contains the tumor purity for each sample in \code{Y}. Order must match that given in the columns of \code{Y}.
}
  \item{borrow4SD}{
Indicator (TRUE/FALSE) of whether or not to borrow information across cell types in order to estimate the standard deviations in pure reference samples.
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
The following model is designed to compute ICeD-T model fits when purified reference sample data is available and when one does not wish to use previously defined references. In the vernacular established above, the fixed cell type is typically a tumor cell type and fixedCT_rho thus represents tumor purity.

In order to utilize this function, one must create at least two dummy "pure" tumor samples and append them to X with an appropriate addition to fixedCT. 
}
\value{
A list item containing:

- rho: A matrix of estimated immune cell type abundances. One row per subject, one column per cell type.

- fixedCT: The label defining the fixed cell type.

- sigma2M: The estimated variance parameter for consistent marker genes. One parameter for each subject in the order specified in rho rows. 

- sigma2A: The estimated variance parameter for Aberrant genes. One parameter for each subject in the order specified in rho rows. If sigma2M > sigma2A for a subject, this indicates that these parameters must be switched, so does the values of P_Consistent and PP_Consistent.

- Z: The utilized reference matrix with one row per gene and one column per immune cell type in fixedCT. 

- CT_var: The variance profiles of the established immune cell types from purified reference samples.

- P_Consistent: The estimated proportion of consistent genes, one value for each subject in the order established by rho. If sigma2M > sigma2A and one is interested in Aberrant genes, "P_Consistent" must be replaced by 1-P_Consistent for such subjects. 

- PP_Consistent: A matrix of conditional probabilities that a given gene is Aberrant given the model estimates (one row per gene, one column per subject). For subjects where sigma2M > sigma2A, the corresponding column of PP_Consistent must be replaced by 1-PP_Consistent.

}
\references{
Wilson, Douglas R; Jin, Chong; Ibrahim, Joseph G.; Sun, Wei "ICED-T Provides Accurate Estimates of Immune Cell Abundance in Tumor Samples by Allowing for Aberrant Gene Expression Patterns" <<In Review>>, 2018.
}
\author{
Douglas Roy Wilson, Jr., Wei Sun
}
\note{
Please note the addition of "dummy" pure tumor expression profiles. Additionally, note the proper treatment of $P_Consistent and $PP_Consistent if sigma2M > sigma2A for any subject. 
}

\keyword{ Expression Deconvolution}
