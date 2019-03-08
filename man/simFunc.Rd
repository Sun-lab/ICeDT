\name{simFunc}
\alias{simFunc}
\title{
Simulation of gene expression data while allowing aberrant genes.
}
\description{
Simulate gene expression data of both purified samples and bulk smaples of mixed cell types. 
}
\usage{
simFunc(nS, nG, nrep, nCT, pctAb, meanVar_Rel, abMech = 1)
}
\arguments{
  \item{nS}{
  number of bulk samples to be simulated. 
}
  \item{nG}{
  number of genes. 
}
  \item{nrep}{
  number of purified samples per cell type. 
}
  \item{nCT}{
  number of cell types, including tumor cell type. 
}
  \item{pctAb}{
  the percentage of genes that are aberrant. 
}
  \item{meanVar_Rel}{
  a data frame storing the mean/var relationship used for simulation, with element x and y being mean and squareroot of sd, respectively. 
}
  \item{abMech}{
  The mechanism which defines abberancy. currently only value 1 is allowed. 
}
}
\details{
  With abMech=1, pctAb of genes are randomly selected as abberant, and their expression are up/down regulated or allow these genes to be expressed in tumor samples. 
}
\value{
A list of the following elements
\item{MGene_Exp }{ matrix of size nG x nS which contains the gene expression in mixture samples. One Row per gene and one column per sample. }

\item{PGene_Exp }{ matrix of size nG x (nCT*nrep) which contains the gene expression for purified reference samples. }
\item{cellType }{ cell type for each subject. Order matches columns in PGene_Exp. }
\item{aberrant }{ a vector indicating whether or not a gene was simulated as an aberrant gene. }
\item{rho }{ matrix of nS x nCT: cell type proportions. }
\item{tumor_purity }{ Vector of length nS containing the tumor purities for each subject. Order matches columns in MGene_Exp. }
}
\examples{
set.seed(666)
data(mean_var_relation)
simData = simFunc(nS=135, nG=250, nrep=5, nCT=5, pctAb=20,
                  meanVar_Rel=mean_var_relation)
names(simData)
simData$tumor_purity[1:5]
simData$PGene_Exp[1:2,1:7]
}
\references{
Wilson, Douglas R; Jin, Chong; Ibrahim, Joseph G.; Sun, Wei "ICED-T Provides Accurate Estimates of Immune Cell Abundance in Tumor Samples by Allowing for Aberrant Gene Expression Patterns" <<In Review>>, 2018.
}
\author{
Douglas Roy Wilson, Jr., Wei Sun
}

\keyword{ utilities }
