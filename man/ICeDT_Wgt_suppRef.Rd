\name{ICeDT_Wgt_suppRef}
\alias{ICeDT_Wgt_suppRef}
\title{
  ICeD-T Model Fit (Weights, Supplied Reference)
}
\description{
  This function uses a supplied immune reference and mixed expression profiles to estimate the cell type composition of tumor samples. This function can utilize information regarding tumor purity if available.
}
\usage{
  ICeDT_Wgt_suppRef(Y,fixedCT_rho,useRho=FALSE,RefMat,RefVar,
						  varLog = FALSE, limitWgt = TRUE, limitQuant = c(0.10,0.90),
						  userWgt=NULL, useIQR=FALSE,maxIter_prop = 100,
              maxIter_PP=100,RhoConv_CO = 1e-4, Subj_CO)
}
\arguments{
  \item{Y}{
    A matrix of normalized mixture expression profiles (e.g. FPKM, TPM) of dimensions nG by nS. This matrix contains one row per gene and one column per mixture expression profile.
  }
  \item{fixedCT_rho}{
    Vector of length nS which contains the cell type proportion for fixedCT for each sample. Order must match that given in the columns of Y.
  }
  \item{useRho}{
    useRho is a boolean variable (TRUE = use fixed cell type proportion information, FALSE = estimate fixed cell type proportion).
  }
  \item{RefMat}{
    A reference matrix containing one row per gene and one column for each assessed immune cell type (dimensions nG by K). 
  }
   \item{RefVar}{
    A reference variance matrix containing one row per gene and one column for each assessed immune cell type (dimensions nG by K). Different from RefMat which characterizes the "average expression", this reference contains the "variance in expression." Order of columns and rows must match that in RefMat.
  }
  \item{varLog}{
	A boolean variable (TRUE=RefVar is computed for log-scale expression, FALSE=RefVar is computed from normal-scale expression).
  }
  \item{limitWgt}{
    A boolean variable to specify whether or not to control variance weights for robustness considerations (TRUE = constrain model weights, FALSE = do not constrain model weights)
  }
  \item{limitQuant}{
	A vector of length two specifying the upper and lower quantiles used to limit model weights. The first element specifies the lower quantile to use for weight replacment (Bottom limitQuant[1]*100\% replaced). The second element specifies the upper quantile to use for weight replacement (Weights above limitQuant[2]*100\% replaced).
  }
  \item{userWgt}{
	If user specifies useIQR = 6, then user specified model weights are used. This variable contains these weights.
  }
  \item{useIQR}{
    A categorical variable (0 or 6) which gives the form of the weights to use. useIQR = 0 utilizes the weight structure suggested in the ICeD-T manuscript. useIQR = 6 utilizes user supplied weights.
  }
  \item{maxIter_prop}{
    Maximum number of iterations used for updating proportion parameters within a single EM iteration.
  }
  \item{maxIter_PP}{
    Maximum number of EM iterations allowed.
  }
  \item{RhoConv_CO}{
    The maximum difference in estimated cell type proportions for defining convergence.
  }
  \item{Subj_CO}{
    The number of subjects who must meet the RhoConv_CO cutoff in order to stop iterating the model.
  }
}
\details{
  The following model is designed to compute ICeD-T model fits using no weights and a supplied reference. In the vernacular established above, the fixed cell type is typically a tumor cell type and fixedCT_rho thus represents tumor purity.
}
\value{
  A list item containing:
    
  - $IC_Abundance : One row per subject, one column per cell type containing ICeD-T estimates.
  
  - $Fixed_CellType: The label defining the fixed cell type.
  
  - $Sigma2M : The estimated variance parameter for consistent marker genes. One parameter for each subject in the order specified in IC_Abundance rows. 
  
  - $Sigma2A : The estimated variance parameter for Aberrant genes. One parameter for each subject in the order specified in IC_Abundance rows. If Sigma2M > Sigma2A for a subject, this indicates that these parameters must be switched, as must P_Consistent and PP_Consistent.
  
  - $P_Consistent: The estimated proportion of consistent genes, one value for each subject in the order established by IC_Abundance. If Sigma2M > Sigma2A and one is interested in Aberrant genes, "P_Consistent" must be replaced by 1-P_Consistent for such subjects. 
  
  - $PP_Consistent: A matrix of conditional probabilities that a given gene is Aberrant given the model estimates (one row per gene, one column per subject). For subjects where Sigma2M > Sigma2A, the corresponding column of PP_Consistent must be replaced by 1-PP_Consistent.
  
}
\references{
  Wilson, Douglas R; Ibrahim, Joseph G.; Sun, Wei ''ICED-T Provides Accurate Estimates of Immune Cell Abundance in Tumor Samples by Allowing for Aberrant Gene Expression Patterns'' <<In Review>>, 2018.
}
\author{
  Douglas Roy Wilson, Jr.
}
\note{
  Additionally, note the proper treatment of $P_Consistent and $PP_Consistent if Sigma2M > Sigma2A for any subject. 
}

\keyword{ Expression Deconvolution }

