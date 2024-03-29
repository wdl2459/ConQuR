% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supporting_functions.R
\name{PERMANOVA_R2}
\alias{PERMANOVA_R2}
\title{PERMANOVA R2 of batch and variable of interest}
\usage{
PERMANOVA_R2(TAX, batchid, covariates, key_index)
}
\arguments{
\item{TAX}{The taxa read count table, samples (row) by taxa (col).}

\item{batchid}{The batch indicator, must be a factor.}

\item{covariates}{The data.frame contains the key variable of interest and other covariates.}

\item{key_index}{An integer, location of the variable of interest in \code{covariates}.}
}
\value{
A list
\itemize{
  \item tab_count - A table summarizing PERMANOVA R2 computed on the original taxa read count table in Bray-Curtis dissimilarity.
  \item tab_rel - A table summarizing PERMANOVA R2 computed on the corresponding relative abundance table in Euclidean dissimilarity (Aitchison dissimilarity).
}
}
\description{
PERMANOVA R2 of batch and variable of interest
}
\details{
Three PERMANOVA R2 will be computed: (1) the standard one (adnois), (2) on euclidified dissimilarities (adonis2, sqrt.dist=T), and (3) with a constant added to the non-diagonal dissimilarities such that all eigenvalues are non-negative in the underlying PCoA (adonis2, add=T).
}
\references{
\itemize{
  \item Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
}
}
