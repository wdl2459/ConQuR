% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supporting_functions.R
\name{RF_Pred_Multiclass}
\alias{RF_Pred_Multiclass}
\title{Predict multiclass variables based on a taxa read count table by random forest}
\usage{
RF_Pred_Multiclass(TAX, factor, fold = 5, seed = 2020)
}
\arguments{
\item{TAX}{The taxa read count table, samples (row) by taxa (col).}

\item{factor}{The multiclass variable to predict, e.g., the key variable, never smoker/former smoker/current smoker, must be a factor.}

\item{fold}{The number of folds; default is 5.}

\item{seed}{The seed to generate fold indices for samples; default is 2020.}
}
\value{
A list
\itemize{
  \item pred - A table summarizing the predicted probabilities and true labels for all samples.
  \item cross_entropy_across_fold - mean cross-entropy across folds.
}
}
\description{
Predict multiclass variables based on a taxa read count table by random forest
}
