
#' Example data, a taxa read count table, with batchid, key variable and covariates
#'
#' A dataset containing 100 taxa from 3 batches, key variable is sbp, with covariates, sex, race and age
#'
#' @format A taxa read count (273 samples by 100 taxa), batchid and the metadata:
#' \describe{
#'   \item{batchid}{factor, with levels 0, 1, 2}
#'   \item{sbp}{key variable, systolic blood pressure, continuous variable}
#'   \item{sex}{covariate 1, binary variable}
#'   \item{race}{covariate 2, binary variable}
#'   \item{age}{covariate 3, continuous variable}
#' }
"Sample_Data"
