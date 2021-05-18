
#' Remove batch effects from a taxa read count table
#'
#' @import quantreg
#' @import cqrReg
#' @import glmnet
#' @import dplyr
#' @import doParallel
#' @import gplots
#' @import vegan
#' @import ade4
#' @import compositions
#' @import randomForest
#' @import ROCR
#' @import ape
#' @import GUniFrac
#'
#' @param tax_tab The taxa read count table, samples (row) by taxa (col).
#' @param batchid The batch indicator, must be a factor.
#' @param covariates The data.frame contains the key variable of interest and other covariates, e.g., data.frame(key, x1, x2).
#' @param batch_ref A character, the name of the reference batch, e.g.,``2''.
#' @param logistic_lasso A logical value, TRUE for L1-penalized logistic regression, FALSE for standard logistic regression; default is FALSE.
#' @param quantile_type A character, ``standard'' for standard quantile regression, ``lasso'' for L1-penalized quantile regression, ``composite'' for composite quantile regression; default is ``standard''.
#' @param simple_match A logical value, TRUE for using the simple quantile-quantile matching, FALSE for not; default is FALSE.
#' @param lambda_quantile A character, the penalization parameter in quantile regression if \code{quantile_type}=``lasso'' or ``composite''; only two choices ``2p/n'' or ``2p/logn'', where p is the number of expanded covariates and n is the number of non-zero read count; default is ``2p/n''.
#' @param interplt A logical value, TRUE for using the data-driven linear interpolation between zero and non-zero quantiles to stablize border estimates, FALSE for not; default is FALSE.
#' @param delta A real constant in (0, 0.5), determing the size of the interpolation window if interplt=TRUE, a larger delta leads to a narrower interpolation window; default is 0.4999.
#' @param taus A sequence of quantile levels, determing the ``precision'' of estimating conditional quantile functions; default is seq(0.005, 0.995, by=0.005).
#' @param num_core A real constant, the number of cores used for computing; default is 2.
#'
#' @details
#' \itemize{
#'   \item Choose \code{batch_ref} based on prior knowledge, or try several options, there is no default.
#'   \item The option ``composite'' of \code{quantile_type} is aggressive, use with caution.
#'   \item If choose \code{simple_match}=TRUE, \code{logistic_lasso}, \code{quantile_type}, \code{lambda_quantile}, \code{interplt} and \code{delta} won't take effect.
#'   \item Always use a fine grid of \code{taus} if the size of data is adequate.
#' }
#'
#' @return The corrected taxa read count table, samples (row) by taxa (col).
#'
#' @references
#' \itemize{
#'   \item Ling, W. et al. (2021+). ConQuR: batch effects removal for microbiome data in large-scale epidemiology studies via conditional quantile regression.
#'   \item Ling, W. et al. (2020+). Statistical inference in quantile regression for zero-inflated outcomes. Statistica Sinica.
#'   \item Machado, J.A.F., Silva, J.S. (2005). Quantiles for counts. Journal of the American Statistical Association 100(472), 1226–1237.
#'   \item Koenker, R. & Bassett Jr, G. (1978). Regression quantiles. Econometrica: journal of the Econometric Society, 33-50.
#'   \item Koenker, R. (2005). Econometric Society Monographs: Quantile Regression. New York: Cambridge University.
#'   \item Zou, H. & Yuan, M. (2008). Composite quantile regression and the oracle model selection theory. The Annals of Statistics 36, 1108-1126.
#' }
#'
#' @export

ConQuR <- function(tax_tab, batchid, covariates,
                   batch_ref,
                   logistic_lasso=F, quantile_type="standard", simple_match=F,
                   lambda_quantile="2p/n", interplt=F,
                   delta=0.4999, taus=seq(0.005, 0.995, by=0.005), num_core=2){

  # relevel batch id
  batchid = relevel(batchid, ref=batch_ref)


  #### by simple quantile-quantile matching is chosen ####

  if (simple_match == T){

    registerDoParallel(num_core)

    tax_new = foreach (ll=1:ncol(tax_tab), .combine=cbind) %do%{
      y = tax_tab[, ll]
      simple_QQ(y=y, batchid=batchid, batch_ref=batch_ref, taus=taus)
    }
    tax_new[tax_new < 0] = 0

    rownames(tax_new) = rownames(tax_tab)
    colnames(tax_new) = colnames(tax_tab)
    return(tax_new)
  }


  #### otherwise, correct data via regression ####

  ### process data ###
  X = data.frame(covariates, batchid)
  X_span = model.matrix(~., X)[,-1]

  X_correct = X
  X_correct$batchid = batch_ref
  X_correct$batchid = factor(X_correct$batchid)

  X_span_correct = X_span
  X_span_correct[, grep( "batchid", colnames(X_span_correct) )] = 0

  ### correct each of the taxa ###
  registerDoParallel(num_core)

  tax_new = foreach (ll=1:ncol(tax_tab), .combine=cbind) %do%{
    y = tax_tab[, ll]
    ConQuR_each(y=y, X=X, X_span=X_span, X_correct=X_correct, X_span_correct=X_span_correct,
                delta=delta, taus=taus, logistic_lasso=logistic_lasso, quantile_type=quantile_type, lambda_quantile=lambda_quantile, interplt=interplt)
  }
  tax_new[tax_new < 0] = 0

  rownames(tax_new) = rownames(tax_tab)
  colnames(tax_new) = colnames(tax_tab)
  return(tax_new)

}



#' Tune over variations of ConQuR
#'
#' @param tax_tab The taxa read count table, samples (row) by taxa (col).
#' @param batchid The batch indicator, must be a factor.
#' @param covariates The data.frame contains the key variable of interest and other covariates, e.g., data.frame(key, x1, x2).
#' @param batch_ref_pool A vector of characters, the candidates for reference batch, e.g., c(``0'', ``2'').
#' @param logistic_lasso_pool A vector of logical values, whether or not using the L1-penalized logistic regression, e.g., c(T, F).
#' @param quantile_type_pool A vector of characters, the candidates for quantile regression type, e.g., c(``standard'', ``lasso'').
#' @param simple_match_pool A vector of logical values, whether or not using the simple quantile-quantile matching, e.g., c(T, F).
#' @param lambda_quantile_pool A vector of characters, the candidates for the penalization parameter in quantile regression (``lasso'' or ``composite''), e.g., c(NA, ``2p/n'', ``2p/logn'').
#' @param interplt_pool A vector of logical values, whether or not using the data-driven linear interpolation between zero and non-zero quantiles, e.g., c(T, F).
#' @param frequencyL A real constant between 0 and 1, the lower bound of prevalence that needs tuning.
#' @param frequencyU A real constant between 0 and 1, the upper bound of prevalence that needs tuning.
#' @param cutoff A real constant, the grid size of prevalence for tuning; default is 0.1.
#' @param delta A real constant in (0, 0.5), determing the size of the interpolation window if interplt=TRUE, a larger delta leads to a narrower interpolation window; default is 0.4999.
#' @param taus A sequence of quantile levels, determing the ``precision'' of estimating conditional quantile functions; default is seq(0.005, 0.995, by=0.005).
#' @param num_core A real constant, the number of cores used for computing; default is 2.
#'
#' @details
#' \itemize{
#'   \item ``original'', i.e., the original data without correction is always a default candidate.
#'   \item If ``standard'' is one candidate for \code{quantile_type_pool}, always include NA as one candidate for \code{lambda_quantile_pool}.
#'   \item Be cautious with candidate ``composite'' for \code{quantile_type_pool}, the underlying assumption is strong and the computation might be slow.
#'   \item The tuning procedure finds the local optimal in each cutoff. If \code{frequencyL}=0.2, \code{frequencyU}=0.5 and \code{cutoff}=0.1, the functions determines the combination achieving maximum removal of batch variations on taxa present in 20\%-30\%, ..., 40\%-50\% of the samples, respectively.
#'   \item The same reference batch is used across taxa in the final optimal corrected table.
#' }
#'
#' @return A list
#' \itemize{
#'   \item tax_final - The optimal corrected taxa read count table, samples (row) by taxa (col).
#'   \item method_final - A table summarizing variations of ConQuR chosen for each prevalence cutoff.
#'}
#'
#' @references
#' \itemize{
#'   \item Ling, W. et al. (2021+). ConQuR: batch effects removal for microbiome data in large-scale epidemiology studies via conditional quantile regression
#'   \item Ling, W. et al. (2020+). Statistical inference in quantile regression for zero-inflated outcomes. Statistica Sinica.
#'   \item Machado, J.A.F., Silva, J.S. (2005). Quantiles for counts. Journal of the American Statistical Association 100(472), 1226–1237.
#'   \item Koenker, R. & Bassett Jr, G. (1978). Regression quantiles. Econometrica: journal of the Econometric Society, 33-50.
#'   \item Koenker, R. (2005). Econometric Society Monographs: Quantile Regression. New York: Cambridge University.
#'   \item Zou, H. & Yuan, M. (2008). Composite quantile regression and the oracle model selection theory. The Annals of Statistics 36, 1108-1126.
#'   \item Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#' }
#'
#' @export

Tune_ConQuR <- function(tax_tab, batchid, covariates,
                        batch_ref_pool,
                        logistic_lasso_pool,
                        quantile_type_pool,
                        simple_match_pool,
                        lambda_quantile_pool,
                        interplt_pool,
                        frequencyL,
                        frequencyU,
                        cutoff=0.1, delta=0.4999, taus=seq(0.005, 0.995, by=0.005), num_core=2){

  # prepare the method list, taxa pool corrresponding to the grid of frequency, and result table
  method1 = expand.grid(logistic=logistic_lasso_pool, quantile=quantile_type_pool, lambda=lambda_quantile_pool, interplt=interplt_pool)
  method1 = method1[-c( which(method1$quantile=="standard" & !is.na(method1$lambda)), which(method1$quantile!="standard" & is.na(method1$lambda)) ), ]
  method2 = any(simple_match_pool == T)

  freq_grid = seq(from=frequencyL, to=frequencyU, by=cutoff)
  prate = apply(tax_tab, 2, function(z){length( which(z > 0) ) / nrow(tax_tab) })

  cut_list = NULL
  for (ii in 1:(length(freq_grid)-1)){
    start = freq_grid[ii]
    end = freq_grid[ii+1]

    cut_list[[ii]] = which(prate > start & prate <= end)
  }

  R2_initial = data.frame(method=c("original", apply(method1, 1, paste, collapse="_")), batch_R2=NA)
  if (method2 == T) R2_initial = rbind(R2_initial, c("simple", NA))

  # greedy experiment on reference batch and on each cutoff of frequency
  tax_new_list = vector(mode="list", length=length(batch_ref_pool))
  method_chosen_list = vector(mode="list", length=length(batch_ref_pool))

  names(tax_new_list) = names(method_chosen_list) = batch_ref_pool

  # search on reference batch
  for (current_ref in 1:length(batch_ref_pool)){

    tax_new = matrix(ncol=ncol(tax_tab), nrow=nrow(tax_tab))
    method_chosen = NULL

    # search on each cutoff of frequency
    for (current_cutoff in 1:length(cut_list)){

      R2 = R2_initial
      tab_list = NULL

      if (length( cut_list[[current_cutoff]] ) == 0) next

      # benchmark -- on original data
      tab = tax_tab[, cut_list[[current_cutoff]]]
      tab_list[[1]] = tab

      if (length( cut_list[[current_cutoff]] ) == 1) tab = matrix(tab, nrow=nrow(tax_tab))
      index = which( apply(tab, 1, sum) > 0 )
      R2[1, 2] = adonis(formula = tab[index, ] ~ batchid[index])$aov.tab[1, 5]

      # do correction for all combination, record results and compute PERMANOVA R2
      for (current_method in 1:nrow(method1)){
        tab_list[[1+current_method]] = ConQuR(tax_tab=tab, batchid=batchid, covariates=covariates,
                                              batch_ref=as.character(batch_ref_pool[current_ref]),
                                              logistic_lasso=method1[current_method, 'logistic'],
                                              quantile_type=as.character(method1[current_method, 'quantile']),
                                              simple_match=F,
                                              lambda_quantile=as.character(method1[current_method, 'lambda']),
                                              interplt=method1[current_method, 'interplt'],
                                              delta=delta, taus=taus, num_core=num_core)

        if (length( cut_list[[current_cutoff]] ) == 1) tab_list[[1+current_method]] = matrix(tab_list[[1+current_method]], nrow=nrow(tax_tab))
        index = which( apply(tab_list[[1+current_method]], 1, sum) > 0 )
        R2[1+current_method, 2] = adonis(formula = tab_list[[1+current_method]][index, ] ~ batchid[index])$aov.tab[1, 5]
      }

      if (method2 == T){
        tab_list[[1+nrow(method1)+1]] = ConQuR(tax_tab=tab, batchid=batchid, covariates=covariates,
                                               batch_ref=as.character(batch_ref_pool[current_ref]),
                                               logistic_lasso=F,
                                               quantile_type="standard",
                                               simple_match=T,
                                               lambda_quantile="2p/n",
                                               interplt=F,
                                               delta=delta, taus=taus, num_core=num_core)

        if (length( cut_list[[current_cutoff]] ) == 1) tab_list[[1+nrow(method1)+1]] = matrix(tab_list[[1+nrow(method1)+1]], nrow=nrow(tax_tab))
        index = which( apply(tab_list[[1+nrow(method1)+1]], 1, sum) > 0 )
        R2[1+nrow(method1)+1, 2] = adonis(formula = tab_list[[1+nrow(method1)+1]][index, ] ~ batchid[index])$aov.tab[1, 5]
      }

      # find the optimal choice
      index_optimal = which.min(R2$batch_R2)

      if (length( cut_list[[current_cutoff]] ) == 1){
        tax_new[, cut_list[[current_cutoff]]] = tab_list[[index_optimal]]
      } else{
        for (col in 1:length( cut_list[[current_cutoff]] )){
          tax_new[, cut_list[[current_cutoff]][col]] = tab_list[[index_optimal]][, col]
        }
      }

      method_chosen[current_cutoff] = R2$method[index_optimal]

    }

    # correct the remaining taxa by the default method -- standard logistic + standard quantile + no interpolation
    cut_list_remaining = which(prate <= frequencyL | prate > frequencyU)

    if (!length( cut_list_remaining ) == 0){

      tab = tax_tab[, cut_list_remaining]
      if (length( cut_list_remaining ) == 1) tab = matrix(tab, nrow=nrow(tax_tab))

      tax_new_remaining = ConQuR(tax_tab=tab, batchid=batchid, covariates=covariates,
                                 batch_ref=as.character(batch_ref_pool[current_ref]),
                                 logistic_lasso=F,
                                 quantile_type="standard",
                                 simple_match=F,
                                 lambda_quantile="2p/n",
                                 interplt=F,
                                 delta=delta, taus=taus, num_core=num_core)

      for (col in 1:length( cut_list_remaining )){
        tax_new[, cut_list_remaining[col]] = tax_new_remaining[, col]
      }

    }

    tax_new_list[[current_ref]] = tax_new
    method_chosen_list[[current_ref]] = method_chosen

  }

  # determine the reference batch
  R2_ref = NULL
  for (current_ref in 1:length(batch_ref_pool)){
    index = which( apply(tax_new_list[[current_ref]], 1, sum) > 0 )
    R2_ref[current_ref] = adonis(formula = tax_new_list[[current_ref]][index, ] ~ batchid[index])$aov.tab[1, 5]
  }

  index_final = which.min(R2_ref)

  tax_final = tax_new_list[[index_final]]
  rownames(tax_final) = rownames(tax_tab)
  colnames(tax_final) = colnames(tax_tab)

  method_final = matrix(ncol=length(cut_list), nrow=7)
  rownames(method_final) = c("batch_ref", "original", "simple_match", "logistic_lasso", "quantile_type", "lambda", "interplt")
  colnames(method_final) = paste0(freq_grid[-length(freq_grid)], "-", freq_grid[-1])

  method_temp = method_chosen_list[[index_final]]
  for (ii in 1:length(cut_list)){
    if (is.na(method_temp[ii])){
      method_final[, ii] = rep(NA, 7)
    } else if (method_temp[ii] == "original"){
      method_final[, ii] = c(NA, T, rep(NA, 5))
    } else if (method_temp[ii] == "simple"){
      method_final[, ii] = c(names(method_chosen_list)[index_final], F, T, rep(NA, 4))
    } else{
      method_final[, ii] = c(names(method_chosen_list)[index_final], F, F, unlist( strsplit(method_temp[ii], '_') ) )
    }
  }

  return(list(tax_final=tax_final, method_final=method_final))

}

