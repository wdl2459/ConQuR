
########## subsequent analyses ##########


#' PERMANOVA R2 of batch and variable of interest
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param batchid The batch indicator, must be a factor.
#' @param covariates The data.frame contains the key variable of interest and other covariates.
#' @param key_index An integer, location of the variable of interest in \code{covariates}.
#'
#' @details Three PERMANOVA R2 will be computed: (1) the standard one (adnois), (2) on euclidified dissimilarities (adonis2, sqrt.dist=T), and (3) with a constant added to the non-diagonal dissimilarities such that all eigenvalues are non-negative in the underlying PCoA (adonis2, add=T).
#'
#' @return A list
#' \itemize{
#'   \item tab_count - A table summarizing PERMANOVA R2 computed on the original taxa read count table.
#'   \item tab_rel - A table summarizing PERMANOVA R2 computed on the corresponding relative abundance table.
#'}
#'
#' @references
#' \itemize{
#'   \item Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#' }
#'
#' @export

PERMANOVA_R2 <- function(TAX, batchid, covariates, key_index){

  tab_count = tab_rel = matrix(ncol=3, nrow=2)
  colnames(tab_count) = colnames(tab_rel) = c("standard", "sqrt.dist=T", "add=T")
  rownames(tab_count) = rownames(tab_rel) = c("batch", "key")

  # count
  tab_count[1,1] = adonis(formula = TAX ~ batchid)$aov.tab[1, 5]
  tab_count[1,2] = adonis2(formula = TAX ~ batchid, sqrt.dist=TRUE)$R2[1]
  tab_count[1,3] = adonis2(formula = TAX ~ batchid, add=TRUE)$R2[1]

  tab_count[2,1] = adonis(formula = TAX ~ ., data=data.frame(covariates[, key_index]))$aov.tab[1, 5]
  tab_count[2,2] = adonis2(formula = TAX ~ ., data=data.frame(covariates[, key_index]), sqrt.dist=TRUE)$R2[1]
  tab_count[2,3] = adonis2(formula = TAX ~ ., data=data.frame(covariates[, key_index]), add=TRUE)$R2[1]

  # rel
  TAX_temp = t( apply(TAX, 1, function(z){
    z / sum(z)
  }) )

  tab_rel[1,1] = adonis(formula = TAX_temp ~ batchid)$aov.tab[1, 5]
  tab_rel[1,2] = adonis2(formula = TAX_temp ~ batchid, sqrt.dist=TRUE)$R2[1]
  tab_rel[1,3] = adonis2(formula = TAX_temp ~ batchid, add=TRUE)$R2[1]

  tab_rel[2,1] = adonis(formula = TAX_temp ~ ., data=data.frame(covariates[, key_index]))$aov.tab[1, 5]
  tab_rel[2,2] = adonis2(formula = TAX_temp ~ ., data=data.frame(covariates[, key_index]), sqrt.dist=TRUE)$R2[1]
  tab_rel[2,3] = adonis2(formula = TAX_temp ~ ., data=data.frame(covariates[, key_index]), add=TRUE)$R2[1]

  return(list(tab_count=tab_count, tab_rel=tab_rel))

}


### PCoA plots: Bray-Curtis, GUniFrac (unweighted, weighted, generalized), Aitchinson

#' Stratified PCoA plots
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param factor The variable for stratification, e.g., batchid or the variable of interest, must be a factor.
#' @param sub_index A vector of sample indices, to restrict the analysis to a subgroup of samples, e.g., c(1:5, 15:20); default is NULL.
#' @param dissimilarity The dissimilarity type, ``Bray'' for Bray-Curtis dissimilarity, ``Aitch'' for Aitchison dissimilarity, ``GUniFrac'' for generalized UniFrac dissimilarity; default is ``Bray''.
#' @param GUniFrac_type The generalized UniFrac type, ``d_1'' for weighted UniFrac, ``d_UW'' for unweighted UniFrac, ``d_VAW'' for variance adjusted weighted UniFrac, ``d_0'' for generalized UniFrac with alpha 0, ``d_0.5'' for generalized UniFrac with alpha 0.5; default is ``d_0.5''.
#' @param tree The rooted phylogenetic tree of R class ``phylo'', must be provided when \code{dissimilarity}=``GUniFrac''; default is NULL.
#' @param main The title of plot; default is NULL.
#' @param aa A real number, the character size for the title.
#'
#' @return Print a PCoA plot.
#'
#' @references
#' \itemize{
#'   \item Chen, J., & Chen, M. J. (2018). Package ‘GUniFrac’. The Comprehensive R Archive Network (CRAN).
#' }
#'
#' @export

Plot_PCoA <- function(TAX, factor, sub_index=NULL, dissimilarity="Bray", GUniFrac_type="d_0.5", tree=NULL, main=NULL, aa=1.5){

  nfactor = length(table(factor))
  if (is.null(sub_index)){
    sub_index = seq(ncol(TAX))
  }

  if (dissimilarity == "Bray"){

    index = which( apply(TAX[, sub_index], 1, sum) > 0 )
    bc =  vegdist(TAX[index, sub_index])
    MDS = cmdscale(bc, k=4)
    s.class(MDS, fac = as.factor(factor[index]), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else if (dissimilarity == "Aitch"){

    Z = as.matrix(clr(as.matrix(TAX[, sub_index])+0.5))
    MDS = cmdscale(vegdist(Z, method = "euclidean"), k=4)
    s.class(MDS, fac = as.factor(factor), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else if (dissimilarity == "GUniFrac"){

    index = which( apply(TAX[, sub_index], 1, sum) > 0 )
    unifracs = GUniFrac(TAX[index, sub_index], tree, alpha=c(0, 0.5, 1))$unifracs

    if (!GUniFrac_type %in% c("d_1", "d_UW", "d_VAW", "d_0", "d_0.5")) stop("Please use one of d_1, d_UW, d_VAW, d_0, d_0.5 for GUniFrac dissimilarity.")

    d = unifracs[, , GUniFrac_type]
    MDS = cmdscale(d, k=4)
    s.class(MDS, fac = as.factor(factor[index]), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else{
    stop("Please use one of Bray, Aitch or GUniFrac as the dissimilarity.")
  }

}


### RF - predict key variable

# Predict binary key variable using random forest, k-fold cross-validation (out-of-bag is not stable), AUC/ROC of the data accumulated from k fold

#' Predict binary variables based on a taxa read count table by random forest
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param factor The binary variable to predict, e.g., the key variable, case/control.
#' @param fold The number of folds; default is 5.
#' @param main The title of plot; default is NULL.
#' @param seed The seed to generate fold indices for samples; default is 2020.
#'
#' @return A list
#' \itemize{
#'   \item Print a ROC curve of predictions accumulated from the folds, e.g., on all samples.
#'   \item pred - A table summarizing the predicted probabilities and true labels for all samples.
#'   \item auc_across_fold - AUC of the ROC curves across folds.
#'   \item auc_on_all - AUC of the ROC curve on all samples (the printed).
#'}
#'
#' @export

RF_Pred <- function(TAX, factor, fold=5, main=NULL, seed=2020){

  set.seed(seed)
  ss = sample(1:fold,size=nrow(TAX),replace=T,prob=rep(1/fold, fold)) # assign fold for samples

  temp = as.matrix(TAX)
  rownames(temp) = NULL
  colnames(temp) = NULL

  # do prediction across folds
  record_tab = NULL
  record_auc = NULL
  for (ii in 1:fold){
    trainy = factor[ss!=ii]
    testy = factor[ss==ii]

    trainx = temp[ss!=ii, ]
    testx = temp[ss==ii, ]

    rf_classifier = randomForest(trainy ~ ., data=trainx, importance=TRUE)

    prediction_for_roc_curve = predict(rf_classifier,testx,type="prob")
    pred = prediction(prediction_for_roc_curve[,2],testy)
    record_tab = rbind(record_tab, data.frame(prob=prediction_for_roc_curve[,2],testy))

    auc.perf = performance(pred, measure = "auc")
    record_auc[ii] = auc.perf@y.values[[1]]
  }

  # ROC on the accumulated prediction, e.g., on all samples
  pred = prediction(record_tab[,1],record_tab[,2])
  perf = performance(pred, "tpr", "fpr")

  auc.perf = performance(pred, measure = "auc")
  auc_value = auc.perf@y.values[[1]]

  plot(perf, main=main, col=2)
  abline(coef = c(0,1))
  text(0.85, 0.1, paste0("AUC=", round(auc_value, digits=2)))

  return(list(pred=record_tab, auc_across_fold=record_auc, auc_on_all=auc_value))

}


# Predict continuous key variable using random forest, k-fold cross-validation (out-of-bag is not stable), rmse from k fold

#' Predict continuous variables based on a taxa read count table by random forest
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param variable The continuous variable to predict.
#' @param fold The number of folds; default is 5.
#' @param main The title of plot; default is NULL.
#' @param seed The seed to generate fold indices for samples; default is 2020.
#'
#' @return A list
#' \itemize{
#'   \item Print a boxplot of RMSEs across folds.
#'   \item pred - A table summarizing the predicted and true values for all samples.
#'   \item rmse_across_fold - RMSEs across folds.
#'}
#'
#' @export

RF_Pred_Regression <- function(TAX, variable, fold=5, main=NULL, seed=2020){

  set.seed(seed)
  ss = sample(1:fold,size=nrow(TAX),replace=T,prob=rep(1/fold, fold)) # assign fold for samples

  temp = as.matrix(TAX)
  rownames(temp) = NULL
  colnames(temp) = NULL

  # do prediction across folds
  record_tab = NULL
  record_rmse = NULL
  for (ii in 1:fold){
    trainy = variable[ss!=ii]
    testy = variable[ss==ii]

    trainx = temp[ss!=ii, ]
    testx = temp[ss==ii, ]

    rf_regression = randomForest(trainy ~ ., data=trainx, importance=TRUE)

    prediction_for_rmse <- predict(rf_regression,testx)
    record_tab = rbind(record_tab, data.frame(pred=prediction_for_rmse,testy))

    record_rmse[ii] = sqrt(mean( (prediction_for_rmse - testy)^2 ))
  }

  boxplot(record_rmse, main=main)

  return(list(pred=record_tab, rmse_across_fold=record_rmse))

}

