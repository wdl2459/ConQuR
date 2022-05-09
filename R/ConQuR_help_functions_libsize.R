
########## ConQuR-libsize: help functions ##########


### fast computation of zero-inflated conditional quantile function

# fast computation (not by exact tau)
proposed.method.fast.libsize <- function(logistic, xl, quantile, qmat, n, libsize, delta, taus, logistic_lasso, interplt){

  # estimate the probability of being positive
  if (logistic_lasso == T){
    p_hat = predict(logistic, newx=xl, type="response")
  } else p_hat = predict(logistic, newdata=xl, type="response")

  # estimate the 3(2) parts, zero, transition (optional), positive
  Taus.s = P1 = P2 = list()
  for (ii in 1:nrow(xl)){
    if (interplt == T){
      P1[[ii]] = length( which(taus < (1 - p_hat[ii])) )
      P2[[ii]] = ( taus[ which(taus >= (1 - p_hat[ii]) & taus <= (1 - p_hat[ii] + n^(-delta))) ] - 1 + p_hat[ii] ) * n^delta
      part3 = ( taus[ which(taus > (1 - p_hat[ii] + n^(-delta))) ] - 1 + p_hat[ii] ) / p_hat[ii]

      Taus.s[[ii]] = c(n^(-delta) / p_hat[ii], part3)
    } else{
      P1[[ii]] = length( which(taus <= (1 - p_hat[ii])) )
      Taus.s[[ii]] = ( taus[ which(taus > (1 - p_hat[ii])) ] - 1 + p_hat[ii] ) / p_hat[ii]
    }
  }

  # estimate and map the positive part, construct conditional quantile function with 3(2) parts
  quant = matrix(0, ncol=length(taus), nrow=nrow(xl))
  for (ii in 1:nrow(xl)){
    fittedy = as.numeric(c(1, qmat[ii, ])) %*% as.matrix(quantile)

    ###### libsize specific ######
    fittedy = exp(fittedy) # transform back to relative abundance

    if (length(Taus.s[[ii]]) > 0 & Taus.s[[ii]][1] < 1){
      location = unlist(lapply(Taus.s[[ii]], locate.tau, taus))
      fit = fittedy[location]
      if (any(location == 0)){
        fit = c(rep(0, length(which(location == 0))), fit)
      }

      ###### libsize specific ######
      if (interplt == T){
        temp_qf = c( rep(0, P1[[ii]]), fit[1]*P2[[ii]], fit[-1] )
      } else{
        temp_qf = c( rep(0, P1[[ii]]), fit )
      }
      # temp_qf_sf = stepfun(taus, c(0, temp_qf))
      # quant[ii, ] = rearrange(temp_qf_sf)(taus)
      quant[ii, ] = temp_qf
    } else quant[ii, ] = rep(0, length(taus))
  }

  return(quant)

}

### simple quantile-quantile matching

simple_QQ_libsize <- function(y, batchid, batch_ref, taus){

  # info about batchid
  tab_batch = table(batchid)
  batchN = length(tab_batch)
  name_batch = names(tab_batch)

  # empirical quantile functions of reference and other levels of batch
  ziqr.fast = vector(mode="list", length=batchN)
  for (mm in 1:batchN){
    ziqr.fast[[mm]] = quantile(y[which(batchid == name_batch[mm])], probs=taus, type=1)
  }
  ziqr.fast.correct = ziqr.fast[[1]]

  # quantile-quantile matching
  match_count = unlist( lapply(1:length(y), function(kk){
    if (batchid[kk] == batch_ref){
      value = y[kk]
    } else{
      ref = ziqr.fast[[which( name_batch == batchid[kk] )]]

      loc = which(ref == y[kk])
      value = mean( ziqr.fast.correct[loc] ) # if multiple quantiles equal to obs, take the mean of corrected

      if (length(loc) == 0){
        loc_temp = which(ref < y[kk])
        if (length(loc_temp) == 0) loc_temp = 1
        value = ziqr.fast.correct[max(loc_temp)] # if nothing equal to obs, take the max of those smaller than obs: quantile function is left-continuous
      }
    }

    value
  }) )

  return(match_count)

}


### ConQuR for each taxon

ConQuR_each_libsize <- function(y, X, X_span, X_correct, X_span_correct, batch_ref=batch_ref,
                                libsize,
                                delta, taus, logistic_lasso, quantile_type, lambda_quantile, interplt, logistic_thres=8){

  to_skip <<- F


  ### logistic regression

  sy = 1*(y > 0)
  if ( sum(sy) <= 1 ){ y_new = y; return(y_new) } # keep original count when there is no or only 1 meaningful observation

  if ( sum( table(sy) <= logistic_thres ) | length(table(sy)) == 1 | logistic_lasso == F ){

    # choose or force (not balanced or too few in either categories) to run standard logistic
    logistic_lasso = F

    # original design matrix + corrected design matrix with non-ref levels of batch == 0
    X_logistic = cbind(X, libsize)
    X_logistic_correct = cbind(X_correct, libsize)

    logistic_fit = glm(sy ~ ., family = "binomial", data = X_logistic)

  } else{

    # standardize original design matrix for lasso logistic
    X_logisic_span = cbind(X_span, libsize)
    x_logistic = scale(X_logisic_span, center=T, scale=T)

    # lasso logistic
    logistic_cv = cv.glmnet(x_logistic, sy, family = "binomial", alpha = 1, lambda = NULL)
    logistic_fit = glmnet(x_logistic, sy, family = "binomial", alpha = 1, lambda = logistic_cv$lambda.min)

    # original design matrix (standardized) + corrected design matrix (standardized) with non-ref levels of batch == (standardized) 0
    X_logistic = x_logistic

    X_logistic_correct = x_logistic
    batch_null_values = apply( X_logistic_correct[, grep( "batchid", colnames(X_logistic_correct) ), drop=F], 2, min )
    X_logistic_correct[,  grep( "batchid", colnames(X_logistic_correct) )] = t(replicate(nrow(X_logistic_correct), batch_null_values))

  }


  ### quantile regression

  index = which(y > 0)
  y_sub = y[index]
  data_sub = X[index, ]

  # arrange design matrices according to quantile regression type
  if (quantile_type == "standard"){

    # check contrasts, i.e., distributions of factor covariates + restrict to useful covariates space
    todel = NULL
    for (kk in 1:ncol(data_sub)){
      if (is.factor(data_sub[, kk]) & sum(table(data_sub[, kk]) > 0) == 1){ todel = c(todel, kk) }
    }
    if (!is.null(todel)){ data_sub = data_sub[, -todel, drop=F] }

    # check whether data_sub becomes empty
    if (ncol(data_sub) == 0){
      y_new = simple_QQ_libsize(y=y, batchid=batchid, batch_ref=batch_ref, taus=taus)
      return(y_new)
    } # use simple match when data_sub becomes empty

    # check singularity, i.e, whether p >= n
    cov_num = ncol( model.matrix(~., data_sub) )
    if (cov_num >= length(index)){ y_new = y; return(y_new) } # keep original count when p >= n

    standard_name = data.frame(standname = c("(Intercept)", colnames(X_span)))

    X_quantile = X_span
    X_quantile_correct = X_span_correct

  } else if (quantile_type == "composite" | quantile_type == "lasso"){

    # relevel factors in data_sub and X according to the remaining levels
    X_consistent_to_data_sub = X

    for (kk in 1:ncol(data_sub)){
      if ( is.factor(data_sub[, kk]) & table(data_sub[, kk])[1] == 0 ){
        ref_new = names(table(data_sub[, kk]))[which( table(data_sub[, kk]) != 0 )[1]]
        data_sub[, kk] = relevel(data_sub[, kk], ref=ref_new)
        X_consistent_to_data_sub[, kk] = relevel(X_consistent_to_data_sub[, kk], ref=ref_new)
      }
    }

    # sub design matrix for composite or lasso quantile
    x = model.matrix(~., data_sub)[, -1]

    # subtract and divide the design matrix according to sub design matrix
    x_entire = model.matrix(~., X_consistent_to_data_sub)[, -1]
    x_entire = sweep(x_entire, 2, apply(x, 2, mean), '-')
    x_quantile_entire = sweep(x_entire, 2, apply(x, 2, sd), '/')
    x_quantile_entire[, which( is.na( colSums(x_quantile_entire))  ) ] = 0 # fill 0 in NA just to avoid problem, values won't matter

    # standardize sub design matrix + restrict to useful covariates space
    x = scale(x, center=T, scale=T)
    x_quantile = x[, colSums(!is.na(x)) > 0, drop=F]

    # check whether x_quantile becomes empty
    if (ncol(x_quantile) == 0){
      y_new = simple_QQ_libsize(y=y, batchid=batchid, batch_ref=batch_ref, taus=taus)
      return(y_new)
    } # use simple match when x_quantile becomes empty

    # original design matrix (standardized according to sub design matrix) + corrected design matrix (standardized according to sub design matrix) with non-ref levels of batch == (standardized) 0
    standard_name = data.frame(standname = c("(Intercept)", colnames(x_quantile_entire)))

    X_quantile = x_quantile_entire

    X_quantile_correct = x_quantile_entire
    batch_null_values = apply( X_quantile_correct[, grep( "batchid", colnames(X_quantile_correct) ), drop=F], 2, min )
    X_quantile_correct[, grep( "batchid", colnames(X_quantile_correct) )] = t(replicate(nrow(X_quantile_correct), batch_null_values))

  } else{
    stop("Currently, ConQuR does not support other types of quantile regression.")
  }

  # specify lambda for composite or lasso quantile
  if (quantile_type == "composite" | quantile_type == "lasso"){

    if (lambda_quantile == "2p/n"){
      lambda_q = 2*ncol(x_quantile) / nrow(x_quantile)
    } else if (lambda_quantile == "2p/logn"){
      lambda_q = 2*ncol(x_quantile) / log(nrow(x_quantile))
    } else{
      stop("Please specify lambda for composite or lasso quantile regression, 2p/n or 2p/logn.")
    }

  }

  # fit quantile model
  tryCatch({

    ###### libsize specific ######
    y_sub = log(y_sub) # model log relative abundance <=> model log count with library size as an offset

    # quantile regression
    if (quantile_type == "standard"){
      quantile_fit = rq(y_sub ~ ., data=data_sub, tau=taus)
    } else if (quantile_type == "composite"){
      quantile_fit = cqr.fit.lasso(x_quantile, y=y_sub, tau=taus, lambda=lambda_q, method="cd")
    } else if (quantile_type == "lasso"){
      quantile_fit = rq(y_sub ~ ., data=data.frame(x_quantile), tau=taus, lambda=lambda_q, method="lasso")
    }

  }, error=function(e){ to_skip <<- T })
  if (to_skip == T){ y_new = y; return(y_new) } # keep original count when singularity occurs


  ### construct conditional quantile functions

  # arrange the coefficients
  if (quantile_type == "composite"){
    mat = matrix( rep(quantile_fit$beta, length(taus)), ncol=length(taus), byrow=F )
    mat = rbind(quantile_fit$b, mat)
    rownames(mat) = c("(Intercept)", colnames(x_quantile))
    colnames(mat) = paste("tau=", taus)

    quant_coef = data.frame( coefname=rownames(mat), mat )
  } else if (quantile_type == "standard" | quantile_type == "lasso"){
    quant_coef = data.frame( coefname=rownames(quantile_fit$coefficients), quantile_fit$coefficients )
  }

  quant_coef = standard_name %>%
    left_join(quant_coef, by=c("standname" = "coefname"))
  quant_coef = quant_coef[, -1]
  quant_coef[is.na(quant_coef)] = 0

  # piecewise estimation of conditional quantile functions
  ziqr.fast = proposed.method.fast.libsize(logistic=logistic_fit, xl=X_logistic,
                                           quantile=quant_coef, qmat=X_quantile,
                                           n=length(y), libsize=libsize,
                                           delta=delta, taus=taus, logistic_lasso=logistic_lasso, interplt=interplt)

  ziqr.fast.correct = proposed.method.fast.libsize(logistic=logistic_fit, xl=X_logistic_correct,
                                                   quantile=quant_coef, qmat=X_quantile_correct,
                                                   n=length(y), libsize=libsize,
                                                   delta=delta, taus=taus, logistic_lasso=logistic_lasso, interplt=interplt)


  ### quantile-quantile matching

  y_new = unlist( lapply(1:length(y), function(kk){
    loc = which(ziqr.fast[kk, ] == y[kk])
    value = mean( ziqr.fast.correct[kk, loc] ) # if multiple quantiles equal to obs, take the mean of corrected

    if (length(loc) == 0){
      loc_temp = which(ziqr.fast[kk, ] < y[kk])
      if (length(loc_temp) == 0) loc_temp = 1
      value = ziqr.fast.correct[kk, max(loc_temp)] # if nothing equal to obs, take the max of those smaller than obs: quantile function is left-continuous
    }

    value
  }) )


  return(y_new)

}

