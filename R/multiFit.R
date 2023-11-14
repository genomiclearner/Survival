#'@title Fit models on multiple outcomes.
#'@param xmat Matrix of predictors, each row is an observation vector.
#'@param ymat Matrix of outcomes.
#' Quantitative for family = \strong{"gaussian"}.
#' A factor of two levels for family = \strong{"binomial"}.
#' A survival object for family = \strong{"survival"}.
#'@param xmat.list The user defines a numeric vector to specify the number of subset columns of the xmat.
#'@param method Method for fitting models.
#'It can be one base learner function for all outcomes or a list of base learner functions for each outcome.
#'The list of all base learners can be obtained by \code{list.learners()}.
#'@param family Response type for each response.
#'If all response variable are within the same family it can be \strong{"gaussian"}, \strong{"binomial"} or \strong{"survival"},
#'otherwise it is a vector with elements \strong{"gaussian"}, \strong{"binomial"} and \strong{"survival"} to indicate each response family.
#'@param dist1 Assumed distribution for response variable.
#' If the argument is a character string, then it is assumed to name an element from survreg.distributions.
#' These include \strong{"weibull"}, \strong{"exponential"}, \strong{"gaussian"}, \strong{"logistic"},\strong{"lognormal"} and \strong{"loglogistic"}. Otherwise, it is assumed to be a user defined list conforming to the format described in \code{"survreg.distributions"}.
#'@param weights Weight for response.
#'@return It returns a multiFit object. It is a list of 6 parameters containing information about the fitted models and fitted values for each outcome.
#'@description This function fit individual models to predict each outcome separately.
#'@export multiFit



multiFit <- function(xmat, ymat,  xmat.list,
                     method,family, dist1=dist1,weights=weights)
{

  nx <- ncol(xmat)
  ny <- length(family)

  # check family input
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial or survival")
    }
    if (family == "gaussian") {
      family = rep("gaussian", ny)
    } else if (family == "binomial") {
      family = rep("binomial", ny)
    }else if (family== "survival"){
      family = rep("survival", ny)
    }
  }


  # if (length(family) != ny) {
  #   stop("length of family must be consistent with response")
  # }
  #
  # if (sum(family %in% c("gaussian", "binomial")) != ny) {
  #   stop("each family must be gaussian or binomial")
  # }
  #
  # check family method consistency
  if (length(method) == 1) {
    method <- rep(list(method),ny)
  }
  # if (length(method) != ny) {
  #   stop("length of method.step1 must be 1 or the same as response column")
  # }
  # for (ii in 1:ny) {
  #   if (!check.match(family[ii], FUN=method[[ii]])) {
  #     stop("method.step1 must be consistent with response category")
  #   }
  # }

  # initialize matrix to save fitted results

  y.fitted <- matrix(NA, nrow=nrow(xmat), ncol=ny)
  models <- vector("list", ny)
  colnames(y.fitted) <- names(models) <- colnames(ymat)
  fit <- vector("list", ny)
  mybest <- vector("list", ny)
  # colnames of predictor matrix are needed to avoid problem when calling function "predict"
  colnames(xmat) <- paste0("X", 1:nx)
  xmat0 <- xmat

  for (kk in 1:ny){
    xmat <- xmat0
    xmat <- xmat[,xmat.list[[kk]]]
    xmat <- as.data.frame(xmat)
    names(xmat) <- names(xmat0)[xmat.list[[kk]]]
    if(is.null(dist1)){
      if(family[kk]%in%'survival'){
        fit[[kk]] <- method[[kk]](xmat, ymat[[kk]],family=family[kk],weights = weights[,kk])}
      else fit[[kk]] <- method[[kk]](xmat, ymat[,kk],family=family[kk],weights = weights[,kk])
      models[[kk]] <- fit[[kk]]$model
      y.fitted[,kk] <- fit[[kk]]$y.fitted
      mybest[[kk]] <- fit[[kk]]$mybest## elasticnet best params
    }else{
      if(family[kk]%in%'survival'){
        fit[[kk]] <- method[[kk]](xmat, ymat[[kk]], family=family[kk],weights = weights[,kk],dist1=dist1)}
      else fit[[kk]] <- method[[kk]](xmat, ymat[,kk], family=family[kk],weights = weights[,kk],dist1=dist1)
      models[[kk]] <- fit[[kk]]$model
      y.fitted[,kk] <- fit[[kk]]$y.fitted
      mybest[[kk]] <- fit[[kk]]$mybest ## elasticnet best params
    }


  }


  multiFit.fits <- list(fit=fit, xmat.list=xmat.list,
                        y.fitted=y.fitted,
                        model=models,
                        mybest=mybest,#elasticnet params
                        method=method,
                        family=family)
  class(multiFit.fits) <- "multiFit"
  return(multiFit.fits)
}

#'@method predict multiFit
#'@export predict.multiFit
#'@export
predict.multiFit <- function(object, newdata, xmat.list,...)
{
  ny <- length(object$model)
  newdata <- as.data.frame(newdata)
  pred <- matrix(NA, nrow(newdata), ny)
  mtd <- object$method

  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  newdata0 <- newdata
  for(ii in 1:ny)
  {
    newdata <- newdata0
    newdata <- as.data.frame(newdata[,object$xmat.list[[ii]]])
    names(newdata) <- names(newdata0)[xmat.list[[ii]]]
    model <- object$model[[ii]]
    pred[,ii] <- object$fit[[ii]]$predFun(model, newdata)

    if(object$family[ii]=="binomial")
    { # predicted are probability should be within [0,1]
      pred[,ii][pred[,ii]>1] <- 1
      pred[,ii][pred[,ii]<0] <- 0
    }
  }
  return(pred)
}



