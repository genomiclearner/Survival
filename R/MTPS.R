#'@title Fit Models using Revised Stacking Algorithm.
#'@aliases MTPS
#'@description Fit a model using standard stacking algorithm or revised stacking algorithms to simultaneous predict multiple outcomes.
#'@param xmat Predictor matrix, each row is an observation vector.
#'@param ymat Responses matrix. Quantitative for family = \strong{"gaussian"}. A factor of two levels for family = \strong{"binomial"}. A survival object for family = \strong{"survival"}.
#'@param xmat.list The user defines a numeric vector to specify the number of subset columns of the xmat.
#'@param family Response type for each response. If all response variable are within the same family it can be "gaussian", "binomial" or "survival", otherwise it is a vector with elements \strong{"gaussian"}, \strong{"binomial"} and \strong{"survival"} to indicate each response family.
#'@param cv Logical, indicate if use Cross-Validation Stacking algorithm.
#'@param residual Logical, indicate if use Residual Stacking algorithm.
#'@param nfold Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5.
#'@param method.step1 Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.
#'@param method.step2 Base Learners for fitting models in Step 2 of Stacking Algorithm. (see above)
#'@param resid.type The residual type for Residual Stacking.
#'@param resid.std Logical, whether or not use standardized residual.
#'@param dist1 Assumed distribution for response variable. If the argument is a
#'
#'character string, then it is assumed to name an element from
#'
#'\code{"survreg.distributions"}. These include \strong{"weibull"},
#'
#'\strong{"exponential"}, \strong{"gaussian"}, \strong{"logistic"},
#'
#'\strong{"lognormal"} and \strong{"loglogistic"}. Otherwise, it is assumed
#'
#' to be a user defined list conforming to the format described in
#'
#' \code{"survreg.distributions"}.
#'@param weights Weight for response.
#'@return It returns a MTPS object. It is a list of 4 parameters containing
#'
#'information about step 1 and step 2 models and the revised stacking algorithm
#'
#' method.
#'@export MTPS









MTPS <- function(xmat, ymat, xmat.list, family ,
                 cv = FALSE, residual = TRUE,
                 nfold=5, method.step1, method.step2,
                 resid.type=c("deviance", "pearson", "raw"),
                 resid.std=FALSE, dist1= NULL, weights=NULL,weight2=FALSE) {


  #resid.type <- match.arg(resid.type)
  ################if 'survival', it is no use of resid.type

  # if family are different, can not neglect the total number
  ny <- ncol(ymat)
  # if(sum(family %in% c("gaussian", "binomial", "survival"))==1&&sum(family %in% c("survival"))==1){
  #   ny <- ny/2
  # }

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

  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }
  if (sum(family %in% c("gaussian", "binomial", "survival")) != ny) {
    stop("family must be gaussian or binomial or survival or their combination")
  }

  # check family method consistency
  if (length(method.step1) == 1) {
    method.step1 <- rep(list(method.step1),ny)
  }
  if (length(method.step2) == 1) {
    method.step2 <- rep(list(method.step2),ny)
  }
  if (length(method.step1) != ny) {
    stop("length of method.step1 must be 1 or the same as response column")
  }
  if (length(method.step2) != ny) {
    stop("length of method.step2 must be 1 or the same as response column")
  }
  # for (ii in 1:ny) {
  #   if (!check.match(family[ii], FUN=method.step1[[ii]])) {
  #     stop("method.step1 must be consistent with response category")
  #   }
  # }
  # if (!residual) {
  #   for (ii in 1:ny) {
  #     if (!check.match(family[ii], FUN=method.step2[[ii]])) {
  #       stop("method.step2 must be consistent with response category")
  #     }
  #   }
  # } else {
  #   for (ii in 1:ny) {
  #     if (!check.match("gaussian", FUN=method.step2[[ii]])) {
  #       stop("residual stacking does not allow binary outcome model in second step")
  #     }
  #   }
  # }

  # step 1
  if (cv) {

    fit1 <-  cv.multiFit(xmat=xmat, ymat=ymat, nfold=nfold,
                             method=method.step1,
                             family=family, dist1=dist1,weights = weights)
  } else {
    fit1 <-  multiFit(xmat=xmat, ymat=ymat, xmat.list,
                          method=method.step1,
                          family=family, dist1=dist1,weights = weights)
    #  fit <- list(fit1, fit2, fit3, fit4, fit5, fit6)

  }
  pred1 <- fit1$y.fitted


  # step 2
  if(sum(family%in%c('survival'))==ny) {
    ymat.step2 <- matrix(NA, nrow=nrow(xmat), ncol=ny)
    for(tt in 1:ny){
      ymat.step2[,tt] <- ymat[[tt]][,1]
    }
    family <- rep("gaussian", ny)
    pred1 <- log(pred1)
    ymat.step2 <- log(ymat.step2)}else ymat.step2 <- ymat
  #
  if(weight2){
    w2 <- cbind( ymat[[1]][,2],(ymat.step2-pred1)[,1], ymat[[2]][,2],(ymat.step2-pred1)[,2])
    w2 <- as.data.frame(w2)
    names(w2) <- c("statusph01","ageph01","statusph02","ageph02")
    weight1 <- find_km_weights_mat(w2)
  }else {
    weight1 <- NULL
  }

  if (residual) {
    fit2 <- rs.multiFit(yhat=pred1,ymat=ymat.step2,xmat=xmat,
                        family=family,
                        resid.type=resid.type,
                        resid.std= resid.std,
                        method=method.step2,dist1=NULL,weights=weight1) #no weight in step 2
  } else {

    fit2 <-  multiFit(xmat=pred1, ymat=ymat.step2, xmat.list,
                          method=method.step2,family=family,dist1=NULL,weights=weight1)
  }
  #}
  fit <- list(fit1 = fit1, fit2 = fit2,
              cv=cv, residual=residual)
  class(fit) <- "MTPS"
  return(fit)
}


#'@method predict MTPS
#'@export predict.MTPS
#'@export
predict.MTPS <- function(object, newdata, xmat.list,...) {

  if (object$cv) {
    pred1 <- predict(object$fit1, newdata)
  } else {
    pred1 <- predict(object$fit1, newdata, xmat.list)
  }

  # reverse the log-transformation
  if(sum(object$fit$family%in%c('survival'))==length(object$fit$family)) pred1 <- log(pred1)


  if (object$residual) {
    pred2 <- predict(object$fit2, pred1, object$fit1$mybest)
  } else {
    pred2 <- predict(object$fit2, pred1, object$fit1$mybest)
  }
  # reverse the log-transformation
  if(sum(object$fit$family%in%c('survival'))==length(object$fit$family)) pred2 <- exp(pred2)

  return(pred2)

}




