#'@title Evaluation using Cross-Validation
#'@description Use cross-validation to evaluate model performance.
#'@details The most frequently used evaluation metric of survival models is the concordance index (C-index).
#'It is a measure of rank correlation between predicted risk scores. Only use uncensored observations
#'@param xmat Matrix of predictors, each row is an observation vector
#'@param ymat Matrix of outcomes.
#' Quantitative for family = \strong{"gaussian"}.
#' A factor of two levels for family = \strong{"binomial"}.
#' A survival object for family = \strong{"survival"}.
#'@param family Response type for each response.
#'If all response variable are within the same family it can be \strong{"gaussian"}, \strong{"binomial"} or \strong{"survival"},
#'otherwise it is a vector with elements \strong{"gaussian"}, \strong{"binomial"} and \strong{"survival"} to indicate each response family.
#'@param nfolds Integer, number of folds for Cross-Validation to evaluate the performance of stacking algorithms.
#'@param cv Logical, indicate if use Cross-Validation Stacking algorithm.
#'@param residual Logical, indicate if use Residual Stacking algorithm.
#'@param cv.stacking.nfold Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5.
#'@param method.step1 Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.
#'@param method.step2 Base Learners for fitting models in Step 2 of Stacking Algorithm. (see above)
#'@param resid.type The residual type for Residual Stacking.
#'@param resid.std Logical, whether or not use standardized residual.
#'@param dist1 Assumed distribution for survival response.
#' If the argument is a character string, then it is assumed to name an element from survreg.distributions.
#' These include \strong{"weibull"}, \strong{"exponential"}, \strong{"gaussian"}, \strong{"logistic"},\strong{"lognormal"} and \strong{"loglogistic"}.
#' Otherwise, it is assumed to be a user defined list conforming to the format described in \code{"survreg.distributions"}.
#'@param weights Weight for survival response.
#'@export cv.MTPS
#'@return It returns the mean squared error of continuous outcomes. AUC, accuracy, recall and precision for binary outcomes of predictions using cross-validation.



cv.MTPS <- function(xmat, ymat, family,
                               nfolds = 5,
                               cv = FALSE, residual = TRUE,
                               cv.stacking.nfold = 5, method.step1, method.step2,
                               resid.type=c("deviance", "pearson", "raw"),
                               resid.std=FALSE,dist1= NULL, weights=NULL)
{

  resid.type <- match.arg(resid.type)

  ny <- ncol(ymat)
  # check family input
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
  if (sum(family %in% c("gaussian", "binomial")) != ny) {
    stop("family must be gaussian or binomial or their combination")
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
  # metrics
  cindex <- which(family=="gaussian")
  bindex <- which(family=="binomial")
  # if (length(cindex) > 0) {
    metrics.ctn <- matrix(NA, nrow = 1, ncol = length(cindex))
    colnames(metrics.ctn) <- colnames(ymat[, cindex])
    rownames(metrics.ctn) <- "MSE"
    metrics.ctn <- as.data.frame(metrics.ctn)
  # }
  # if (length(bindex) > 0) {
    metrics.bin <- matrix(NA, nrow = 4, ncol = length(bindex))
    colnames(metrics.bin) <- colnames(ymat[, bindex])
    rownames(metrics.bin) <- c("AUC", "Accuracy", "Recall", "precision")
    metrics.bin <- as.data.frame(metrics.bin)
  # }

  idx.cv <- createFolds(rowMeans(xmat), k=nfolds, list=F)
  pred <- ymat; pred[!is.na(pred)] <- NA
  for (i.fold in 1:nfolds) {
    # make train and test data for i.fold-th fold
    y.train <- ymat[idx.cv!=i.fold, ]
    y.test  <- ymat[idx.cv==i.fold, ]
    x.train <- xmat[idx.cv!=i.fold, ]
    x.test  <- xmat[idx.cv==i.fold, ]
    fit <- MTPS(xmat = x.train, ymat = y.train, family = family,
                           cv = cv, residual = residual,
                           nfold = cv.stacking.nfold,
                           method.step1 = method.step1,
                           method.step2 = method.step2,
                           resid.type = resid.type, resid.std = resid.std,dist1= dist1, weights=weights)


    pred[idx.cv==i.fold, ] <- predict(fit, x.test)
  }
  # metrics
  if (length(cindex) > 0) {
    metrics.ctn["MSE",] <- apply((pred[, cindex] - ymat[, cindex])^2, 2, mean)
  }
  for (jj in bindex) {
    metrics.bin["AUC", which(jj==bindex)] <- AUC(pred[,jj], outcome=ymat[,jj])
    table <- table((pred[,jj] > 0.5) * 1, ymat[,jj])
    metrics.bin["Accuracy", which(jj==bindex)] <- (table[1,1] + table[2,2]) / sum(table)
    metrics.bin["Recall", which(jj==bindex)] <- table[2,2] / (table[2,2] + table[1,2])
    metrics.bin["precision", which(jj==bindex)] <- table[2,2] / (table[2,2] + table[2,1])
  }
  metrics <- list(continuous = metrics.ctn,
                  binary = metrics.bin)
  return(metrics)
}






