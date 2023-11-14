#'@title Area Under Curve
#'@description The AUC function calculates the numeric value of area under the ROC curve \(AUC\) with the trapezoidal rule and optionally plots the  ROC curve
#'@details The ROC curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at different threshold settings.
#'By default the total area under the curve is computed, but a truncated AUC statistics can be specified with the cutoff argument.
#'It specifies the bounds of FPR. The common choice of cutoff can be 1 (i.e. no truncate) or 0.2 (i.e. specificity > 0.8)
#'@param prob A numeric vector of predicted probability
#'@param outcome A numeric vector of observed binary outcome
#'@param cutoff Number between 0 and 1 to specify where threshold of ROC curve should be truncated. The default value is 1 \(no truncation\)
#'@param ROC.plot Logical. Whether or not to plot ROC curve
#'@export AUC
#'@return The value of the area under the curve.
#'@importFrom graphics abline
#'@examples
#' set.seed(1)
#' # simulate predictors
#' x1 <- rnorm(200)
#' x2 <- rnorm(200)
#' # simulate outcome
#' pr <- 1/(1+exp(-(3 * x1 + 2 * x2 + 1)))
#' y <- rbinom(200, 1, pr)
#' df <- data.frame(y = y,x1 = x1, x2 = x2)
#' # fit logistic regression model on the first 100 observation
#' lg.model <- glm(y ~ x1 + x2, data = df[1 : 100, ], family="binomial")
#' # predict outcome for the last 100 observation
#' prob <- predict(lg.model, df[101:200, c("x1", "x2")], type = "response")
#' # calculate AUC and plot thr ROC Curve
#' AUC(prob, y[101:200], ROC=TRUE)
#' # calculate AUC and plot thr ROC Curve with cutoff
#' AUC(prob, y[101:200], cutoff=0.2, ROC=TRUE)


AUC <- function(prob, outcome, cutoff = 1, ROC.plot = FALSE)
{
	NN <- length(outcome)
	if(length(prob)!=NN) stop("prob and binary should be the same length")
	if((max(prob)>1) | (min(prob)<0)) stop("prob values should be in [0,1]")

	# remove missing values, which cannot contribute to AUC calculation
	idx.retain <- (!is.na(outcome)) & (!is.na(prob))
	prob <- prob[idx.retain]
	outcome  <-  outcome[idx.retain]

	# sort binary outcome by order of predicted probability
	ord <- order(prob, decreasing=T)
	prob <- prob[ord]
	outcome <- outcome[ord]

	CP <- sum(outcome)    # condition positive
	CN <- NN-sum(outcome) # condition negative
	TP <- cumsum(outcome) # true positive
	FP <- (1:NN) - cumsum(outcome) # false Positve
	TPR <- TP/CP		 # True positive rate (sensitivity)
	FPR <- FP/CN		 # False positive rate (1- specificity)
	if(ROC.plot) {plot(FPR, TPR); abline(v=cutoff)}
	auc <- trapezoid(FPR, TPR, cutoff)

	return(auc)
}

trapezoid<- function(FPR,TPR,cc)
{
  ordFPR=order(FPR)
  FPR=FPR[ordFPR]
  TPR=TPR[ordFPR]
  if (max(FPR)<cc)
  {
    FPRcut=FPR
    TPRcut=TPR
  }else if(min(FPR)>cc)
  {
    stop("threshold smaller than smallest allowed value")
  }else
  {
    FPRlcc=sum(FPR<cc)
    FPR1=FPR[FPRlcc]; FPR2=FPR[FPRlcc+1]
    TPR1=FPR[FPRlcc]; TPR2=TPR[FPRlcc+1]
    TPRcc=( TPR1*(cc-FPR1) + TPR2*(FPR2-cc) ) / (FPR2-FPR1)
    FPRcut=c(FPR[1:FPRlcc],cc)
    TPRcut=c(TPR[1:FPRlcc],TPRcc)
  }
  return(sum( diff(FPRcut) * ( TPRcut[-1] + head(TPRcut,-1) ) )/2)
}

