\name{predict.MTPS}
\alias{predict.MTPS}
\title{Make predictions from a "MTPS" model}
\description{
This function makes predictions from a revised stacking model.}
\usage{
\method{predict}{MTPS}(object, newdata, ...)
}
\arguments{
  \item{object}{A fitted object from \code{"MTPS"}}

  \item{newdata}{Matrix of new predictors at which predictions are to be made}

  \item{...}{additional arguments affecting the predictions produced}
}
\value{
  The predicted value from new predictors.
}

\examples{
data("HIV")
set.seed(1)
xmat <- as.matrix(XX)
ymat <- as.matrix(YY)
id <- createFolds(rowMeans(XX), k=5, list=FALSE)
training.id <- id != 1
y.train <- ymat[training.id, ]
y.test  <- ymat[!training.id, ]
x.train <- xmat[training.id, ]
x.test  <- xmat[!training.id, ]
# Cross-Validation Residual Stacking
fit.rs <- MTPS(xmat = x.train, ymat = y.train,
  family = "gaussian",cv = FALSE, residual = TRUE,
  method.step1 = rpart1, method.step2 = lm1)
pred.rs <- predict(fit.rs, x.test)

# Residual Stacking for Survival Analysis
  set.seed(1)
  data("simdat_mtps")

 # prepare training and test set
 id.train <- sample(1:100, 80)
 xmat.train <- xmat[id.train,]
 xmat.test <- xmat[-id.train,]
 ymat.train <- cbind(list(survival::Surv(ymat[id.train,"time01"],ymat[id.train,"status01"])),
 list(survival::Surv(ymat[id.train,"time02"],ymat[id.train,"status02"])))
 # Produce the Kaplan-Meier estimator
 weights <- find_km_weights_mat(ymat[id.train,],num_outcome = 2)
 # fit Residual Stacking Model for Survival Data
 fit <- MTPS(xmat.train, ymat.train, family = 'survival', cv=FALSE,
 residual = TRUE, nfold=5, method.step1 = surv,
 method.step2 = lm1, dist1 = "lognormal", weights = weights)
 # predict the survival time on test set
 pre <- predict.MTPS(fit, xmat.test)



}
