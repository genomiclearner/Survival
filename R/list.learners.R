#'@title List Available Base Learners
#'@name list.learners
#'@description This function lists all base learners provided in the package.
#'@details lm1: linear regression

#'@details glm1: generalized linear models

#'@details glmnet1: Does k-fold cross-validation to chose best alpha and lambda for generalized linear models via penalized maximum likelihood.

#'@details glmnet.lasso: LASSO, lambda is chose by k-fold cross-validation for glmnet

#'@details glmnet.ridge: Ridge regression, lambda is chose by k-fold cross-validation for glmnet

#'@details rpart1: regression tree

#'@details lda1: linear discriminant analysis

#'@details qda1: quadratic discriminant analysis

#'@details KNN1: k-nearest neighbour classification, k is chose by cross-validation

#'@details svm1: support vector machine

#'@details surv: parametric survival regression

#'@details ela1: AFT model with Elasticnet, weights is the weight of survival time

#'@details lasso1: AFT model with Lasso, weights is the weight of survival time



#'@export list.learners
#'@return The name of all base learners provided in the package
#'@examples
#'list.learners()


list.learners <- function(){

  learners <- c("lm1", "glm1", "glmnet1", "glmnet.lasso", "glmnet.ridge", "rpart1",
                "lda1", "qda1", "KNN1", "svm1","ela1","surv","lasso1")
  print("Models that can be chosen are:")
  print(learners)

}
