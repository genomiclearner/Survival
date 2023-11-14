#'@importFrom utils head
#'@importFrom stats coef fitted glm kmeans lm model.matrix predict quantile time
#'@importFrom caret trainControl train
#'@importFrom glmnet glmnet cv.glmnet
#'@importFrom rpart rpart prune
#'@importFrom MASS lda qda
#'@importFrom e1071 svm
#'@importFrom class knn knn.cv
#'@importFrom survival Surv survreg

#'@export modify.parameter
modify.parameter <- function(FUN, ...) {
  if (!is.function(FUN)) {
    stop("FUN must be a valid function")
  }
  .FUN <- FUN
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  .FUN
}

#'@export glmnet1
glmnet1 <- function(xmat, ymat,family,alpha=seq(0,10,by=2)/10, ...) {
  tmp0 <- data.frame(yy=ymat, xmat)
  if(family=="binomial") ymat <- factor(ymat)
  foldid <- createFolds(ymat, k=5, list=FALSE)
  model <-  cv.glmnet2(xmat, ymat, alpha=alpha, foldid=foldid, family=family, ...)
  coef.mat <- as.numeric(coef(model,s="lambda.1se"))
  y.fitted <- predict(model, xmat, s="lambda.1se", type="response")
  predFun <- function(model,xnew){
    predict(model, newx=as.matrix(xnew), s="lambda.1se", type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export glmnet.lasso
glmnet.lasso <- modify.parameter (glmnet1, alpha=1)
#'@export glmnet.ridge
glmnet.ridge <- modify.parameter (glmnet1, alpha=0)

#'@export surv
surv <- function(xmat, ymat, family, dist1, ...) {
  tmp0 <- data.frame(yy=Surv(ymat[,1],ymat[,2]), xmat)
  model <- survreg(yy~., data=tmp0, dist=dist1)
  y.fitted <- fitted(model)
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export lm1
lm1 <- function(xmat, ymat, family,weights, ...) {
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- lm(yy~., data=tmp0,weights=weights, ...)
  y.fitted <- fitted(model)
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew), type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export ela1
ela1 <- function(xmat,ymat,weights, ...){
  ymat<-ymat[,1]

  cv_5 = trainControl(method = "cv", number = 5)
  tmp0 <- data.frame(yy=ymat, xmat)

  myfit_elnet<-train(log(yy) ~ .,data = tmp0,method = "glmnet",trControl = cv_5)

  myid<-order(myfit_elnet$result$RMSE)[1]
  mybest<-myfit_elnet$results[myid, ]

  xmat2<- model.matrix(~.,xmat)[,-1]
  model <- glmnet(xmat2, as.numeric(log(ymat)),  alpha=mybest$alpha, lambda=mybest$lambda, weights = weights)

  y.fitted <- exp(predict(model, s=mybest$lambda, newx=xmat2))

  predFun<-function(model,xnew){
    xnew<- model.matrix(~.,xnew)[,-1]
    mybest<-model$mybest
    exp(predict(model, s=mybest$lambda,newx=xnew))
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun, mybest=mybest))
}

#'@export lasso1
lasso1<-function(xmat,ymat,weights,alpha=1,lambda=10^seq(10, -2, length=100), ...){

  xmat2<- model.matrix(~.,xmat)[,-1]
  lasso.cvout<-cv.glmnet(xmat2, as.numeric(log(ymat[,1])),weights = weights, alpha = 1, lambda = lambda)
  mybest<-lasso.cvout$lambda.min

  model<-glmnet(xmat2, as.numeric(log(ymat[,1])), weights = weights, alpha = 1, lambda = mybest)

  y.fitted<-exp(predict(model,s=mybest,xnew = xmat2))

  predFun<-function(model,xnew){
    xmat3<-model.matrix(~.,xnew)[,-1]
    exp(predict(model,s=model$mybest,xnew = xmat3))
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun, mybest=mybest))

}

#'@export glm1
glm1 <- function(xmat, ymat, family, ...) {
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- glm(yy~., data=tmp0, family=family, ...)
  y.fitted <- fitted(model)
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew), type="response")
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}


#'@export rpart1
rpart1 <- function(xmat, ymat, family, ...){
  tmp0 <- data.frame(yy=ymat, xmat)
  fit0 <- rpart(yy~., data=tmp0, ...)
  model <- prune(fit0, cp=fit0$cptable[which.min(fit0$cptable[,"xerror"]),"CP"])
  y.fitted <- predict(model, newdata=data.frame(xmat))
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export lda1
lda1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- lda(yy~., data=tmp0, ...)
  y.fitted <- predict(model, newdata=data.frame(xmat))$posterio[,"1"]
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))$posterio[,"1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export qda1
qda1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- qda(yy~., data=tmp0, ...)
  y.fitted <- predict(model, newdata=data.frame(xmat))$posterio[,"1"]
  predFun <- function(model,xnew){
    predict(model, newdata=data.frame(xnew))$posterio[,"1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export KNN1
KNN1 <- function(xmat, ymat, family, ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  # Select k applying cross-validation
  knn_acc <- rep(NA, sqrt(length(ymat)))
  for (i in 1:sqrt(length(ymat))){
    # print(i)
    results<-knn.cv(train=xmat, cl=ymat, k = i, prob = TRUE)
    knn_acc[i]<-(table(results,ymat)[1,1]+table(results,ymat)[2,2])/(sum(table(results,ymat)))
  }
  k <- which.max(knn_acc)
  model <- list(model=knn.cv(train=xmat, cl=ymat, k = k, prob = TRUE),
                  train=xmat,cl=ymat,k=k)
  y.fitted <- as.numeric(as.character(model$model)) * attributes(model$model)$prob +
    (1-as.numeric(as.character(model$model))) * (1-attributes(model$model)$prob)
  predFun <- function(model,xnew){
    knn.fit <- knn(train=model$train, test=xnew, cl=model$cl, k=model$k, prob = TRUE)
    as.numeric(as.character(knn.fit)) * attributes(knn.fit)$prob +
      (1-as.numeric(as.character(knn.fit))) * (1-attributes(knn.fit)$prob)
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}

#'@export svm1
svm1 <- function(xmat, ymat, family, kernel = "linear", ...){
  if(family=="binomial") ymat <- factor(ymat)
  tmp0 <- data.frame(yy=ymat, xmat)
  model <- svm(yy~.,data=tmp0, kernel = kernel, probability = TRUE, ...)
  pred <- predict(model, xmat, probability=TRUE)
  y.fitted <- attr(pred, "probabilities")[, "1"]
  predFun <- function(model,xnew){
    pred <- predict(model, xnew, probability=TRUE)
    attr(pred, "probabilities")[, "1"]
  }
  return(list(model=model,y.fitted=y.fitted, predFun=predFun))
}
