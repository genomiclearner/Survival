
cv.multiFit <- function(xmat, ymat,
                        nfold=5,
                        method,
                        family=family, dist1=dist1,weights= weights)
{

  ny <- ncol(ymat)
  nx <- ncol(xmat)

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
  if (length(method) == 1) {
    method <- rep(list(method),ny)
  }
  if (length(method) != ny) {
    stop("length of method.step1 must be 1 or the same as response column")
  }
  # for (ii in 1:ny) {
  #   if (!check.match(family[ii], FUN=method[[ii]])) {
  #     stop("method.step1 must be consistent with response category")
  #   }
  # }

  if(sum(family %in% c("survival"))==0){
    ymat.stand <- scale(ymat)
    km1 <- kmeans(ymat.stand, 10, nstart=100)
    # km1 <- kmeans(ymat, 10, nstart=100)
    idx.cv <- createFolds(factor(km1$cluster), k=nfold, list=FALSE)

    # prepare results to be returned
    y.fitted <- ymat; y.fitted[!is.na(y.fitted)] <- NA
    models <- vector("list", nfold)
    names(models) <- paste0("fold",1:nfold)
    fits <- vector("list", nfold)}else{
    ymat.time <- matrix(NA,ncol=ny, nrow=nrow(xmat))
    ymat.status <- matrix(NA,ncol=ny, nrow=nrow(xmat))
    for(i in 1:ny){
      ymat.time[,i] <- ymat[[i]][,1]
      ymat.status[,i] <- ymat[[i]][,2]
    }
    km1 <- kmeans(ymat.status, nfold, nstart=100)
    idx.cv <- createFolds(factor(km1$cluster), k=nfold, list=FALSE)

    # prepare results to be returned
    y.fitted <- ymat.time; y.fitted[!is.na(y.fitted)] <- NA
    models <- vector("list", nfold)
    names(models) <- paste0("fold",1:nfold)
    fits <- vector("list", nfold)
    }

  if(sum(family %in% c("survival"))==0){
    for(ii in 1:nfold){
      #make train and test data for ii-th fold
      y.train <- ymat[idx.cv!=ii, ]
      y.test <-  ymat[idx.cv==ii, ]
      x.train <- xmat[idx.cv!=ii, ]
      x.test <-  xmat[idx.cv==ii, ]
      colnames(x.test) <- paste0("X",1:ncol(x.test))

    ######################should use train data

    fits[[ii]] <- multiFit(xmat=x.train, ymat=y.train,method,family)

    y.fitted[idx.cv==ii,] <- predict.multiFit(fits[[ii]], x.test)
    models[[ii]] <- fits[[ii]]$model
  }
  }else{
    for(ii in 1:nfold){
      # make train and test data for ii-th fold
      y.train.time <- ymat.time[idx.cv!=ii, ]
      y.train.status <- ymat.status[idx.cv!=ii, ]
      y.train <- vector("list",ny)
      for(j in 1:ny){
        y.train[[j]] <- Surv(y.train.time[,j], y.train.status[,j])
      }

      y.test.time <- ymat.time[idx.cv==ii, ]
      y.test.status <- ymat.status[idx.cv==ii, ]
      y.test <- vector("list",ny)
      for(j in 1:ny){
        y.test[[j]] <- Surv(y.test.time[,j], y.test.status[,j])
      }

      x.train <- xmat[idx.cv!=ii, ]
      x.test <-  xmat[idx.cv==ii, ]
      colnames(x.test) <- paste0("X",1:ncol(x.test))

      cvtrain<-cbind(y.train.time,y.train.status,x.train)
      names(cvtrain)[1:(2*ny)]<-c(paste0("ageph0", 1:ny),paste0("statusph0",1:ny))

      if(is.null(weights)) {cvweights<-NULL
      }else cvweights<-find_km_weights_mat(cvtrain)
      if(is.null(dist1)) {fits[[ii]] <- multiFit(xmat=x.train, ymat=y.train,method,family,dist1 = dist1,weights = cvweights)
      }else fits[[ii]] <- multiFit(xmat=x.train, ymat=y.train,method,family, dist1=dist1,weights= cvweights)



      y.fitted[idx.cv==ii,] <- predict.multiFit(fits[[ii]], x.test)
      models[[ii]] <- fits[[ii]]$model
    }

  }

  cv.multiFit.fits <- list(fit=fits,
                        y.fitted=y.fitted,
                        model=models,
                        method=method,
                        family=family)
  class(cv.multiFit.fits) <- "cv.multiFit"
  return(cv.multiFit.fits)

}

#'@method predict cv.multiFit
#'@export predict.cv.multiFit
#'@export
predict.cv.multiFit <- function(object, newdata, ...)
{
  ny <- ncol(object$y.fitted)
  nfold <- length(object$model)
  temp0 <- array(NA, c(nrow(newdata), ny, nfold))
  clas <- object$method

  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  for(ii in 1:ny) for(jj in 1:nfold)
  {
    model <- object$model[[jj]][[ii]]
    temp0[,ii,jj] <- object$fit[[jj]][["fit"]][[ii]]$predFun(model, newdata)
  }
  pred <- apply(temp0, c(1,2), mean)

  bindex<-object$family=="binomial"
  # predicted are probability should be within [0,1] for binary outcome
  pred[,bindex][pred[,bindex]>1] <- 1
  pred[,bindex][pred[,bindex]<0] <- 0

  return(pred)
}






