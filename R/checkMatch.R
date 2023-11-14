
check.match <- function(family, FUN) {

  continuous <- c(lm1,glm1,glmnet1,rpart1,glmnet.lasso,glmnet.ridge)
  binary <- c(glm1,glmnet1,rpart1,glmnet.lasso,glmnet.ridge,lda1,qda1,KNN1,svm1)
  sur <-c(ela1,surv,lasso1)
  clength <- length(continuous)
  blength <- length(binary)
  surlength <- length(sur)
  userDefine <- TRUE

  for (ii in 1:clength) {
    if (identical(FUN, continuous[[ii]]) ) {
      userDefine <- FALSE
      if (family == "gaussian") {
        return(TRUE)
      }
    }
  }
  for (ii in 1:blength) {
    if (identical(FUN, binary[[ii]]) ) {
      userDefine <- FALSE
      if (family == "binomial") {
        return(TRUE)
      }
    }
  }

  for (ii in 1:clength) {
    if (identical(FUN, sur[[ii]]) ) {
      userDefine <- FALSE
      if (family == "survival") {
        return(TRUE)
      }
    }
  }
# 	if (userDefine){
# 		warning("Please make sure that method is consistnet with outcome type")
# 		return(TRUE)
# 	}

  return(FALSE)

}

