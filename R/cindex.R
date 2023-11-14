#'@title Concordance Index
#'@description Compute the C-index in survival analysis.
#'@details The most frequently used evaluation metric of survival models is the concordance index (C-index).
#'It is a measure of rank correlation between predicted risk scores. Only use uncensored observations
#'@param pre A numeric vector of predicted survival time.
#'@param object A dataframe or matrix, the first column is observed survival time the second column is survival status.
#'@export cindex
#'@return The value of C-index.
#'@examples
#' set.seed(1)
#' # simulate predicted survival time
#' x1 <- rnorm(100)
#' # simulate observed survival time
#' x2 <- rnorm(100)
#' # simulate observed survival status
#' x3 <- rep(c(1, 0), 50)
#' # calculate C-index
#' cindex(x1, cbind(x2, x3))

cindex <- function(pre, object){
  total.pairs <- 0
  c <- 0
  test <- cbind(pre, object)
  test <- test[order(test[,2]),]
  for (i in 1:nrow(test)){
    if(i ==nrow(test)) break
    else if(test[i,3]==0) next # only consider uncensored observations
    else{
      for (j in (i+1):nrow(test)){
        total.pairs <- total.pairs+1
        if(test[j,2] == test[i,2]) {if(test[j,1] > test[i,1] & test[j,1] == test[i,1]) c <- c + 0.5}
        #"if(test[j,1] > test[i,1] & test[j,1] == test[i,1]) c <- c + 0.5"
        #if we want add all possible results we can use straightly "c <- c + 0.5"

        if(test[j,2] > test[i,2]){if(test[j,1] > test[i,1]) {c <- c+1}}
        if(test[j,2] > test[i,2]){if(test[j,1] == test[i,1]) {c <- c+0.5}}
      }
    }

  }
  return (c/total.pairs)
}
