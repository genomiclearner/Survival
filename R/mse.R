#'@title Mean Square Error Loss
#'@description Compute the mean square error loss in survival analysis.
#'@details The mean square error is caculated on log-scale, only consider uncensored observations
#'@param pre A numeric vector of predicted survival time
#'@param object A dataframe or matrix, the first column is observed survival time the second column is survival status
#'@export mse
#'@return The value of MSE.
#'@examples
#' set.seed(1)
#' # simulate predicted survival time
#' x1 <- rnorm(100)+5
#' # simulate observed survival time
#' x2 <- rnorm(100)+10
#' # simulate observed survival status
#' x3 <- rep(c(1, 0), 50 )
#' # calculate C-index
#' mse(x1,  cbind(x2, x3))

mse <- function(pre, object){
  # pre_time: the predicted event time
  # obj: 2 columns, time & status
  ind <- which(object[,2] == 1)
  event_time <- object[ind,1]
  pre_event_time <- pre[ind]
  names(event_time) <- names(pre_event_time)
  return (mean((log(event_time) - log(pre_event_time))^2))
}

