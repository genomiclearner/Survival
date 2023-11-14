#'@title Modify Default Parameters For Base Learner
#'@description Modify default parameters for methods provided in the package
#'@param FUN Method
#'@param ... Modified arguments
#'@export modify.parameter
#'@return It returns a new function with modified parameters.
#'@examples
#' glmnet.lasso <- modify.parameter(glmnet1, alpha=1)
#' glmnet.ridge <- modify.parameter(glmnet1, alpha=0)




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
