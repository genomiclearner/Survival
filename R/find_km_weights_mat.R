#'@title Calculate Kaplan-Meier Weights
#'@aliases find_km_weights_mat
#'@description Produce the Kaplan-Meier weights
#'@details The Kaplan-Meier weights were defined as
#'\deqn{ w_{n1}=\frac{d_{(1)} }{n} , w_{ni}=\frac{d_{(i)} }{n-i+1} \prod^{i-1}_{j=1} (\frac{n-j}{n-j+1} )^{d_{(j)} }, i=2,...,n.}
#'@param multi_dat A a data frame that contains the time and status of
#'
#' multi-variate survival outcomes for every subject.
#'
#' The name of each column corresponding to the time or status must contain the
#'
#' characters “time” or “status”, respectively.
#'@param num_outcome Number of the outcomes
#'@export find_km_weights_mat
#'@return A matrix of Kaplan-Meier weights for each survival outcome.
#'@importFrom dplyr `%>%` mutate arrange
#'@references Huang J, Ma S, Xie H. Regularized estimation in the accelerated failure time model with high-dimensional covariates. Biometrics. 2006;62(3):813-820.
#'@examples
#' #multi-outcome
#' data(simdat_meps)
#' find_km_weights_mat(ymat,2)


find_km_weights_mat = function(multi_dat,num_outcome) {
  num_disease <- num_outcome
  dat <- multi_dat

  n = nrow(dat)

  weights_mat = matrix(nrow = n, ncol = num_disease)
  for (ii in 1:num_disease) {
    # rearrange dat -> ordered time of the current disease
    columns = which(is.element(colnames(dat),
                               c(paste0("time0", ii), paste0("status0", ii))))
    temp = dat[, columns] %>% mutate(order = c(1:n))
    names(temp)[1] <- ifelse(grepl("time", names(temp)[1]),"time","status")
    names(temp)[2] <-ifelse(grepl("status", names(temp)[2]),"status","time")

    # find weights
    temp = temp %>% arrange(time)
    # find K-M weights
    weights_vec = c()
    weights_vec[1] = temp$status[1] / n

    for (cc in 2:n) {
      coef_cc = temp$status[cc] / (n - cc + 1)
      product_cc = 1
      for (jj in 1:(cc-1)) {
        product_cc = product_cc * (((n - jj) / (n - jj + 1))^temp$status[jj])
      }
      weights_vec[cc] = coef_cc * product_cc
    }

    temp$weight = weights_vec

    temp = temp %>% arrange(order)
    weights_mat[, ii] = temp$weight
  }

  return(weights_mat)
}
