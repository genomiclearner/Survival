% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mse.R
\name{mse}
\alias{mse}
\title{Mean Square Error Loss}
\usage{
mse(pre, object)
}
\arguments{
\item{pre}{A numeric vector of predicted survival time}

\item{object}{A dataframe or matrix, the first column is observed survival time the second column is survival status}
}
\value{
The value of MSE.
}
\description{
Compute the mean square error loss in survival analysis.
}
\details{
The mean square error is caculated on log-scale, only consider uncensored observations
}
\examples{
set.seed(1)
# simulate predicted survival time
x1 <- rnorm(100)+5
# simulate observed survival time
x2 <- rnorm(100)+10
# simulate observed survival status
x3 <- rep(c(1, 0), 50 )
# calculate C-index
mse(x1,  cbind(x2, x3))
}
