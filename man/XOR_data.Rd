\name{XOR_data}
\alias{XOR_data}
\docType{data}
\title{
Synthetic XOR Pattern Dataset
}
\description{
Simulated classification dataset containing 400 observations with 5 features demonstrating XOR patterns, linear class differences, and random noise.
}
\usage{data("XOR_data")}
\format{
A data frame with 400 rows and 6 variables:
\describe{
  \item{class}{Binary class labels (1 or 2)}
  \item{Variable_A}{Normally distributed with subtle class difference (delta mu=0.25)}
  \item{Variable_B}{High-variance normal distribution (sigma=3) with moderate class separation (delta mu=-0.7)}
  \item{Variable_C}{XOR pattern component 1 (mu=3 vs 10 between classes)}
  \item{Variable_D}{XOR pattern component 2 (mu=3 vs 10 between classes)}
  \item{Variable_E}{Uniform noise (1-10)}
}
}
\source{
Synthetic data generated with rnorm() and runif()
}
\examples{
data(XOR_data)
str(XOR_data)
summary(XOR_data)
}
\keyword{datasets}
