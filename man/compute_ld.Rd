% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_proxies.R
\name{compute_ld}
\alias{compute_ld}
\title{Compute two commonly used linkage disequilibrium statistics.}
\usage{
compute_ld(x, y)
}
\arguments{
\item{x}{a numeric vector of ones and zeros}

\item{y}{a numeric vector of ones and zeros}
}
\value{
A list with two items:
\describe{
 \item{R.squared}{Squared Pearson correlation coefficient.}
 \item{D.prime}{Coefficient of linkage disequilibrium D divided by the
  theoretical maximum.}
}
}
\description{
Compute R.squared and D.prime for two binary numeric vectors.
}
\details{
Find more details here:
\url{https://en.wikipedia.org/wiki/Linkage_disequilibrium}
}
\examples{
compute_ld(c(0,0,0,1,1,1), c(1,1,1,1,0,0))

}

