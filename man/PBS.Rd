\name{PBS}
\alias{PBS}
\title{Parametric bootstrap procedure for constructing prediction intervals of network meta-analysis}
\usage{
PBS(y, S, B=2000)
}
\arguments{
\item{y}{Contrast-based summary data of the outcome measure}

\item{S}{Covariance estimates of \code{y}}
\item{B}{Number of bootstrap resampling (default: 2000).}
}
\value{
The parametric bootstrap prediction intervals for network meta-analysis.
\itemize{
\item \code{Estimates}: Restricted maximum likelihood (REML) estimates, their SE, and 95\% Wald-type confidence intervals.
\item \code{Between-studies_SD}: Between-studies SD estimate.
\item \code{95\%PI}: 95\% prediction intervals by the parametric bootstrap.
}
}
\description{
Parametric bootstrap procedure for constructing prediction intervals of network meta-analysis.
}
\references{
Noma, H., Hamura, Y., Sugasawa, S. and Furukawa, T. A. (2022+).
Improved methods to construct prediction intervals for network meta-analysis. Forthcoming.
}
\examples{
data(dstr)
attach(dstr)

# Transforming the arm-level data to the contrast-based summaryies
edat <- data.edit(study,trt,d,n)

y <- edat$y
S <- edat$S

PBS(y,S,B=10)   # Results of the NMA analysis (log OR scale); B is recommended to be >= 1000.
}
