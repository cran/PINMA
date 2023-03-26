\name{KR}
\alias{KR}
\title{Kenward-Roger-type adjustment for constructing prediction intervals of network meta-analysis}
\usage{
KR(y, S)
}
\arguments{
\item{y}{Contrast-based summary data of the outcome measure}

\item{S}{Covariance estimates of \code{y}}
}
\value{
Results of the Kenward-Roger-type adjustment for inference of multivariate random-effects model and prediction intervals for network meta-analysis.
\itemize{
\item \code{Estimates}: Restricted maximum likelihood (REML) estimates, their SE, and Wald-type 95\% confidence intervals by the Kenward-Roger-type adjustment.
\item \code{Between-studies_SD}: Between-studies SD estimate.
\item \code{95\%PI}: 95\% prediction intervals by the Kenward-Roger-type adjustment.
}
}
\description{
Kenward-Roger-type adjustment for constructing prediction intervals of network meta-analysis.
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

KR(y,S)    # Results of the NMA analysis (log OR scale)
}
