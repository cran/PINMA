\name{tPI}
\alias{tPI}
\title{The ordinary t-approximation for constructing prediction intervals of network meta-analysis}
\usage{
tPI(y, S)
}
\arguments{
\item{y}{Contrast-based summary data of the outcome measure}

\item{S}{Covariance estimates of \code{y}}
}
\value{
The ordinary t-approximation prediction intervals for network meta-analysis.
\itemize{
\item \code{Estimates}: Restricted maximum likelihood (REML) estimates, their SE, and Wald-type 95\% confidence intervals.
\item \code{Between-studies_SD}: Between-studies SD estimate.
\item \code{95\%PI}: 95\% prediction intervals by the ordinary t-approximation.
}
}
\description{
The ordinary t-approximation for constructing prediction intervals of network meta-analysis.
}
\references{
Cooper, H., Hedges, L. V., and Valentine, J. C. (2009). \emph{The Handbook of Research Synthesis and Meta-Analysis}, 2nd edition. New York: Russell Sage Foundation.

Chaimani, A., and Salanti, G. (2015). Visualizing assumptions and results in network meta-analysis: the network graphs package. \emph{Stata Journal} \strong{15}, 905-920.
}
\examples{
data(dstr)
attach(dstr)

# Transforming the arm-level data to the contrast-based summaryies
edat <- data.edit(study,trt,d,n)

y <- edat$y
S <- edat$S

tPI(y,S)   # Results of the NMA analysis (log OR scale)
}
