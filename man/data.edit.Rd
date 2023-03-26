\name{data.edit}
\alias{data.edit}
\title{Transforming arm-level data to contrast-based summary statistics}
\usage{
data.edit(study,trt,d,n)
}
\arguments{
\item{study}{Study ID}
\item{trt}{Numbered treatment (=1,2,...)}
\item{d}{Number of events}
\item{n}{Sample size}
}
\value{
Contrast-based summary statistics are generated.
\itemize{
\item \code{y}: Contrast-based summary estimates.
\item \code{S}: Vectored within-study covariance matrix.
}
}
\description{
Transforming arm-level data to contrast-based summary statistics.
}
\examples{
data(dstr)
attach(dstr)

edat <- data.edit(study,trt,d,n)
}
