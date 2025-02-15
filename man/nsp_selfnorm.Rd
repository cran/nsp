% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsp_selfnorm.R
\name{nsp_selfnorm}
\alias{nsp_selfnorm}
\title{Self-normalised Narrowest Significance Pursuit algorithm with general covariates and user-specified threshold}
\usage{
nsp_selfnorm(
  y,
  x,
  M,
  lambda,
  power = 1/2,
  min.size = 20,
  eps = 0.03,
  c = exp(1 + 2 * eps),
  overlap = FALSE
)
}
\arguments{
\item{y}{A vector containing the data sequence being the response in the linear model y_t = beta(t) x_t + z_t.}

\item{x}{The design matrix in the regression model above, with the regressors as columns.}

\item{M}{The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
are used.}

\item{lambda}{The threshold parameter for measuring the significance of non-constancy (of the linear regression parameters), for use with the
self-normalised multiscale supremum-type deviation measure described in the paper.}

\item{power}{A parameter for the (rough) estimator of the global sum of squares of z_t; the span of the moving window in that estimator is
\code{min(n, max(round(n^power), min.size))}, where \code{n} is the length of \code{y}.}

\item{min.size}{(See immediately above.)}

\item{eps}{Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.}

\item{c}{Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.}

\item{overlap}{If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that
interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.}
}
\value{
A list with the following components:
\item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end,
respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are their midpoints.}
\item{threshold.used}{The threshold \code{lambda}.}
}
\description{
This function runs the self-normalised Narrowest Significance Pursuit (NSP) algorithm on data sequence \code{y} and design matrix \code{x}
to obtain localised regions (intervals) of the domain in which the parameters of the linear regression model y_t = beta(t) x_t + z_t significantly
depart from constancy (e.g. by containing change-points). For any interval considered by the algorithm,
significant departure from parameter constancy is achieved if the self-normalised multiscale
deviation measure (see Details for the literature reference) exceeds \code{lambda}. This function is
used by the higher-level function \code{\link{nsp_poly_selfnorm}}
(which estimates a suitable \code{lambda} so that a given global significance level is guaranteed), and human users may prefer to use that function
if \code{x} describe polynomial covariates; however, \code{nsp_selfnorm} can also be run directly, if desired.
The function assumes independence, symmetry and finite variance of the errors z_t, but little else; in particular they do not need to have a constant
variance across t.
}
\details{
The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
models", preprint.
}
\examples{
set.seed(1)
g <- c(rep(0, 100), rep(10, 100), rep(0, 100))
x.g <- g + stats::rnorm(300) * seq(from = 1, to = 4, length = 300)
wn003 <- sim_max_holder(100, 500, .03)
lambda <- as.numeric(stats::quantile(wn003, .9))
nsp_selfnorm(x.g, matrix(1, 300, 1), 100, lambda)
}
\seealso{
\code{\link{nsp_poly}}, \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_poly_selfnorm}}
}
\author{
Piotr Fryzlewicz, \email{p.fryzlewicz@lse.ac.uk}
}
