#' Self-normalised Narrowest Significance Pursuit algorithm with general covariates and user-specified threshold
#'
#' This function runs the self-normalised Narrowest Significance Pursuit (NSP) algorithm on data sequence \code{y} and design matrix \code{x}
#' to obtain localised regions (intervals) of the domain in which the parameters of the linear regression model y_t = beta(t) x_t + z_t significantly
#' depart from constancy (e.g. by containing change-points). For any interval considered by the algorithm, 
#' significant departure from parameter constancy is achieved if the self-normalised multiscale
#' deviation measure (see Details for the literature reference) exceeds \code{lambda}. This function is
#' used by the higher-level function \code{\link{nsp_poly_selfnorm}}
#' (which estimates a suitable \code{lambda} so that a given global significance level is guaranteed), and human users may prefer to use that function
#' if \code{x} describe polynomial covariates; however, \code{nsp_selfnorm} can also be run directly, if desired.
#' The function assumes independence, symmetry and finite variance of the errors z_t, but little else; in particular they do not need to have a constant
#' variance across t.
#' 
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param y A vector containing the data sequence being the response in the linear model y_t = beta(t) x_t + z_t.
#' @param x The design matrix in the regression model above, with the regressors as columns.
#' @param M The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
#' are used.
#' @param lambda The threshold parameter for measuring the significance of non-constancy (of the linear regression parameters), for use with the 
#' self-normalised multiscale supremum-type deviation measure described in the paper.
#' @param power A parameter for the (rough) estimator of the global sum of squares of z_t; the span of the moving window in that estimator is 
#' \code{min(n, max(round(n^power), min.size))}, where \code{n} is the length of \code{y}.
#' @param min.size (See immediately above.)
#' @param eps Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.
#' @param c Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.
#' @param overlap If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that 
#' interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.
#' @return A list with the following components:
#' \item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end, 
#'	respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are their midpoints.}
#' \item{threshold.used}{The threshold \code{lambda}.}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp_poly}}, \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_poly_selfnorm}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(10, 100), rep(0, 100))
#' x.g <- g + stats::rnorm(300) * seq(from = 1, to = 4, length = 300)
#' wn003 <- sim_max_holder(100, 500, .03)
#' lambda <- as.numeric(stats::quantile(wn003, .9))
#' nsp_selfnorm(x.g, matrix(1, 300, 1), 100, lambda)
#' @importFrom stats rnorm quantile
#' @export




nsp_selfnorm <- function(y, x, M, lambda, power = 1/2, min.size = 20, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
	
	d <- dim(x)
	
	x.c <- cbind(y, x)
	
	x.c.ads <- all_dyadic_scans_array(x.c)

	Vn2est <- est_var(y, x, power, min.size, TRUE)
	
	res <- ircs2sas(c(1, d[1]), x.c.ads, M, lambda, Vn2est, eps, c, overlap)
	
	intervals <- data.frame(t(order_chron(res)))
	colnames(intervals) <- c("starts", "ends", "values")
	intervals$midpoints <- floor((intervals$starts+intervals$ends)/2)
	
	list(intervals=intervals, threshold.used=lambda)
}
