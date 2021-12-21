#' Self-normalised Narrowest Significance Pursuit algorithm for piecewise-polynomial signals
#'
#' This function runs the Narrowest Significance Pursuit (NSP) algorithm on a data sequence \code{y} believed to follow the model
#' y_t = f_t + z_t, where f_t is a piecewise polynomial of degree \code{deg}, and z_t is noise. It returns localised regions (intervals) of the 
#' domain, such that each interval must contain a change-point in the parameters of the polynomial f_t
#' at the global significance level \code{alpha}.
#' For any interval considered by the algorithm, 
#' significant departure from parameter constancy is achieved if the multiscale
#' deviation measure (see Details for the literature reference) exceeds a threshold, which is either provided as input
#' or determined from the data (as a function of \code{alpha}). The function assumes independence, symmetry and finite variance of the 
#' errors z_t, but little else; in particular they do not need to have a constant variance across t.
#' 
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param y A vector containing the data sequence.
#' @param M The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
#' are used.
#' @param thresh.val Numerical value of the significance threshold (lambda in the paper); or \code{NULL} if the threshold is to be determined from
#' the data.
#' @param power A parameter for the (rough) estimator of the global sum of squares of z_t; the span of the moving window in that estimator is 
#' \code{min(n, max(round(n^power), min.size))}, where \code{n} is the length of \code{y}.
#' @param min.size (See immediately above.)
#' @param alpha Desired maximum probability of obtaining an interval that does not contain a change-point (the significance threshold will be 
#' determined as a function of this parameter).
#' @param deg The degree of the polynomial pieces in f_t (0 for the piecewise-constant model; 1 for piecewise-linearity, etc.).
#' @param eps Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.
#' @param c Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.
#' @param overlap If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that 
#' interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.
#' @return A list with the following components:
#' \item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end, 
#'	respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are their midpoints.}
#' \item{threshold.used}{The threshold value.}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp_poly}}, \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_selfnorm}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(10, 100), rep(0, 100))
#' x.g <- g + stats::rnorm(300) * seq(from = 1, to = 4, length = 300)
#' nsp_poly_selfnorm(x.g, 100)
#' @importFrom stats quantile rnorm
#' @export





nsp_poly_selfnorm <- function(y, M = 1000, thresh.val = NULL, power = 1/2, min.size = 20, alpha = 0.1, deg = 0, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
	
	# Self-normalised NSP when x is believed to be a piecewise polynomial plus (possibly heterogeneous and/or heavy-tailed) noise.
	# x - data (referred to in the paper as Y).
	# M - number of intervals to draw; this implementation uses a deterministic equispaced grid for drawing intervals.
	# thresh - threshold value (lambda_alpha in the paper), best left at NULL but specified implicitly via alpha.
	# power - parameter for estimating the global RSS, best left at 1/2.
	# min.size - parameter for estimating the global RSS, best left at 20.
	# alpha - desired maximum probability of obtaining an interval that does not cover a true change-point.
	# deg - degree of the underlying polynomial
	# eps - epsilon from the paper.
	# c - c (linked to epsilon) from the paper.
	# overlap - FALSE means no overlap, TRUE means an overlap as specified in Section 4 of the paper.
	#
	# Returns:
	# object of a class with the following two fields -
	#
	# intervals - data frame whose first two columns are start- and end-points (respectively) of the detected intervals of significance;
	#    the third column are the corresponding self-normalised scans; the fourth are the mid-points of the intervals.
	# threshold.used - thresh.	
	
	n <- length(y)
	
	x.c <- matrix(y, n, deg+1)
	
	for (i in 1:(deg+1)) {
		
		x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
		
		
	}
	
	if (is.null(thresh.val)) {
		
		if (is.null(alpha)) alpha <- 0.1
		varname <- paste("wiener.holder_", as.character(eps), sep="")
		if (exists(varname)) wh <- get(varname) else stop("This value of eps does not correspond to pre-computed threshold values; use sim_max_holder first.")
		thresh.val <- as.numeric(stats::quantile(wh, 1-alpha))
		
	}

	nsp_selfnorm(y, x.c, M, thresh.val, power, min.size, eps, c, overlap)
	
}
