#' Narrowest Significance Pursuit algorithm for piecewise-polynomial signals with autoregression
#'
#' This function runs the Narrowest Significance Pursuit (NSP) algorithm on a data sequence \code{y} believed to follow the model
#' Phi(B)y_t = f_t + z_t, where f_t is a piecewise polynomial of degree \code{deg}, Phi(B) is a characteristic polynomial of autoregression of order
#' \code{ord} with unknown coefficients, and z_t is noise. The function returns localised regions (intervals) of the domain, such that each interval
#' must contain a change-point in the parameters of the polynomial f_t, or in the autoregressive parameters,
#' at the global significance level \code{alpha}.
#' For any interval considered by the algorithm, 
#' significant departure from parameter constancy is achieved if the multiscale
#' deviation measure (see Details for the literature reference) exceeds a threshold, which is either provided as input
#' or determined from the data (as a function of \code{alpha}). The function works best when the errors z_t are independent and
#' identically distributed Gaussians.
#' 
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint. For how to determine the \code{"univ"} threshold, see Kabluchko, Z. (2007) "Extreme-value analysis of standardized Gaussian increments".
#' Unpublished.
#' 
#' @param y A vector containing the data sequence.
#' @param ord The assumed order of the autoregression.
#' @param M The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
#' are used.
#' @param thresh.type \code{"univ"} if the significance threshold is to be determined as in Kabluchko (2007); \code{"sim"} for the degree-dependent 
#' threshold determined
#' by simulation (this is only available if the length of \code{y} does not exceed 2150; for longer sequences obtain a suitable threshold by running
#' \code{cov_dep_multi_norm_poly} first).
#' @param thresh.val Numerical value of the significance threshold (lambda in the paper); or \code{NULL} if the threshold is to be determined from
#' the data (see \code{thresh.type}).
#' @param sigma The standard deviation of the errors z_t; if \code{NULL} then will be estimated from the data via the MOLS estimator described in the paper.
#' @param alpha Desired maximum probability of obtaining an interval that does not contain a change-point (the significance threshold will be 
#' determined as a function of this parameter).
#' @param deg The degree of the polynomial pieces in f_t (0 for the piecewise-constant model; 1 for piecewise-linearity, etc.).
#' @param power A parameter for the MOLS estimator of sigma; the span of the moving window in the MOLS estimator is \code{min(n, max(round(n^power), min.size))},
#' where \code{n} is the length of \code{y} (minus \code{ord}).
#' @param min.size (See immediately above.)
#' @param overlap If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that 
#' interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.
#' @param buffer A non-negative integer specifying how many observations to leave out immediately to the left and to the right of a detected interval of
#' significance before recursively continuing the search
#' for the next interval.
#' @return A list with the following components:
#' \item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end, 
#'	respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are their midpoints.}
#' \item{threshold.used}{The threshold value.}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp}}, \code{\link{nsp_poly}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_selfnorm}}, \code{\link{nsp_poly_selfnorm}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(10, 100), rep(0, 100))
#' nsp_poly_ar(stats::filter(g + 2 * stats::rnorm(300), .5, "recursive"), thresh.type="sim")
#' @importFrom stats quantile filter
#' @export



nsp_poly_ar <- function(y, ord = 1, M = 1000, thresh.type = "univ", thresh.val = NULL, sigma = NULL, alpha = 0.1, deg = 0, power = 1/2, min.size = 20, overlap = FALSE, buffer = ord) {
		
	n <- length(y)
	
	x.c <- matrix(y, n, deg+1+ord)
	
	for (i in 1:(deg+1)) {
		
		x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
		
		
	}

	if (ord) for (i in 1:ord) x.c[(1+i):n,deg+1+i] <- y[1:(n-i)]

	x.c <- x.c[(ord+1):n,]
	
	y <- y[(ord+1):n]
	
	if (is.null(thresh.val)) {
		
		if (is.null(sigma)) sigma <- est_var(y, x.c, power, min.size)
		if (is.null(alpha)) alpha <- 0.1

		if (thresh.type == "univ") base.thresh <- thresh_kab(n-ord, alpha)
		if (thresh.type == "sim") base.thresh <- as.numeric(stats::quantile(nearest_simulated_norms(n-ord, deg), 1-alpha))
		
		thresh.val <- sigma * base.thresh
		
	}

	res <- nsp(y, x.c, M, thresh.val, overlap, buffer)

	res$intervals[, c(1, 2, 4)] <- res$intervals[, c(1, 2, 4)] + ord
	
	res

}


