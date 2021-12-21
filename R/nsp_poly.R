#' Narrowest Significance Pursuit algorithm for piecewise-polynomial signals
#'
#' This function runs the Narrowest Significance Pursuit (NSP) algorithm on a data sequence \code{y} believed to follow the model
#' y_t = f_t + z_t, where f_t is a piecewise polynomial of degree \code{deg}, and z_t is noise. It returns localised regions (intervals) of the 
#' domain, such that each interval must contain a change-point in the parameters of the polynomial f_t
#' at the global significance level \code{alpha}.
#' For any interval considered by the algorithm, 
#' significant departure from parameter constancy is achieved if the multiscale supremum-type
#' deviation measure (see Details for the literature reference) exceeds a threshold, which is either provided as input
#' or determined from the data (as a function of \code{alpha}). The function works best when the errors z_t are independent and
#' identically distributed Gaussians.
#' 
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint. For how to determine the \code{"univ"} threshold, see Kabluchko, Z. (2007) "Extreme-value analysis of standardized Gaussian increments".
#' Unpublished.
#' 
#' @param y A vector containing the data sequence.
#' @param M The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
#' are used.
#' @param thresh.type \code{"univ"} if the significance threshold is to be determined as in Kabluchko (2007); \code{"sim"} for the degree-dependent 
#' threshold determined
#' by simulation (this is only available if the length of \code{y} does not exceed 2150; for longer sequences obtain a suitable threshold by running
#' \code{cov_dep_multi_norm_poly} first).
#' @param thresh.val Numerical value of the significance threshold (lambda in the paper); or \code{NULL} if the threshold is to be determined from
#' the data (see \code{thresh.type}).
#' @param sigma The standard deviation of the errors z_t; if \code{NULL} then will be estimated from the data via Median Absolute Deviation (for i.i.d.
#' Gaussian sequences) of the first difference.
#' @param alpha Desired maximum probability of obtaining an interval that does not contain a change-point (the significance threshold will be 
#' determined as a function of this parameter).
#' @param deg The degree of the polynomial pieces in f_t (0 for the piecewise-constant model; 1 for piecewise-linearity, etc.).
#' @param overlap If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that 
#' interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.
#' @return A list with the following components:
#' \item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end, 
#'	respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are their midpoints.}
#' \item{threshold.used}{The threshold value.}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_selfnorm}}, \code{\link{nsp_poly_selfnorm}}
#' @examples
#' set.seed(1)
#' f <- c(1:100, 100:1, 1:100)
#' y <- f + stats::rnorm(300) * 15
#' nsp_poly(y, 100, deg = 1)
#' @importFrom stats mad quantile rnorm
#' @export


nsp_poly <- function(y, M = 1000, thresh.type = "univ", thresh.val = NULL, sigma = NULL, alpha = 0.1, deg = 0, overlap = FALSE) {
	
	n <- length(y)
	
	x.c <- matrix(y, n, deg+1)
	
	for (i in 1:(deg+1)) {
		
		x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
		
	}
	
	if (is.null(thresh.val)) {
		
		if (is.null(sigma)) sigma <- stats::mad(diff(y/sqrt(2)))
		if (is.null(alpha)) alpha <- 0.1

		if (thresh.type == "univ") base.thresh <- thresh_kab(n, alpha)
		if (thresh.type == "sim") base.thresh <- as.numeric(stats::quantile(nearest_simulated_norms(n, deg), 1-alpha))
		
		thresh.val <- sigma * base.thresh
		
	}

	nsp(y, x.c, M, thresh.val, overlap)
	
}
