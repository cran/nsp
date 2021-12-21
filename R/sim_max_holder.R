#' Simulate Holder-like norm of the Wiener process for use in self-normalised NSP
#'
#' This function simulates a sample of size \code{N} of values of the Holder-like norm of the Wiener process discretised with step 1/\code{n}.
#' The sample can then be used to find a suitable threshold for use with the self-normalised NSP.
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param n Number of equispaced sampling points for the Wiener process on \code{[0,1]}.
#' @param N Desired number of simulated values of the norm.
#' @param eps Parameter of the self-normalisation statistic as described in the paper.
#' @param c Parameter of the self-normalisation statistic as described in the paper; use default if unsure how to set.
#' @return Sample of size \code{N} containing the simulated norms.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp_selfnorm}}, \code{\link{nsp_poly_selfnorm}}, \code{\link{cov_dep_multi_norm}}, \code{\link{cov_dep_multi_norm_poly}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(10, 100), rep(0, 100))
#' x.g <- g + stats::rnorm(300) * seq(from = 1, to = 4, length = 300)
#' wn003 <- sim_max_holder(100, 500, .03)
#' lambda <- as.numeric(stats::quantile(wn003, .9))
#' nsp_poly_selfnorm(x.g, M = 100, thresh.val = lambda)
#' @importFrom stats rnorm quantile
#' @export


sim_max_holder <- function(n, N, eps, c = exp(1 + 2 * eps)) {
	
	# Simulate a sample of size N of values of the Holder-like norm of the Wiener process discretised with step 1/n.
	# See the "Example" in the description of the function "nsp_selfnorm".
	
	max.holder.sample <- rep(0, N)
	
	for (i in 1:N) {
		
		e <- stats::rnorm(n)
		
		max.holder.sample[i] <- max_holder(e, eps, c)
						
	}
	
	max.holder.sample
			
}
