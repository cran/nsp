#' Simulate covariate-dependent multiscale sup-norm for use in NSP, for piecewise-polynomial models
#'
#' This function simulates the multiscale sup-norm adjusted for the form of the covariates, as described in Section 5.3
#' of the paper, for piecewise-polynomial models of degree \code{deg}. This is done for i.i.d. N(0,1) innovations.
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param n The data length (for which the multiscale norm is to be simulated)
#' @param deg The degree of the polynomial model (0 for the piecewise-constant model; 1 for piecewise-linearity, etc.).
#' @param N Desired number of simulated values of the norm.
#' @return Sample of size \code{N} containing the simulated norms.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{cov_dep_multi_norm}}, \code{\link{sim_max_holder}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(2, 100))
#' x.g <- g + stats::rnorm(200)
#' mscale.norm.200 <- cov_dep_multi_norm_poly(200, 0, 100)
#' nsp_poly(x.g, 100, thresh.val = stats::quantile(mscale.norm.200, .95))
#' @importFrom stats rnorm quantile
#' @export




cov_dep_multi_norm_poly <- function(n, deg, N = 10000) {
		
	x.c <- matrix(0, n, deg+1)
	
	for (i in 1:(deg+1)) {
		
		x.c[,i] <- seq(from=0, to=1, length=n)^(i-1)
				
	}

	cov_dep_multi_norm(x.c, N)	
	
}
