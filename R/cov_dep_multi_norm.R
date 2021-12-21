#' Simulate covariate-dependent multiscale sup-norm for use in NSP
#'
#' This function simulates the multiscale sup-norm adjusted for the form of the covariates, as described in Section 5.3
#' of the paper. This is done for i.i.d. N(0,1) innovations.
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param x The design matrix with the regressors (covariates) as columns.
#' @param N Desired number of simulated values of the norm.
#' @return Sample of size \code{N} containing the simulated norms.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{cov_dep_multi_norm_poly}}, \code{\link{sim_max_holder}}
#' @examples
#' set.seed(1)
#' g <- c(rep(0, 100), rep(2, 100))
#' x.g <- g + stats::rnorm(200)
#' mscale.norm.200 <- cov_dep_multi_norm(matrix(1, 200, 1), 100)
#' nsp_poly(x.g, 100, thresh.val = stats::quantile(mscale.norm.200, .95))
#' @importFrom stats rnorm quantile
#' @export




cov_dep_multi_norm <- function(x, N = 1000) {
	
	# Covariate-dependent multiresolution norm simulation
	
	d <- dim(x)

	res <- rep(0, N)
	
	for (i in 1:N) {
	
		z <- stats::rnorm(d[1])
	
		x.c <- cbind(z, x)
	
		x.c.ads <- all_dyadic_scans_array(x.c)
		
		res[i] <- check_interval_array(c(1, d[1]), x.c.ads, 0)
		
	}

	res
		
}
