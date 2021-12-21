#' Compute the theoretical threshold for the multiscale sup-norm if the underlying distribution is standard normal
#'
#' This function computes the theoretical threshold, corresponding to the given significance level \code{alpha}, for the multiscale sup-norm
#' if the underlying distribution is standard normal.
#'
#' For the underlying theory, see Z. Kabluchko (2007) Extreme-value analysis of standardized Gaussian increments. Unpublished.
#' 
#' @param n The sample size.
#' @param alpha The significance level.
#' @param method "asymp" for the asymptotic method; "bound" for the Bonferroni method.
#' @return The desired threshold.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{cov_dep_multi_norm}}, \code{\link{cov_dep_multi_norm_poly}}, \code{\link{sim_max_holder}}
#' @examples
#' set.seed(1)
#' f <- c(1:100, 100:1, 1:100)
#' y <- f + stats::rnorm(300) * 15
#' x <- matrix(0, 300, 2)
#' x[,1] <- 1
#' x[,2] <- seq(from = 0, to = 1, length = 300)
#' nsp(y, x, 100, 15 * thresh_kab(300, .1))
#' @importFrom stats rnorm
#' @export



thresh_kab <- function(n, alpha = 0.1, method = "asymp") {
	
	an <- sqrt(2 * log(n)) + (1/2*log(log(n)) + log(0.8197466 / 2 /sqrt(pi))) / sqrt(2 * log(n))
	
	bn <- 1 / sqrt(2 * log(n))

	if (method == "bound") beta <- alpha/2
	else if (method == "asymp") beta <- 1 - sqrt(1-alpha)

	an + bn * log(1/log(1/(1 - beta)))
	
}
