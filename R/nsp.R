#' Narrowest Significance Pursuit algorithm with general covariates and user-specified threshold
#'
#' This function runs the bare-bones Narrowest Significance Pursuit (NSP) algorithm on data sequence \code{y} and design matrix \code{x}
#' to obtain localised regions (intervals) of the domain in which the parameters of the linear regression model y_t = beta(t) x_t + z_t significantly
#' depart from constancy (e.g. by containing change-points). For any interval considered by the algorithm, 
#' significance is achieved if the multiscale supremum-type
#' deviation measure (see Details for the literature reference) exceeds \code{lambda}. This function is
#' mainly to be used by the higher-level functions \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}} and \code{\link{nsp_tvreg}}
#' (which estimate a suitable \code{lambda} so that a given global significance level is guaranteed), and human users may prefer to use those functions
#' instead; however, \code{nsp} can also be run directly, if desired.
#' The function works best when the errors z_t in the linear regression formulation y_t = beta(t) x_t + z_t are independent and
#' identically distributed Gaussians.
#' 
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param y A vector containing the data sequence being the response in the linear model y_t = beta(t) x_t + z_t.
#' @param x The design matrix in the regression model above, with the regressors as columns.
#' @param M The minimum number of intervals considered at each recursive stage, unless the number of all intervals is smaller, in which case all intervals
#' are used.
#' @param lambda The threshold parameter for measuring the significance of non-constancy (of the linear regression parameters), for use with the multiscale
#' supremum-type deviation measure described in the paper.
#' @param overlap If \code{FALSE}, then on discovering a significant interval, the search continues recursively to the left and to the right of that 
#' interval. If \code{TRUE}, then the search continues to the left and to the right of the midpoint of that interval.
#' @param buffer A non-negative integer specifying how many observations to leave out immediately to the left and to the right of a detected interval of
#' significance before recursively continuing the search
#' for the next interval.
#' @return A list with the following components:
#' \item{intervals}{A data frame containing the estimated intervals of significance: \code{starts} and \code{ends} is where the intervals start and end, 
#'	respectively; \code{values} are the values of the deviation measure on each given interval; \code{midpoints} are the midpoints of the intervals.}
#' \item{threshold.used}{The threshold \code{lambda}.}
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_selfnorm}}, \code{\link{nsp_poly_selfnorm}}
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



nsp <- function(y, x, M, lambda, overlap = FALSE, buffer = 0) {
		
	d <- dim(x)
	
	x.c <- cbind(y, x)
	
	x.c.ads <- all_dyadic_scans_array(x.c)
	
	res <- iter_random_checks_scan_array(c(1, d[1]), x.c.ads, M, lambda, overlap, buffer)
	
	intervals <- data.frame(t(order_chron(res)))
	colnames(intervals) <- c("starts", "ends", "values")
	intervals$midpoints <- floor((intervals$starts+intervals$ends)/2)
	
	list(intervals=intervals, threshold.used=lambda)
}
