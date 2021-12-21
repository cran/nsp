#' Plot NSP intervals of significance at appropriate places along the graph of data
#'
#' This function plots the intervals of significance returned by one of the \code{nsp*} functions, at appropriate places along the graph of data.
#' It shows them as shaded rectangular areas (hence the name of the function) "attached" to the graph of the data. Note: the data sequence \code{y} 
#' needs to have been plotted beforehand.
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param y The data.
#' @param nsp.obj Object returned by one of the \code{nsp*} functions with \code{y} on input.
#' @param half.height Half-height of each rectangle; if \code{NULL} then set to twice the estimated standard deviation of the data.
#' @param show.middles Whether to display lines corresponding to the midpoints of the rectanlges (rough change-point location estimates).
#' @param col.middles Colour of the midpoint lines.
#' @param lwd Line width for the midpoint lines.
#' @param density Density of the shading.
#' @param col.rects Colour of the shading.
#' @param x.axis.start Time index the x axis starts from. The NSP intervals of significance get shifted by \code{x.axis.start-1} prior to plotting.
#' @return The function does not return a value.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{draw_rects}}, \code{\link{nsp}}
#' @examples
#' set.seed(1)
#' f <- c(rep(0, 100), 1:100, rep(101, 100))
#' x.f <- f + 15 * stats::rnorm(300)
#' x.f.n <- nsp_poly(x.f, 100, "sim", deg=1)
#' stats::ts.plot(x.f)
#' draw_rects_advanced(x.f, x.f.n, density = 3)
#' @importFrom graphics rect
#' @importFrom stats mad ts.plot rnorm
#' @export



draw_rects_advanced <- function(y, nsp.obj, half.height = NULL, show.middles = TRUE, col.middles = "blue", lwd = 3, density = 10, col.rects = "red", x.axis.start = 1) {

	
	loc.est <- round((nsp.obj$intervals[,1] + nsp.obj$intervals[,2])/2)
	
	if (is.null(half.height)) half.height <- 2 * stats::mad(diff(y)/sqrt(2))
	
	centres.y <- y[loc.est]
	
	d <- dim(nsp.obj$intervals)
	if (d[1]) for (i in 1:d[1]) {
				
		graphics::rect(nsp.obj$intervals[i,1]+x.axis.start-1, centres.y[i]-half.height, nsp.obj$intervals[i,2]+x.axis.start-1, centres.y[i]+half.height, density=density, col=col.rects)
		
		if (show.middles) graphics::lines(rep(loc.est[i]+x.axis.start-1, 2), c(centres.y[i]-half.height, centres.y[i]+half.height), col=col.middles, lwd = lwd)
		
	}
	
}
