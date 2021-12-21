#' Draw NSP intervals of significance as shaded rectangular areas on the current plot
#'
#' This function draws intervals of significance returned by one of the \code{nsp*} functions on the current plot. It shows them as shaded
#' rectangular areas (hence the name of the function).
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param nsp.obj Object returned by one of the \code{nsp*} functions.
#' @param yrange Vector of length two specifying the (lower, upper) vertical limit of the rectangles.
#' @param density Density of the shading.
#' @param col Colour of the shading.
#' @param x.axis.start Time index the x axis starts from. The NSP intervals of significance get shifted by \code{x.axis.start-1} prior to plotting.
#' @return The function does not return a value.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{draw_rects_advanced}}, \code{\link{nsp}}
#' @examples
#' set.seed(1)
#' h <- c(rep(0, 150), 1:150)
#' x.h <- h + stats::rnorm(300) * 50
#' x.h.n <- nsp_poly(x.h, 1000, "sim", deg=1)
#' draw_rects(x.h.n, c(-100, 100))
#' @importFrom graphics rect
#' @importFrom stats rnorm
#' @export



draw_rects <- function(nsp.obj, yrange, density = 10, col = "red", x.axis.start = 1) {

	# Draw intervals of significance, as shaded rectangular areas, on the current plot.
	# nsp.obj - quantity returned by one of the nsp_* functions.
	# yrange - vector of length two specifying the (lower, upper) vertical limit of the rectangles.
	# density - density of the shading; try using 10 or 20.
	# col - colour of the shading.
	# x.axis.start - time index the x axis stars from.
	
	d <- dim(nsp.obj$intervals)
	if (d[1]) for (i in 1:d[1]) {
		
		graphics::rect(nsp.obj$intervals[i,1]+x.axis.start-1, yrange[1], nsp.obj$intervals[i,2]+x.axis.start-1, yrange[2], density=density, col=col)		
		
	}
	
}
