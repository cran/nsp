#' Change-point importance (prominence) plot
#'
#' This function produces a change-point prominence plot based on the NSP object provided. The heights of the bars are arranged in non-decreasing
#' order and correspond directly to the lengths of the NSP intervals of significance. Each bar is labelled as s-e where s (e) is the start (end) of the
#' corresponding NSP interval of significance, respectively. The change-points corresponding to the narrower intervals can be seen as more prominent.
#'
#' The NSP algorithm is described in P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' 
#' @param nsp.obj Object returned by one of the \code{nsp*} functions.
#' @return The function does not return a value.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @seealso \code{\link{draw_rects}}, \code{\link{draw_rects_advanced}}
#' @examples
#' set.seed(1)
#' f <- c(rep(0, 100), 1:100, rep(101, 100))
#' x.f <- f + 15 * stats::rnorm(300)
#' x.f.n <- nsp_poly(x.f, 100, "sim", deg=1)
#' cpt_importance(x.f.n)
#' @importFrom graphics barplot
#' @importFrom stats rnorm
#' @export




cpt_importance <- function(nsp.obj) {

	# Change-point prominence plot as described in Section 4 of the paper.
	# nsp.obj - quantity returned by one of the nsp_* functions.
	
	d <- dim(nsp.obj$intervals)
	if (d[1]) {
		heights <- nsp.obj$intervals[,2] - nsp.obj$intervals[,1]
		h.ord <- order(heights)
		labels <- paste(as.character(round(nsp.obj$intervals[h.ord,1])), "-", as.character(round(nsp.obj$intervals[h.ord,2])), sep = "")
		graphics::barplot(heights[h.ord], names.arg=labels)
	}
	else warning("No change-points to arrange in order of importance.")
	
}
