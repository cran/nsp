#' nsp: Narrowest Significance Pursuit: Inference for Multiple Change-points in Linear Models
#'
#' Implementation of Narrowest Significance Pursuit (NSP), a general and
#' flexible methodology for automatically detecting localised regions in data sequences
#' which each must contain a change-point (understood as an abrupt change in the
#' parameters of an underlying linear model), at a prescribed global significance level.
#' NSP works with a wide range of distributional assumptions on the errors, and yields
#' exact desired finite-sample coverage probabilities, regardless of the form or number
#' of the regressors. A good place to start exploring the package are the \code{nsp*} functions.
#' @author Piotr Fryzlewicz, \email{p.fryzlewicz@@lse.ac.uk}
#' @references P. Fryzlewicz (2021) "Narrowest Significance Pursuit: inference for multiple change-points in linear
#' models", preprint.
#' @seealso \code{\link{nsp}}, \code{\link{nsp_poly}}, \code{\link{nsp_poly_ar}}, \code{\link{nsp_tvreg}}, \code{\link{nsp_selfnorm}},
#' \code{\link{nsp_poly_selfnorm}}
#' @aliases nsp-package
#' @keywords internal
"_PACKAGE"
