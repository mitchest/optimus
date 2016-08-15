##' @title Plot sum-of-AIC results
##'
##' @description S3 \code{\link{plot}} method for sum-of-AIC reuslts from \code{\link[optimus]{find_optimal}}.
##'
##' @param x an object of class \code{aicsums}, as produced by \code{\link[optimus]{find_optimal}}.
##' @param ... additional arguments to pass to \code{\link[graphics]{plot}}.
##'
##' @export
##'
##' @rdname plot.aicsums
##'
##' @importFrom graphics plot
##'
##' @examples
##'
##' ## Example of a find_optimal() call and subsequent plot, and how to interpret

plot.aicsums <- function(x, ...) {
  plot(x$sum_aic ~ x$nclusters)
}
