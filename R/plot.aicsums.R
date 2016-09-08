##' @title Plot sum-of-AIC results
##'
##' @description S3 \code{\link{plot}} method for sum-of-AIC reuslts from \code{\link[optimus]{find_optimal}}.
##'
##' @param x an object of class \code{aicsums}, as produced by \code{\link[optimus]{find_optimal}}.
##' @param col point colour
##' @param pch point type
##' @param ... additional arguments to pass to \code{\link[graphics]{plot}}.
##'
##' @return A plot is drawn on the current graphics device
##'
##' @author Mitchell Lyons
##'
##' @export
##'
##' @rdname plot.aicsums
##'
##' @importFrom graphics plot
##'
##' @examples
##'
##' ## see ?find_optimal()

plot.aicsums <- function(x, col = "black", pch = 16, ...) {
  plot(x$sum_aic ~ x$nclusters, xlab = "Number of clusters", ylab = "Sum-of-AIC",
       col = col, pch = pch)

  invisible()
}


##' @title Plot more sum-of-AIC results
##'
##' @description S3 \code{\link{points}} method for sum-of-AIC reuslts from \code{\link[optimus]{find_optimal}}. Implemented to compare multiple outputs from \code{\link[optimus]{find_optimal}}.
##'
##' @param x an object of class \code{aicsums}, as produced by \code{\link[optimus]{find_optimal}}.
##' @param col point colour - random if not specified
##' @param pch point type - random if not specified
##' @param ... additional arguments to pass to \code{\link[graphics]{points}}.
##'
##' @return Points drawn on the current plot
##'
##' @author Mitchell Lyons
##'
##' @export
##'
##' @rdname points.aicsums
##'
##' @importFrom graphics points
##'
##' @examples
##'
##' ## see ?find_optimal()

points.aicsums <- function(x, col = sample(1:20,1), pch = sample(c(1:15,17:20),1), ...) {
  points(x$sum_aic ~ x$nclusters, xlab = "Number of clusters", ylab = "Sum-of-AIC",
         col = col, pch = pch)

  invisible()
}
