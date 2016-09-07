#' @title Determine the characteristicvariables (e.g. species) of a clustering solution (e.g. classificaiton)
#'
#' @description \code{get_characteristic} takes a clustering solution, fits models based on the underlying multivariate data, and calculates delta AIC value for each variable. In Ecology, particularly vegetation science, this is the process of determining characteristic (or diagnostic/indicator) species of a classification.
#'
#' @param data a data frame (or object that can be coerced by \code{\link[base]{as.data.frame}} containing the "raw" multivariate data. This is not necessarily the data used by the clustering algorithm - it is the data on which you are testing the predictive ablity of the clustering solution.
#' @param clustering a clustering solution for \code{data}, that is, a vector of cluster labels (that can be coerced by \code{\link[base]{as.factor}}). The number of cluster labels must match the number of rows of the object supplied in the \code{data} argument. The solution could for exmaple come form a call to \code{\link[stats]{cutree}}, see Examples
#' @param family a character string denoting the error distribution to be used for model fitting. The options are similar to those in \code{\link[stats]{family}}, but are more limited - see Details.
#' @param K number of trials in binomial regression. By default, K=1 for presence-absence data (logistic regression).
#'
#' @details \code{get_characteristic} is built on the premise that a \emph{good} clustering solution (i.e. a classification) should provide information about the composition and abundance of the multivariate data it is classifying. A natural way to formalize this is with a predictive model, where group membership (clusters) is the predictor, and the multivariate data (site by variables matrix) is the response. \code{get_characteristic} fits linear models to each variable, and calculates the delta AIC (that is, to the corresponding null model). The larger the delta AIC, the \emph{more} characteristic that species is for the given classification. Lyons et al. (2016) provides background, a detailed description of the methodology, and application of delta AIC on both real and simulated ecological multivariate abundance data.
#'
#' At present, \code{get_characteristic} supports the following error distributions for model fitting:
#' \itemize{
#'   \item Gaussian (LM)
#'   \item Negative Binomial (GLM with log link)
#'   \item Poisson (GLM with log link)
#'   \item Binomial (GLM with logit link)
#'   \item Ordinal (Proportional odds model with logit link)
#' }
#'
#' Gaussian LMs should be used for 'normal' data. Negative Binomial and Poisson GLMs shold be used for count data. Binomial GLMs should be used for binary and presence/absence data (when \code{K=1}), or trials data (e.g. frequency scores). If Binomial regression is being used with \code{K>1}, then \code{data} should be numerical values between 0 and 1, interpreted as the proportion of successful cases, where the total number of cases is given by \code{K} (see Details in \code{\link[stats]{family}}). Ordinal regression should be used for ordinal data, for example, cover-abundance scores. LMs fit via \code{\link[mvabund]{manylm}}; GLMs fit via \code{\link[mvabund]{manyglm}}; proportional odds model fit via \code{\link[ordinal]{clm}}.
#'
#' @return a data frame containing the delta AIC values for each variable in \code{data}. The object is of class \code{daic}.
#'
#' Attributes for the data frame are:
#'
#' \describe{
##'   \item{\code{family}}{ which error distribution was used for modelling, see Arguments}
##'   \item{\code{K}}{ number of cases for Binomial regression, see Arguments}
##' }
#'
#' @author Mitchell Lyons
#'
#' @references Lyons et al. 2016. Model-based assessment of ecological community classifications. \emph{Journal of Vegetation Science}, \strong{27 (4)}: 704--715.
#'
#' @seealso \code{\link[optimus]{find_optimal}}, S3 print function for 'daic' class, S3 residual plotting function
#'
#' @keywords characteristic
#'
#' @examples
#'
#' ## Find characteristic species in a classificaiton of the swamps data
#' ## ==================================================================
#'
#' ## load data, remove transect column and cluster
#' data(swamps)
#' swamps <- swamps[,-1]
#' ## perhaps not the best clustering option, but this is base R
#' swamps_clust <- hclust(d = dist(x = swamps, method = "minkowski"),
#'                        method = "average")
#'
#' # calculate delta aics for
#' swamps_daic <- get_characteristic(data = swamps, clustering = cutree(swamps_clust, 10),
#' family = "poisson")
#' swamps_daic
#'
#' # print(swamps_clust_aics)
#'
#'
#' ## 2) find an optimal partioning with find_optimal(),
#' ## then call get_characteristic() on it
#'
#' @export


get_characteristic <- function(data, clustering, family, K = 1) {
  data <- as.data.frame(data)

  # test specified family is supported
  supported_fams <- c("gaussian", "negative.binomial", "poisson", "binomial", "ordinal")
  if (!family %in% supported_fams) {
    stop(paste0("family specified is not valid (typo?) or not yet supported, please choose from: ",
                paste(supported_fams, collapse = ", ")))
  }

  # test multivar data input (i.e. is char/num/factor)
  if (family != "ordinal" & any(unlist(lapply(data, class)) %in% c("factor", "character"))) {
    stop("some of the input data are factors or characters, are you looking for ordinal regression?")
  }

  # as.character, to be safe
  cluster_labels <- as.character(clustering)

  # get nclusters
  nclusters <- length(unique(cluster_labels))

  # fit models to calculate delta aic
  if (family == "gaussian") {
    delta_aics <- gaussian_char(clusterSolution = cluster_labels, data = data, nclusters = nclusters)
  }

  if (family == "negative.binomial") {
    delta_aics <- negbin_char(clusterSolution = cluster_labels, data = data)
  }

  if (family == "poisson") {
    delta_aics <- poisson_char(clusterSolution = cluster_labels, data = data)
  }

  if (family == "binomial") {
    delta_aics <- binomial_char(clusterSolution = cluster_labels, data = data, K=K)
  }

  if (family == "ordinal") {
    delta_aics <- ordinal_char(clusterSolution = cluster_labels, data = data)
  }

  # attributes/class for delta_aics
  attr(delta_aics, "family") <- family
  attr(delta_aics, "K") <- K
  class(delta_aics) <- c("daic","data.frame")

  # return data frame ready for printing
  delta_aics
}

