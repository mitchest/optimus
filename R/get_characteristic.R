#' @title Determine the characteristic variables (e.g. species) of a clustering solution (e.g. classificaiton)
#'
#' @description \code{get_characteristic} takes a clustering solution, fits models based on the underlying multivariate data, and determines 'important' variables for the clustering solution. In Ecology, particularly vegetation science, this is the process of determining characteristic (or diagnostic/indicator) species of a classification.
#'
#' @param data a data frame (or object that can be coerced by \code{\link[base]{as.data.frame}} containing the "raw" multivariate data. This is not necessarily the data used by the clustering algorithm - it is the data on which you are testing the predictive ablity of the clustering solution.
#' @param clustering a clustering solution for \code{data}, that is, a vector of cluster labels (that can be coerced by \code{\link[base]{as.factor}}). The number of cluster labels must match the number of rows of the object supplied in the \code{data} argument. The solution could for example come form a call to \code{\link[stats]{cutree}}, see Examples
#' @param family a character string denoting the error distribution to be used for model fitting. The options are similar to those in \code{\link[stats]{family}}, but are more limited - see Details.
#' @param type a character string, one of \code{"per.cluster"} or \code{"global"}, denoting the type of characteristic variables (species). See Details.
#' @param signif logical, denoting whether \emph{significance} should be returned also when \code{"type=per.cluster"}. Ignored if \code{"type=global"}. See Details. Minimal additional overhead is required if \code{TRUE}.
#' @param K number of trials in binomial regression. By default, K=1 for presence-absence data (with cloglog link).
#'
#' @details \code{get_characteristic} is built on the premise that a \emph{good} clustering solution (i.e. a classification) should provide information about the composition and abundance of the multivariate data it is classifying. A natural way to formalize this is with a predictive model, where group membership (clusters) is the predictor, and the multivariate data (site by variables matrix) is the response. \code{get_characteristic} fits linear models to each variable. If \code{type = "per.cluster"} the coefficients corresponding to each level of the clustering solution for each variable are used to define the characteristic variables for each cluster level. If \code{type = "global"}, characteristic variables are determined (via delta AIC - larger values = more important) for the overall classification. If \code{signif = TRUE}, delta AIC (that is, to the corresponding null model) and the coefficient standard errors are also retuend with the per-cluster characteristic variables. We loosely define that the larger the coefficient (with larger delta AIC values and smaller standard errors guiding \emph{significance}), the \emph{more} characteristic that varaible (species) is. Lyons et al. (2016) provides background, a detailed description of the methodology, and application of delta AIC on both real and simulated ecological multivariate abundance data.
#'
#' At present, \code{get_characteristic} supports the following error distributions for model fitting:
#' \itemize{
#'   \item Gaussian (LM)
#'   \item Negative Binomial (GLM with log link)
#'   \item Poisson (GLM with log link)
#'   \item Binomial (GLM with cloglog link for binary data, logit link otherwise)
#'   \item Ordinal (Proportional odds model with logit link)
#' }
#'
#' Gaussian LMs should be used for 'normal' data. Negative Binomial and Poisson GLMs shold be used for count data. Binomial GLMs should be used for binary and presence/absence data (when \code{K=1}), or trials data (e.g. frequency scores). If Binomial regression is being used with \code{K>1}, then \code{data} should be numerical values between 0 and 1, interpreted as the proportion of successful cases, where the total number of cases is given by \code{K} (see Details in \code{\link[stats]{family}}). Ordinal regression should be used for ordinal data, for example, cover-abundance scores. For ordinal regression, data should be supplied as either 1) factors, with the appropriate ordinal level order specified (see \code{\link[base]{levels}}) or 2) numeric, which will be coerced into a factor with levels ordered in numerical order (e.g. cover-abundance/numeric response scores). LMs fit via \code{\link[mvabund]{manylm}}; GLMs fit via \code{\link[mvabund]{manyglm}}; proportional odds model fit via \code{\link[ordinal]{clm}}.
#'
#' @return either a list of sorted characteristic variables for each cluster (of class \code{perclustchar}) or a data frame containing the delta AIC values for each variable (of class \code{globalchar}). If \code{signif=} is not \code{"none"}, then the corresponding significance metrics are appended.
#'
#' Attributes for the object are:
#'
#' \describe{
#'   \item{\code{family}}{ which error distribution was used for modelling, see Arguments}
#'   \item{\code{type}}{the type of characteristic variables calcualted, see Arguments}
#'   \item{\code{K}}{ number of cases for Binomial regression, see Arguments}
#' }
#'
#' @author Mitchell Lyons
#'
#' @references Lyons et al. 2016. Model-based assessment of ecological community classifications. \emph{Journal of Vegetation Science}, \strong{27 (4)}: 704--715.
#'
#' @seealso \code{\link[optimus]{find_optimal}}, S3 for print 'top-n' variables for each cluster, S3 for residual plots (at some stage)
#'
#' @keywords characteristic, diagnostic, indicator
#'
#' @examples
#'
#' ## Prep the 'swamps' data
#' ## ======================
#'
#' data(swamps) # see ?swamps
#' swamps <- swamps[,-1]
#'
#' ## Find characteristic species in a classificaiton of the swamps data
#' ## ==================================================================
#'
#' ## perhaps not the best clustering option, but this is base R
#' swamps_hclust <- hclust(d = dist(x = log1p(swamps), method = "canberra"),
#'                        method = "complete")
#'
#' # calculate per cluster characteristic species
#' swamps_char <- get_characteristic(data = swamps,
#' clustering = cutree(tree = swamps_hclust, k = 10), family = "poisson",
#' type = "per.cluster")
#'
#' # look at the top 10 characteristic species for cluster 1
#' head(swamps_char[[1]], 10)
#'
#' # calculate global characteristic species
#' swamps_char <- get_characteristic(data = swamps,
#' clustering = cutree(tree = swamps_hclust, k = 10), family = "poisson",
#' type = "global")
#'
#' # top 10 characteristic species for the whole classification
#' head(swamps_char, 10)
#'
#' ## See vignette for more explanation than this example
#' ## ============================================================
#'
#' @export


get_characteristic <- function(data, clustering, family, type="per.cluster", signif=TRUE, K = 1) {
  data <- as.data.frame(data)

  # test specified family is supported
  supported_fams <- c("gaussian", "negative.binomial", "poisson", "binomial", "ordinal")
  if (!family %in% supported_fams) {
    stop("family specified is not valid (typo?) or not yet supported, please choose from: ",
                paste(supported_fams, collapse = ", "))
  }

  # test multivar data input (i.e. is char/num/factor)
  if (family != "ordinal" & any(unlist(lapply(data, class)) %in% c("factor", "character"))) {
    stop("some of the input data are factors or characters, are you looking for ordinal regression?")
  }

  if (family == "ordinal") {
    if (all(unlist(lapply(data, is.factor)))) {
      message("All data are factors, ordinal regression will use factor levels as is - ensure they are correct")
    }
    if (!all(unlist(lapply(data, is.factor))) & any(unlist(lapply(data, is.factor)))) {
      stop("Some data are factors, some are not - don't know how to proceed with ordinal regression")
    }
    if (all(unlist(lapply(data, is.numeric)))) {
      message("All data are numeric - will coerce to factor and levels will be in numeric order")
    }
  }

  # as.character, to be safe
  cluster_labels <- as.character(clustering)

  # get nclusters
  nclusters <- length(unique(cluster_labels))

  # fit models to calculate coefs for each varaible/cluster
  if (type == "per.cluster") {
    if (family == "gaussian") {
      cluster_coefs <- gaussian_char(clusterSolution = cluster_labels, data = data, nclusters = nclusters, type = type, signif = signif)
    }

    if (family == "negative.binomial") {
      cluster_coefs <- negbin_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    if (family == "poisson") {
      cluster_coefs <- poisson_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    if (family == "binomial") {
      cluster_coefs <- binomial_char(clusterSolution = cluster_labels, data = data, K = K, type = type, signif = signif)
    }

    if (family == "ordinal") {
      cluster_coefs <- ordinal_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    # attributes/class for cluster_coefs
    attr(cluster_coefs, "family") <- family
    attr(cluster_coefs, "K") <- K
    attr(cluster_coefs, "type") <- type
    class(cluster_coefs) <- c("perclustchar","list")

    # return data frame ready for printing
    return(cluster_coefs)
  }

  # fit models to calculate delta aic
  if (type == "global") {
    if (family == "gaussian") {
      global_chars <- gaussian_char(clusterSolution = cluster_labels, data = data, nclusters = nclusters, type = type, signif = signif)
    }

    if (family == "negative.binomial") {
      global_chars <- negbin_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    if (family == "poisson") {
      global_chars <- poisson_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    if (family == "binomial") {
      global_chars <- binomial_char(clusterSolution = cluster_labels, data = data, K=K, type = type, signif = signif)
    }

    if (family == "ordinal") {
      global_chars <- ordinal_char(clusterSolution = cluster_labels, data = data, type = type, signif = signif)
    }

    # attributes/class for global_chars
    attr(global_chars, "family") <- family
    attr(global_chars, "K") <- K
    attr(global_chars, "type") <- type
    class(global_chars) <- c("globalchar","data.frame")

    # return data frame ready for printing
    return(global_chars)
  }
}

