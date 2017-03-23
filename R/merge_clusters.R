#' @title Iteratively merges clusters in a way that improves predictive performance
#'
#' @description \code{merge_clusters} takes a clustering solution, generates all possible pairwise combinations of clusters, fits models to each combination, and merges the pair with the lowest delta AIC. The process is repeated iteratively
#'
#' @param data a data frame (or object that can be coerced by \code{\link[base]{as.data.frame}} containing the "raw" multivariate data. This is not necessarily the data used by the clustering algorithm - it is the data on which you are testing the predictive ablity of the clustering solution.
#' @param clustering an initial clustering solution (to be interatively merged) for \code{data}, that is, a vector of cluster labels (that can be coerced by \code{\link[base]{as.factor}}). The number of cluster labels must match the number of rows of the object supplied in the \code{data} argument. The solution could for example come form a call to \code{\link[stats]{cutree}}, see Examples
#' @param family a character string denoting the error distribution to be used for model fitting. The options are similar to those in \code{\link[stats]{family}}, but are more limited - see Details.
#' @param n.iter the number of merging iterations to perform, by default it will merge down to 2 clusters
#' @param K number of trials in binomial regression. By default, K=1 for presence-absence data (with cloglog link)
#' @param quietly suppress messages during merging procedure
#'
#' @details \code{merge_clusters} is built on the premise that a \emph{good} clustering solution (i.e. a classification) should provide information about the composition and abundance of the multivariate data it is classifying. A natural way to formalize this is with a predictive model, where group membership (clusters) is the predictor, and the multivariate data (site by variables matrix) is the response. \code{merge_clusters} fits linear models to each pariwise combination of a given set of clusters, and calculates their delta sum-of-AIC (that is, to the corresponding null model). The smallest delta AIC is taken to be the cluster pair that is \emph{most} similar, so it is merged, and the process is repeated. Lyons et al. (2016) provides background, a detailed description of the methodology, and application of delta AIC on both real and simulated ecological multivariate abundance data.
#'
#' At present, \code{merge_clusters} supports the following error distributions for model fitting:
#' \itemize{
#'   \item Gaussian (LM)
#'   \item Negative Binomial (GLM with log link)
#'   \item Poisson (GLM with log link)
#'   \item Binomial (GLM with cloglog link for binary daya, logit link otherwise)
#'   \item Ordinal (Proportional odds model with logit link)
#' }
#'
#' Gaussian LMs should be used for 'normal' data. Negative Binomial and Poisson GLMs shold be used for count data. Binomial GLMs should be used for binary and presence/absence data (when \code{K=1}), or trials data (e.g. frequency scores). If Binomial regression is being used with \code{K>1}, then \code{data} should be numerical values between 0 and 1, interpreted as the proportion of successful cases, where the total number of cases is given by \code{K} (see Details in \code{\link[stats]{family}}). Ordinal regression should be used for ordinal data, for example, cover-abundance scores. For ordinal regression, data should be supplied as either 1) factors, with the appropriate ordinal level order specified (see \code{\link[base]{levels}}) or 2) numeric, which will be coerced into a factor with levels ordered in numerical order (e.g. cover-abundance/numeric response scores). LMs fit via \code{\link[mvabund]{manylm}}; GLMs fit via \code{\link[mvabund]{manyglm}}; proportional odds model fit via \code{\link[ordinal]{clm}}.
#'
#' @return a list containing the clustering solution (vector) at each merge iteration. The object is of class \code{dsumaic}, and can be directly passed to \code{\link[optimus]{find_optimal}}.
#'
#' Attributes for the data frame are:
#'
#' \describe{
#'   \item{\code{family}}{ which error distribution was used for modelling, see Arguments}
#'   \item{\code{K}}{ number of cases for Binomial regression, see Arguments}
#' }
#'
#' @author Mitchell Lyons
#'
#' @references Lyons et al. 2016. Model-based assessment of ecological community classifications. \emph{Journal of Vegetation Science}, \strong{27 (4)}: 704--715.
#'
#' @seealso \code{\link[optimus]{find_optimal}}, \code{\link[optimus]{get_characteristic}}, S3 print function for 'daic' class, S3 residual plotting function
#'
#' @keywords merging, reallocation, iterative, pairwise
#'
#' @examples
#'
#' \dontrun{
#' ## Prep the 'swamps' data
#' ## ======================
#'
#' data(swamps) # see ?swamps
#' swamps <- swamps[,-1]
#'
#' ## Merge via AIC and compare to hclust heirarchy
#' ## =============================================
#'
#' ## perhaps not the best clustering option, but this is base R
#' swamps_hclust <- hclust(d = dist(x = log1p(swamps), method = "canberra"),
#'                        method = "complete")
#'
#' ## generate iteratively merged clustering solutions, based on sum-of-AIC
#' clustering_aicmerge <- merge_clusters(swamps, cutree(tree = swamps_hclust, k = 30),
#' family = "poisson", n.iter = 20)
#'
#' ## compare to hclust heirarchy
#' optimal_aicmerge <- find_optimal(data = swamps, clustering = clustering_aicmerge,
#' family = "poisson")
#'
#' optimal_hclust <- find_optimal(data = swamps, clustering = swamps_hclust,
#' family = "poisson", cutreeLevels = 10:30))
#'
#' plot(optimal_aicmerge)
#' points(optimal_hclust, col = "red", pch = 16)
#' }
#'
#' @export


# to begin merge based on overall sum-of-aic, but in future add option of merging based on 'rank' (proportion of variables with significant daic)

merge_clusters <- function(data, clustering, family, n.iter = NULL, K = 1, quietly = FALSE) {
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

  # set n.iter to default
  if (is.null(n.iter)) {n.iter <- nclusters - 2}
  if (n.iter > (nclusters - 2)) {
    n.iter <- nclusters - 2
    message("More iterations specified than possible, defaulting to merge down to 2 clusters")
  }

  # pre-allocate and add in initial clustering
  clustering_list <- vector(mode = "list", length = (n.iter + 1)) #init list
  clustering_list[[1]] <- cluster_labels

  # begin mering process
  if (quietly) {(message("Starting merging, no progress messages will be given"))}
  for (i in 1:n.iter) {
    if (!quietly) {message("Starting merging iteration ", i, ":")}

    pairs <- generate_pairwise(cluster_labels)

    fits <- pairwise_fits(data, cluster_labels, pairs, family, K)

    merge.pair <- c(fits$target[fits$dAIC_sum == min(fits$dAIC_sum)],
                    fits$test[fits$dAIC_sum == min(fits$dAIC_sum)])

    cluster_labels[cluster_labels %in% merge.pair] <- (9000 + i)

    if (!quietly) {message("merged cluster: ", merge.pair[1], " and ", merge.pair[2])}

    clustering_list[[(i+1)]] <- cluster_labels
  }

  # attributes/class for clustering_list
  attr(clustering_list, "family") <- family
  if (family == "binomial") {attr(aic_sums, "K") <- K}
  attr(clustering_list, "n.iter") <- n.iter
  class(clustering_list) <- c("dsumaic","list")

  # return data frame ready for plotting
  clustering_list
}

