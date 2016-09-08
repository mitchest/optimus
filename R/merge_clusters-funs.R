

# # fit and get AIC - move this to Examples i think
#
# alloc.list = list()
# for (i in 1:length(2:49)) {
#   alloc.list[[i]] = get(paste0("alloc.test.",c(2:49)[i]))
# }
# save(alloc.list, file="alloc.list.RData")
#
# data = mvabund(pa)
# alloc.sumAIC = numeric(length(alloc.list))
# for (i in 1:length(alloc.list)) {
#   alloc.sumAIC[i] = manyglm(formula=data~as.factor(alloc.list[[i]]), family="binomial")$AICsum
# }
#
# alloc.pairwise.AIC = data.frame(AICsum=alloc.sumAIC, groups=c(2:49))
# load("pa.AICsum.RData")
# pa.AICsum.pairwise = rbind(pa.AICsum, alloc.pairwise.AIC)
# save(pa.AICsum.pairwise, file="pa.AICsum.pairwise.RData")
# load("pa.AICsum.pairwise.RData")
#
#
# # plot pairwise modelling merges
# AIC.hclust = pa.AICsum.pairwise[1:49,]
# AIC.pairwise = pa.AICsum.pairwise[50:97,]




# FUNCTIONS ---------------------------------------------------------------


generate_pairwise <- function(cluster_labels) {
  # generate pairwise tests
  unique_labs <- as.character(unique(cluster_labels))
  target <- character(round(length(unique_labs)^2 / 2))
  test <- character(round(length(unique_labs)^2 / 2))
  # loop across unique_labs to create pairwise comparisons, without reverse tests
  ticker <- 0
  for (itarget in 1:(length(unique_labs) - 1)) {
    for (itest in itarget:(length(unique_labs) - 1)) {
      ticker <- (ticker + 1)
      target[ticker] <- unique_labs[itarget]
      test[ticker] <- unique_labs[(itest+1)]
    }
  }
  target <- target[!target == ""]
  test <-  test[!test == ""]

  data.frame(target = target, test = test, stringsAsFactors = FALSE)
}


pairwise_fits <- function(data, cluster_labels, pairs, family, K) {
  target <- as.character(pairs$target)
  test <- as.character(pairs$test)
  nvars <- numeric(nrow(pairs)) # don't actually do anything with this yet...
  dAIC_sum <- numeric(nrow(pairs))
  #rank <- numeric(nrow(pairs))

  for (i in 1:length(target)) {
    # subset to target sites
    target_rows <- which(cluster_labels == target[i])
    target_data <- data[target_rows,]
    target_cluster <- cluster_labels[target_rows]
    # subset to test sites
    test_rows <- which(cluster_labels == test[i])
    test_data <- data[test_rows,]
    test_cluster <- cluster_labels[test_rows]
    # create data objects for the pariwise test
    combined_data <- as.data.frame(rbind(target_data, test_data))
    combined_cluster <- as.factor(c(as.character(target_cluster), as.character(test_cluster)))

    # fit models
    if (family == "gaussian") {
      nclusters <- length(unique(combined_cluster)) # not really important, since it's a constant, but leaving here for consistency
      combined_data <- mvabund::mvabund(combined_data)
      fit_rss <- apply(X = mvabund::manylm(formula = combined_data ~ as.factor(combined_cluster))$residuals^2,
                       MARGIN = 2, FUN = sum)
      fit_aic <- sum( (2 * (nclusters + 2)) + (nrow(combined_data) * log(fit_rss)) )
      null_rss <- apply(X = mvabund::manylm(formula = combined_data ~ 1)$residuals^2,
                        MARGIN = 2, FUN = sum)
      null_aic <- sum( (2 * (nclusters + 2)) + (nrow(combined_data) * log(null_rss)) )
    }

    if (family == "negative.binomial") {
      combined_data <- mvabund::mvabund(combined_data)
      fit_aic <- mvabund::manyglm(combined_data ~ combined_cluster, family="negative.binomial")$AICsum
      null_aic <- mvabund::manyglm(combined_data ~ 1, family="negative.binomial")$AICsum
    }

    if (family == "poisson") {
      combined_data <- mvabund::mvabund(combined_data)
      fit_aic <- mvabund::manyglm(combined_data ~ combined_cluster, family="poisson")$AICsum
      null_aic <- mvabund::manyglm(combined_data ~ 1, family="poisson")$AICsum
    }

    if (family == "binomial") {
      combined_data <- mvabund::mvabund(combined_data)
      if (K == 1) {
        fit_aic <- mvabund::manyglm(combined_data ~ combined_cluster, family=stats::binomial(link = "cloglog"))$AICsum
        null_aic <- mvabund::manyglm(combined_data ~ 1, family=stats::binomial(link = "cloglog"))$AICsum
      } else {
        fit_aic <- mvabund::manyglm(combined_data ~ combined_cluster, family=stats::binomial(link = "logit"))$AICsum
        null_aic <- mvabund::manyglm(combined_data ~ 1, family=stats::binomial(link = "logit"))$AICsum
      }
    }

    if (family == "ordinal") {
      combined_data <- data.frame(lapply(combined_data, as.factor))
      fit_aic <- sum( manyclm_naked(responses = combined_data, predictor = as.factor(combined_cluster)) )
      null_aic <- sum (manyclm_naked(responses = combined_data, predictor = 1) )
    }

    # store results
    nvars[i] <- ncol(combined_data)
    dAIC_sum[i] <- null_aic - fit_aic
    # # calculate % of species for which dAIC (from null) is >n
    # dAIC <- fit_null$aic-fit$aic
    # rank[i] <- sum(dAIC>4)/length(dAIC)
  }
  data.frame(target = pairs$target, test = pairs$test, #rank = rank
             nvars = nvars, dAIC_sum = dAIC_sum, stringsAsFactors=FALSE)
}


# ## function that calculates the species count after removing species with 0 occurance in a pariwise combination
# PairwiseCombinedSpeciesCount = function(data, cluster_labels, pairs) {
#   target = as.character(pairs$target)
#   test = as.character(pairs$test)
#   for (i in 1:length(target)){
#     # subset to target sites
#     target.rows = which(cluster_labels==target[i])
#     target.data = data[target.rows,]
#     target.cluster = cluster_labels[target.rows]
#     # subset to test sites
#     test.rows = which(cluster_labels==test[i])
#     test.data = data[test.rows,]
#     test.cluster = cluster_labels[test.rows]
#     # create data objects for the pariwise test
#     combined.data = rbind(target.data, test.data)
#     combined.data = combined.data[,colSums(combined.data)>0] # remove species that don't occur in either community
#     print(ncol(combined.data))
#   }
# }


# ## function that calculates lowest ranked pairwise comparison and lowest ranked cluster
# # alloc==allocation used, pairwisemanyglm=object from PairwiseManyglm() using alloc
# PairwiseSummary = function(alloc, pairwisemanyglm) {
#   # print cluster numbers
#   print(table(alloc))
#   # find lowest ranked pairwise comparison
#   print(pairwisemanyglm[pairwisemanyglm$dAIC_sum==min(pairwisemanyglm$dAIC_sum),])
#   # find cluster that has lowest mean ranks
#   #   for (i in unique(alloc)) {
#   #     print(paste0("Cluster ",i," mean delta sum-of-AIC:"))
#   #     print(mean(pairwisemanyglm$dAIC_sum[pairwisemanyglm$target==i | pairwisemanyglm$test==i]))
#   #   }
# }
