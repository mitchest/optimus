# hidden modelling functions for find_optimal() -------------------------


# gaussian lms via mvabund
####### ---> need to check whether the RSS methods is OK for AIC for least squares out of manylm???? Looks like it overfits??
gaussian_loop <- function(cluster_list, data, cutreeLevels) {
  data <- mvabund::mvabund(data)
  get_rss <- function(x, data) {sum( ((mvabund::manylm(formula = data ~ as.factor(x))$residuals)^2) )}
  sum_rss <- unlist(lapply(X = cluster_list, FUN = get_rss, data = data))
  (2 * (cutreeLevels + 2)) + (nrow(data) * log(sum_rss))
}

# negative binomal glms via mvabund
negbin_loop <- function(cluster_list, data) {
  data <- mvabund::mvabund(data)
  sum_aic <- function(x, data) {mvabund::manyglm(formula = data ~ as.factor(x), family = "negative.binomial")$AICsum}
  unlist(lapply(X = cluster_list, FUN = sum_aic, data = data))
}

# poisson glms via mvabund
poisson_loop <- function(cluster_list, data) {
  data <- mvabund::mvabund(data)
  sum_aic <- function(x, data) {mvabund::manyglm(formula = data ~ as.factor(x), family = "poisson")$AICsum}
  unlist(lapply(X = cluster_list, FUN = sum_aic, data = data))
}

# binomal glms (K=1 for logistic regression) via mvabund
binomial_loop <- function(cluster_list, data, K) {
  data <- mvabund::mvabund(data)
  sum_aic <- function(x, data) {mvabund::manyglm(formula = data ~ as.factor(x), family = "binomial", K=K)$AICsum}
  unlist(lapply(X = cluster_list, FUN = sum_aic, data = data))
}

# ordinal regression via clm
ordinal_loop <- function(cluster_list, data) {
  data <- data.frame(lapply(data, as.factor))
  sum_aic <- function(x, data) {manyclm_sum(responses = data, clusters = as.factor(x))}
  unlist(lapply(X = cluster_list, FUN = sum_aic, data = data))
}

manyclm_sum <- function(responses, clusters) {
  sum_clm_aics <- function(x, clusters) {stats::AIC(ordinal::clm(formula = x ~ clusters))}
  sum(unlist(lapply(X = responses, FUN = sum_clm_aics, clusters = clusters)))
}



# more distributions to add -----------------------------------------------

# # cover
# #data = ((cover/100)*(nrow(cover)-1)+0.5)/nrow(cover) # this is a transformation for betareg
# data = cover/100
# # use zero/one inflated beta
# AICsum = numeric(length(groups))
# for (i in 1:length(groups)) {
#   alloc = cutree(clust, groups[i])
#   AICsum[i] = manybeta.AICsum(data, as.factor(alloc))
# }
#
# # use compund Poisson-gamma (i.e. tweedie with power= 1.5)
# AICsum = numeric(length(groups))
# for (i in 1:length(groups)) {
#   Xdata = data.frame(alloc=cutree(clust, groups[i]), dummy=integer(nrow(data)))
#   AICsum[i] = sum(AIC(manyany("glm", data, data~alloc, data=Xdata,
#                               family=tweedie(var.power=1.2, link.power=0), var.power=1.2, composition=FALSE)))
# }
