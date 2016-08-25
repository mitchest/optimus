# hidden modelling functions for get_characteristic() ---------------------



# gaussian lms via mvabund
####### ---> need to check whether the RSS methods is OK for AIC for least squares out of manylm???? Looks like it overfits??
gaussian_char <- function(clusterSolution, data, nclusters) {
  # data <- mvabund::mvabund(data)
  # get_rss <- function(x, data) {sum( ((mvabund::manylm(formula = data ~ as.factor(x))$residuals)^2) )}
  # sum_rss <- unlist(lapply(X = clusterSolution, FUN = get_rss, data = data))
  # (2 * (nclusters + 2)) + (nrow(data) * log(sum_rss))
  stop("Sorry, still working on Gaussian...")
}

# negative binomal glms via mvabund
negbin_char <- function(clusterSolution, data) {
  data <- mvabund::mvabund(data)
  fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = "negative.binomial")
  fit_null <- mvabund::manyglm(formula = data ~ 1, family = "negative.binomial")
  daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC=daic$sort)
}

# poisson glms via mvabund
poisson_char <- function(clusterSolution, data) {
  data <- mvabund::mvabund(data)
  fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = "poisson")
  fit_null <- mvabund::manyglm(formula = data ~ 1, family = "poisson")
  daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC=daic$sort)
}

# binomal glms (K=1 for logistic regression) via mvabund
binomial_char <- function(clusterSolution, data, K) {
  data <- mvabund::mvabund(data)
  # possibly change logistic regression to link='cloglog', but for the meantime use logit link
  fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = binomial(link='logit'), K = K)
  fit_null <- mvabund::manyglm(formula = data ~ 1, family = binomial(link='logit'), K = K)
  daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC=daic$sort)
}

# ordinal regression via clm
ordinal_char <- function(clusterSolution, data) {
  data <- data.frame(lapply(data, as.factor))
  fit <- manyclm(responses = data, predictor = as.factor(clusterSolution))
  fit_null <- manyclm(responses = data, predictor = 1)
  daic <- data.frame(sort(fit_null - fit, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC=daic$sort)
}

manyclm <- function(responses, predictor) {
  if (length(predictor)==1) { # null model
    aics <- unlist(lapply(X = responses,
                          FUN = function(x) {stats::AIC(ordinal::clm(x ~ 1))}))
    variables <- unlist(lapply(X = responses, FUN = names))
  } else {
    aics <- unlist(lapply(X = responses,
                          FUN = function(x, predictor) {stats::AIC(ordinal::clm(x ~ predictor))},
                          predictor = predictor))
    variables <- unlist(lapply(X = responses, FUN = names))
  }
  names(aics) = variables
  aics
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
