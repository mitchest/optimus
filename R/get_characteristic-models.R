# hidden modelling functions for get_characteristic() ---------------------



# gaussian lms via mvabund
####### ---> need to check whether the RSS methods is OK for AIC for least squares out of manylm???? Looks like it overfits??
gaussian_char <- function(clusterSolution, data, nclusters) {
  data <- mvabund::mvabund(data)
  fit_rss <- apply(X = mvabund::manylm(formula = data ~ as.factor(clusterSolution))$residuals^2,
                  MARGIN = 2, FUN = sum)
  fit_aic <- (2 * (nclusters + 2)) + (nrow(data) * log(fit_rss))
  null_rss <- apply(X = mvabund::manylm(formula = data ~ 1)$residuals^2,
                   MARGIN = 2, FUN = sum)
  null_aic <- (2 * (nclusters + 2)) + (nrow(data) * log(null_rss))
  daic <- data.frame(sort(null_aic - fit_aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
}

# negative binomal glms via mvabund
negbin_char <- function(clusterSolution, data) {
  data <- mvabund::mvabund(data)
  fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = "negative.binomial")
  fit_null <- mvabund::manyglm(formula = data ~ 1, family = "negative.binomial")
  daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
}

# poisson glms via mvabund
poisson_char <- function(clusterSolution, data) {
  data <- mvabund::mvabund(data)
  fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = "poisson")
  fit_null <- mvabund::manyglm(formula = data ~ 1, family = "poisson")
  daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
}

# binomal glms (K=1 for logistic regression) via mvabund
binomial_char <- function(clusterSolution, data, K) {
  data <- mvabund::mvabund(data)
  if (K == 1) { # cloglog link is better for pres/abs
    fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = stats::binomial(link='logit'), K = K)
    fit_null <- mvabund::manyglm(formula = data ~ 1, family = stats::binomial(link='cloglog'), K = K)
    daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  } else {
    fit <- mvabund::manyglm(formula = data ~ as.factor(clusterSolution), family = stats::binomial(link='logit'), K = K)
    fit_null <- mvabund::manyglm(formula = data ~ 1, family = stats::binomial(link='logit'), K = K)
    daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
  }
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
}

# ordinal regression via clm
ordinal_char <- function(clusterSolution, data) {
  data <- data.frame(lapply(data, as.factor))
  fit <- manyclm(responses = data, predictor = as.factor(clusterSolution))
  fit_null <- manyclm(responses = data, predictor = 1)
  daic <- data.frame(sort(fit_null - fit, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
}

