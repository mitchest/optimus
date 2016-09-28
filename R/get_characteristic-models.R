# hidden modelling functions for get_characteristic() ---------------------



# gaussian lms via mvabund
####### ---> need to check whether the RSS methods is OK for AIC for least squares out of manylm???? Looks like it overfits??
gaussian_char <- function(clusterSolution, data, nclusters, type) {
  data <- mvabund::mvabund(data)
  clusterSolution <- as.factor(clusterSolution)
  if (type == "per.cluster") {
    fit <- mvabund::manylm(formula = data ~ clusterSolution)
    fit_coefs <- as.data.frame(fit$coefficients)
    cluster_names <- stats::setNames(row.names(fit_coefs), row.names(fit_coefs))
    ret <- lapply(X = cluster_names, FUN = sort_char_coef, fit_coefs)
  }
  if (type == "global") {
  fit_rss <- apply(X = mvabund::manylm(formula = data ~ clusterSolution)$residuals^2,
                  MARGIN = 2, FUN = sum)
  fit_aic <- (2 * (nclusters + 2)) + (nrow(data) * log(fit_rss))
  null_rss <- apply(X = mvabund::manylm(formula = data ~ 1)$residuals^2,
                   MARGIN = 2, FUN = sum)
  null_aic <- (2 * (nclusters + 2)) + (nrow(data) * log(null_rss))
  daic <- data.frame(sort(null_aic - fit_aic, decreasing = TRUE))
  data.frame(variables = row.names(daic),
             dAIC = daic$sort)
  }
  ret
}

# negative binomal glms via mvabund
negbin_char <- function(clusterSolution, data, type) {
  data <- mvabund::mvabund(data)
  clusterSolution <- as.factor(clusterSolution)
  fit <- mvabund::manyglm(formula = data ~ clusterSolution-1, family = "negative.binomial") # -1 for means parameterisaiton
  if (type == "per.cluster") {
    fit_coefs <- as.data.frame(fit$coefficients)
    cluster_names <- stats::setNames(row.names(fit_coefs), row.names(fit_coefs))
    ret <- lapply(X = cluster_names, FUN = sort_char_coef, fit_coefs)
  }
  if (type == "global") {
    fit_null <- mvabund::manyglm(formula = data ~ 1, family = "negative.binomial")
    daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
    ret <- data.frame(variables = row.names(daic),
                      dAIC = daic$sort)
  }
  ret
}

# poisson glms via mvabund
poisson_char <- function(clusterSolution, data, type) {
  data <- mvabund::mvabund(data)
  clusterSolution <- as.factor(clusterSolution)
  fit <- mvabund::manyglm(formula = data ~ clusterSolution-1, family = "poisson")
  if (type == "per.cluster") {
    fit_coefs <- as.data.frame(fit$coefficients)
    cluster_names <- stats::setNames(row.names(fit_coefs), row.names(fit_coefs))
    ret <- lapply(X = cluster_names, FUN = sort_char_coef, fit_coefs)
  }
  if (type == "global") {
    fit_null <- mvabund::manyglm(formula = data ~ 1, family = "poisson")
    daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
    ret <- data.frame(variables = row.names(daic),
                      dAIC = daic$sort)
  }
  ret
}

# binomal glms (K=1 for logistic regression) via mvabund
binomial_char <- function(clusterSolution, data, K, type) {
  data <- mvabund::mvabund(data)
  clusterSolution <- as.factor(clusterSolution)
  # cloglog link is better for pres/abs
  if (K == 1) {
    fit <- mvabund::manyglm(formula = data ~ clusterSolution-1, family = stats::binomial(link='cloglog'), K = K)
    if (type == "per.cluster") {
      fit_coefs <- as.data.frame(fit$coefficients)
      cluster_names <- stats::setNames(row.names(fit_coefs), row.names(fit_coefs))
      ret <- lapply(X = cluster_names, FUN = sort_char_coef, fit_coefs)
    }
    if (type == "global") {
      fit_null <- mvabund::manyglm(formula = data ~ 1, family = stats::binomial(link='cloglog'), K = K)
      daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
      ret <- data.frame(variables = row.names(daic),
                        dAIC = daic$sort)
    }
  # logit link is better proportions/trials
  } else {
    fit <- mvabund::manyglm(formula = data ~ clusterSolution-1, family = stats::binomial(link='logit'), K = K)
    if (type == "per.cluster") {
      fit_coefs <- as.data.frame(fit$coefficients)
      cluster_names <- stats::setNames(row.names(fit_coefs), row.names(fit_coefs))
      ret <- lapply(X = cluster_names, FUN = sort_char_coef, fit_coefs)
    }
    if (type == "global") {
      fit_null <- mvabund::manyglm(formula = data ~ 1, family = stats::binomial(link='logit'), K = K)
      daic <- data.frame(sort(fit_null$aic - fit$aic, decreasing = TRUE))
      ret <- data.frame(variables = row.names(daic),
                        dAIC = daic$sort)
    }
  }
  ret
}

# ordinal regression via clm
ordinal_char <- function(clusterSolution, data, type) {
  data <- data.frame(lapply(data, as.factor))
  clusterSolution <- as.factor(clusterSolution)
  # per cluster coeffs are not well defined for cumulative link models
  # there's a threshold coefficient for each level of the ordinal variable, and then a coef for all the remaining covariate levels
  # thus hard to define the effect of each level with one coefficient value
  # reverting to pres/abs model until the best approach is decided
  if (type == "per.cluster") {
    message("Per-cluster characteristic variables not well defined for cumulative link models.")
    message("Reverting to logistic regression for the mean time - first level of the factor will be used as the threshold. If that's not suitable, recode to binary youeself, then choose the binomial family")
    data = data.frame(lapply(X = data, FUN = ordinal_to_binom))
    ret <- binomial_char(clusterSolution = clusterSolution, data = data, K = 1, type = "per.cluster")
  }
  if (type == "global") {
    fit <- manyclm(responses = data, predictor = clusterSolution)
    fit_null <- manyclm(responses = data, predictor = 1)
    daic <- data.frame(sort(fit_null - fit, decreasing = TRUE))
    ret <- data.frame(variables = row.names(daic),
                      dAIC = daic$sort)
  }
  ret
}
