# hidden utility functions ------------------------------------------------

# Determine if range of vector is FP 0.
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

# manyglm analogy for prop odds model
manyclm <- function(responses, predictor) {
  if (length(predictor)==1) { # null model
    aics <- unlist(lapply(X = responses,
                          FUN = function(x) {stats::AIC(ordinal::clm(formula = x ~ 1))}))
    variables <- unlist(lapply(X = responses, FUN = names))
  } else {
    aics <- unlist(lapply(X = responses,
                          FUN = function(x, predictor) {stats::AIC(ordinal::clm(formula = x ~ predictor))},
                          predictor = predictor))
    variables <- unlist(lapply(X = responses, FUN = names))
  }
  names(aics) <- variables
  aics
}
# manyglm analogy for prop odds model - no variable names
manyclm_naked <- function(responses, predictor) {
  if (length(predictor)==1) { # null model
    aics <- unlist(lapply(X = responses,
                          FUN = function(x) {stats::AIC(ordinal::clm(formula = x ~ 1))}))
  } else {
    aics <- unlist(lapply(X = responses,
                          FUN = function(x, predictor) {stats::AIC(ordinal::clm(formula = x ~ predictor))},
                          predictor = predictor))
  }
  aics
}

# $AICsum analogy form manyglm for prop odds model
manyclm_sum <- function(responses, clusters) {
  sum_clm_aics <- function(x, clusters) {stats::AIC(ordinal::clm(formula = x ~ clusters))}
  sum(unlist(lapply(X = responses, FUN = sum_clm_aics, clusters = clusters)))
}

# rank variable coefficients
sort_char_coef <- function(x, coefs) {
  dat <- setNames(as.numeric(coefs[x, ]), names(coefs[x, ]))
  sorted_coefs <- sort(dat[dat >= 0], decreasing = TRUE)
  data.frame(variables = names(sorted_coefs),
             coef_value = sorted_coefs,
             stringsAsFactors = FALSE)
}

# turn ordinal data into binary (based on first level of factor)
ordinal_to_binom <- function(x) {
  ifelse(x == levels(x)[1], 1, 0)
}
