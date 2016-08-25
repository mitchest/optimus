## Tests for find_optimal() function

# Load packages
library(testthat)
library(optimus)

context("Testing find_optimal() functionality")

# set up
cutreeLevels = 2:10


# mixed data
mixed <- data.frame(a = rpois(n = 42, lambda = 10),
                     b = rpois(n = 42, lambda = 11),
                     c = rpois(n = 42, lambda = 12),
                     d = rpois(n = 42, lambda = 13),
                     e = rep("char", 42),
                     f = as.factor(rep("char", 42)),
                    stringsAsFactors = FALSE)

# dummy counts for tests
counts <- data.frame(a = rpois(n = 42, lambda = 10),
                     b = rpois(n = 42, lambda = 11),
                     c = rpois(n = 42, lambda = 12),
                     d = rpois(n = 42, lambda = 13),
                     e = rpois(n = 42, lambda = 14),
                     f = rpois(n = 42, lambda = 15))
# cluster object that cutree works with
counts_cutree <- hclust(dist(counts))
# list of cluster solutions
counts_list <- lapply(2:10, FUN = function(x, data){cutree(tree = data, k = x)}, data=counts_cutree)


## generic input tests
test_that("find_optimal() stops with wrong data or arguments", {
  # wrong family spec
  expect_error(expr = optimus::find_optimal(data = mixed, clustering = counts_cutree, family = "poisson"), message = "some.*")
  # not cutree with cutree=T
  # not list with cutree=F
  # uneven cluster list lengths
  # nrow data matches length clustering labels
  # non numeric vector supplied to cutreeLevels
})


## count tests
counts_optimal <- optimus::find_optimal(data = counts, clustering = counts_cutree, family = "poisson", cutreeLevels = cutreeLevels)

test_that("find_optimal() returns correct class", {
  expect_equal(class(counts_optimal)[1], "aicsums")
  expect_equal(class(counts_optimal)[2], "data.frame")
})

test_that("find_optimal() returns correct number of rows and columns", {
  expect_equal(ncol(counts_optimal), 2)
  expect_equal(nrow(counts_optimal), length(cutreeLevels))
})

test_that("find_optimal() returns correctly named columns", {
  expect_equal(all(names(counts_optimal) %in% c("sum_aic", "nclusters")), TRUE)
})

# test AIC results are the same if cutree=T and cutree=F generating the same clustering?

