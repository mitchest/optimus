## Tests for find_optimal() function

# Load packages
library(testthat)
library(optimus)

context("Testing find_optimal() functionality")

# set up
cutreeLevels = 2:10


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


## count tests
counts_optimal <- find_optimal(data = counts, clustering = counts_cutree, family = "poisson", cutreeLevels = cutreeLevels)

test_that("find_optimal() returns a data frame", {
  expect_equal(class(counts_optimal), "aicsums")
  expect_equal(class(counts_optimal), "data.frame")
})

test_that("find_optimal() returns correct number of rows and columns", {
  expect_equal(ncol(counts_optimal), 2)
  expect_equal(nrow(counts_optimal), length(cutreeLevels))
})

test_that("find_optimal() returns correctly named columns", {
  expect_equal(names(counts_optimal) %in% c("sum_aic", "nclusters"), TRUE)
})
