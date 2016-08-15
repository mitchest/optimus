## Tests for find_optimal() function

## Load packages
library(testthat)
library(optimus)

context("Testing find_optimal() functionality")

## dummy data for tests
counts <- data.frame(a = rpois(n = 42, lambda = 10),
                     b = rpois(n = 42, lambda = 11),
                     c = rpois(n = 42, lambda = 12),
                     d = rpois(n = 42, lambda = 13),
                     e = rpois(n = 42, lambda = 14),
                     f = rpois(n = 42, lambda = 15))
