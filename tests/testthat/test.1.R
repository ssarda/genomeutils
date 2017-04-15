library(genomeutils)
context("Scaling length preservation")

test_that("Scaling performs as expected on an identical vector", {
  expect_equal(min_max_scale(c(1,1,1)), c(.5,.5,.5))
})

test_that("Scaling preserves length", {
  expect_equal(length(c(1,1,1)), length(min_max_scale(c(1,1,1))))
})

