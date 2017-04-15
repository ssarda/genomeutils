library(genomeutils)
context("Scaling object type expectation")


test_that("Scaling and normalization functions return expected object type", {
  vec = c(1,2,5,6,7,8)
  mat = matrix(runif(5*5), ncol=5)
  df = as.data.frame(mat)
  expect_is(min_max_scale(vec), "numeric")  
  expect_is(row_scale(mat), "matrix")
  expect_is(col_scale(df, c("V1","V2")), "data.frame")
  expect_is(check_inverse(mat), "logical")
  expect_is(plot_pca(mat), "NULL")
  expect_is(plot_mds(mat), "NULL")
  expect_is(plot_svd(mat), "NULL")
})


