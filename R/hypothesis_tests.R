#' @title Hypothesis tests 
#'
#' @description Carries out hypotheses tests testing for significant differences in matched columns
#' of two matrices, and returns the most significant to the least significant result in order.
#' @param x Matrix 1
#' @param y Matrix 2
#' @param tests A vector with one or more of the following values: 1) t (t-test), 2) f (f-test), 3) mean (difference of means), 4) var (difference of variance, 5) wilcox (wilcoxon test), 6) ks (KS test) 
#' @return A data frame with the indices resulting from ordering the vector of, one of 1) differences, 2) statistics, or 3) pvalues computed between matched columns of the two matrices, in decreasing order of magnitude.
#' @export
foo_tests = function(x,y,tests)
{
	df = data.frame()
	if ("t" %in% tests)
	{
		df$t = sort(sapply(1:dim(x)[1], function(j) t.test(x[j, ], y[j, ])$statistic ), decreasing=TRUE, index.return=TRUE)$ix
	}
	if ("f" %in% tests)
	{
		df$f = sort(sapply(1:dim(x)[1], function(j) var.test(x[j, ], y[j, ])$statistic ), decreasing=TRUE, index.return=TRUE)$ix
	}
	if ("mean" %in% tests)
	{
		df$mean = sort(sapply(1:dim(x)[1], function(j) mean(x[j, ]) - mean(y[j, ]) ), decreasing=TRUE, index.return=TRUE)$ix
	}
	if ("var" %in% tests)
	{
		df$var = sort(sapply(1:dim(x)[1], function(j) log2(var(x[j, ])) - log2(var(y[j, ])) ), decreasing=TRUE, index.return=TRUE)$ix
	}
	if ("wilcox" %in% tests)
	{
		df$wilcox = sort(sapply(1:dim(x)[1], function(j) wilcox.test(x[j, ], y[j, ])$statistic ), decreasing=TRUE, index.return=TRUE)$ix
	}
	if ("ks" %in% tests)
	{
		df$ks = sort(sapply(1:dim(x)[1], function(j) ks.test(x[j, ], y[j, ])$statistic ), decreasing=TRUE, index.return=TRUE)$ix
	}
	return (df)
}