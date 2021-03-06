% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_bayes_inference.R
\name{mle_weibull_dist}
\alias{mle_weibull_dist}
\title{Weibull MLE parameters and Model selection}
\usage{
mle_weibull_dist(X, start_k, start_lambda)
}
\arguments{
\item{X}{A Vector of values.}

\item{start_k}{Initial value to begin optimization at for shape of weibull distribution}

\item{start_lambda}{Initial value to begin optimization at for sscale of weibull distribution}
}
\value{
A vector containing the optimized MLE mean and standard deviation as well as the log likelihood and AIC to enable downstream model selection
}
\description{
This function accepts a vector with all values of a given unknown distribution
It estimates the maximum likelihood parameters by fitting a weibull distribution over it.
General-purpose optimization to estimate parameters is based on Nelder–Mead algorithm
The Akaike Information Criterion is also output to enable downstream model selection
}
