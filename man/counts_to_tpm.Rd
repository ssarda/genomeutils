% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_factorize_matrix.R
\name{counts_to_tpm}
\alias{counts_to_tpm}
\title{Convert counts to transcripts per million (TPM).}
\usage{
counts_to_tpm(counts, featureLength, meanFragmentLength)
}
\arguments{
\item{counts}{A numeric matrix of raw feature counts i.e.
fragments assigned to each gene.}

\item{featureLength}{A numeric vector with feature lengths.}

\item{meanFragmentLength}{A numeric vector with mean fragment lengths.}
}
\value{
tpm A numeric matrix normalized by library size and feature length.
}
\description{
Convert a numeric matrix of features (rows) and conditions (columns) with
raw feature counts to transcripts per million.
}
