#' @title Minimum maximum scaling of vector
#'
#' @description This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample
#'
#' @param x A Vector
#' @return A Vector scaled by maximum value.
#' @export
min_max_scale = function(x)
{

    if(max(x) == min(x))
        rep(.5, length(x))
    else
        sapply(1:length(x), function(j) (x[j] - min(x)) / (max(x)-min(x)) )
}



#' @title Row normalization
#'
#' @description This function loads a file as a matrix. It standardizes by substracting the median row value and dividing by the median absolute deviation.
#'
#' @param x A Matrix
#' @return The Scaled Matrix
#' @export
row_scale <- function(z)
{
rowmed <- apply(z, 1, median)
rowmad <- apply(z, 1, mad)  # median absolute deviation
rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
return(rv)
}



#' @title Column normalization
#'
#' @description This function standardizes only certain columns of an existing dataframe (specified by the user).
#' Standardization involves scaling columnwise.
#' @param dat A DataFrame
#' @param select A vector of selected column names for normalization
#' @return A DataFrame just with selected columns standardized
#' @export
col_scale <- function(dat, select)
{
	dat %>% dplyr::mutate_each_(dplyr::funs(scale(.) %>% as.vector), vars=select)
}


#' @title Upper Quartile Normalization 
#' @description Each column is divided by the 75% quantile of the counts for each library.
#' The calculated quantile is then scaled by the median across cells to keep the absolute level of expression relatively consistent.
#' Spike-ins should be excluded from the calculation of total expression in order to correct for total cell RNA content, therefore they should be specified to exclude them from downstream analysis
#' @param expr_mat Input matrix with expression values for genes in the form (Genes x Samples)
#' @param spikes A vector of genes that were spike-ins in the original experiment
#' @return A matrix with columns upper quartile normalized.
#' @export
uq_scale <- function (expr_mat, spikes = NULL) 
{
    UQ <- function(x)
	{
        quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
}




#' @title Compute counts per million (CPM)
#' @description Convert raw counts to counts per million (CPM) by dividing each column by its total then multiplying by 1,000,000.
#' Spike-ins should be excluded from the calculation of total expression in order to correct for total cell RNA content, therefore they should be specified to exclude them from downstream analysis
#' @param expr_mat Input matrix with expression values for genes in the form (Genes x Samples)
#' @param spikes A vector of genes that were spike-ins in the original experiment
#' @return A matrix with CPM expression values
#' @export
counts_to_cpm <- function (expr_mat, spikes = NULL) 
{
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}


							 
#' @title Convert counts to transcripts per million (TPM).
#' 
#' @description Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
#' @export
counts_to_tpm <- function(counts, featureLength, meanFragmentLength)
{
  
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i){featureLength - meanFragmentLength[i] + 1}))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) 
  {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}




#' @title Size Factor normalization as used by DESeq for differential expression. Also called RLE (Relative Log Expression) 
#' @description First the geometric mean of each gene across all cells is calculated. The size factor for each cell is the median across genes of the ratio of the expression to the gene’s geometric mean.
#' Spike-ins should be excluded from the calculation of total expression in order to correct for total cell RNA content, therefore they should be specified to exclude them from downstream analysis
#' @param expr_mat Input matrix with expression values for genes in the form (Genes x Samples)
#' @param spikes A vector of genes that were spike-ins in the original experiment
#' @return A matrix with SF normalized RLE values
#' @export
counts_to_rle <- function (expr_mat, spikes = NULL) 
{
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
        median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
            0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}




#' @title Can the matrix be inverted?
#'
#' @description This function checks if target matrix is invertible
#'
#' @param m A Matrix
#' @return LOGICAL Indicating invertibility of matrix
#' @export
check_inverse <- function(m) class(try(qr.solve(m),silent=T))=="matrix"



#' @title Principal Component Analysis 2D plot
#'
#' @description This function accepts a matrix, computes the principal components using a kernel transformation and plots the first 2 eigen vectors
#'
#' @param X A Gene Matrix
#' @return A plot with the first two eigen vectors plotted on the X and Y axis
#' @export
plot_pca = function(X, ...)
{
    pca_data = t(X)
    pcomponents = kernlab::rotated(kernlab::kpca(pca_data, kernel=kernlab::vanilladot(), ...))
    plot(pcomponents[,1], pcomponents[,2])
}




#' @title Multi-Dimensional Scaling 2D plot
#'
#' @description This function accepts a matrix, computes the pairwise distance matrix, and projects the distances onto a lwer dimensional (2D space)
#'
#' @param X A Matrix
#' @return A 2D MDS plot
#' @export
plot_mds = function(X)
{
    mds_dist = dist(X)
    mds_fit = cmdscale(mds_dist,eig=TRUE, k=2)
    plot(mds_fit$points[,1], mds_fit$points[,2])
}



#' @title Singular Value Decomposition
#'
#' @description This function computes the singular value decomposition of a matrix and plots the first two eigen vectors.
#'
#' @param X A Matrix
#' @return A plot with the first two eigen vectors plotted on the X and Y axis
#' @export
plot_svd = function(X)
{
	xx <- svd(X %*% t(X))
	xxd <- xx$v %*% sqrt(diag(xx$d))
	x1 <- xxd[, 1]
	y1 <- xxd[, 2]
	plot(x1, y1)
}