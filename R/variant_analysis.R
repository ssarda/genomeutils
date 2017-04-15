#' @title GWAS
#' @description Conducts GWAS (per SNP association) for given genotype-phenotype data.
#' @param genodata Specify the entire genotype data object in SnpMatrix format
#' @param phenodata Data frame with a column of sample IDs, corresponding to the row names of genodata, and a columns for the continuous outcome variable. These columns must be named "id" and "phenotype".
#' @param family Specifies the link function for GLM. Default is guassian
#' @param filename Output filename with path
#' @param append LOGICAL whether to append output
#' @param workers Number of cores to use for parallel processing of function
#' @param flip In order to fit the model, genotype data is converted to numeric format using the as function from snpStats. The genotypes of each SNP are then coded as continuous, thereby taking on the value of 0, 1, and 2. To indicate the number of minor alleles instead, a flip.matrix procedure is included, which can be turned on or off using the flip argument.
#' @param select.snps Subset the analysis via a vector of SNP IDs
#' @param nSplits Number of SNP-wise splits that are made to the genotype data. The function runs the GWA analysis on these smaller subsets of the genotype data one at a time. After each subset has finished running the function will print a progress update onto the R console. By default this is set to 10
#' @return Outputs results in file with columns: "SNP", "Estimate", "Std.Error", "t-value", "p-value"
#' @export
GWAS <- function(genodata=genotypes,  phenodata=phenotypes, family = gaussian, filename=NULL, append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE, select.snps=NULL, nSplits=10)
{
    #Check that a filename was specified
    if(is.null(filename)) stop("Must specify a filename for output.")

    #Check that the genotype data is of class 'SnpMatrix'
    if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")

    #Check that there is a variable named 'phenotype' in phenodata table
    if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")

    #Check that there is a variable named 'id' in phenodata table
    if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")

    #If a vector of SNPs is given, subset genotype data for these SNPs
    if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata) %in% select.snps)]

    #Check that there are still SNPs in 'SnpMatrix' object
    if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")

    #Print the number of SNPs to be checked
    cat(paste(ncol(genodata), " SNPs included in analysis.\n"))

    #If append=FALSE than we will overwrite file with column names
    if(!isTRUE(append))
	{
        columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
        write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }

    # Check sample counts
    if (nrow(phenodata) != nrow(genodata))
	{
        warning("Number of samples mismatch. Using subset found in phenodata.")
    }

    # Order genodata rows to be the same as phenodata
    genodata <- genodata[phenodata$id,]

    cat(nrow(genodata), "samples included in analysis.\n")

    # Change which allele is counted (major or minor)
    flip.matrix<-function(x)
	{
        zero2 <- which(x==0)
        two0 <- which(x==2)
        x[zero2] <- 2
        x[two0] <- 0
        return(x)
    }

    nSNPs <- ncol(genodata)
    genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset

    snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
    snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

    cl <- doParallel::makeCluster(workers)

    show(cl)                    
    doParallel::registerDoParallel(cl)

    doParallel::foreach (part=1:nSplits) %do% {
        # Returns a standard matrix of the alleles encoded as 0, 1 or 2
        genoNum <- snpStats::as(genodata[,snp.start[part]:snp.stop[part]], "numeric")

        # Flip the numeric values of genotypes to count minor allele
        if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)

        # For each SNP, concatenate the genotype column to the phenodata and fit a generalized linear model
        rsVec <- colnames(genoNum)
        res <- doParallel::foreach(snp.name=rsVec, .combine='rbind') %dopar% {
            a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
            a$coefficients['snp',]
        }

        # write results so far to a file
        write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)

        cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
    }

    doParallel::stopCluster(cl)

    return(print("Done."))
}