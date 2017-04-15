#' @title Load a File as Matrix
#'
#' @description This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being kept.
#'
#' @param infile Path to the input file
#' @return A matrix of the input file
#' @export
fastload_mat <- function(infile)
{
	in.dt <- data.table::fread(infile, header = TRUE)
	in.dt <- in.dt[!duplicated(in.dt[, 1]), ]
	in.mat <- as.matrix(in.dt[, -1, with = FALSE])
	rownames(in.mat) <- unlist(in.dt[, 1, with = FALSE])
	in.mat
}




#' @title Write a New File as Tab Delimited
#'
#' @description This function writes to a file. It assumes that the object to be written #' is a data.table object. It handles lists or matrix formats.
#' Column identifiers are retained.
#' @param outmat Path to the output object
#' @return NULL
#' @export
fastwrite_tab <- function(outmat, outfile)
{
	data.table::fwrite(outmat, file = outfile, append = FALSE, quote = "auto", sep = "\t", sep2 = c("","|",""), eol = if (.Platform$OS.type=="windows") "\r\n" else "\n", na = "NA", row.names = FALSE, col.names = TRUE, logicalAsInt = FALSE, nThread = getDTthreads())
}




#' @title Read the genome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the [Biostrings]{readDNAStringSet} function.
#' @return A data.table storing the gene id in the first column and the corresponding sequence as string in the second column.
#' @export
read_genome <- function(file, format, ...)
{        
	if(!is.element(format,c("fasta","fastq")))
		stop("Please choose a file format that is supported by this function.")
	geneids <- seqs <- NULL
    if(format == "fasta")
	{
		tryCatch ({
		genome <- Biostrings::readDNAStringSet(filepath = file, format = format, ...)
		genome_names <- as.vector(unlist(sapply(genome@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
        genome.dt <- data.table::data.table(geneids = genome_names,seqs = toupper(as.character(genome)))
		data.table::setkey(genome.dt,geneids)}, error = function(e) {stop(paste0("File ",file, " could not be read properly. \n","Please make sure that ",file," contains only DNA sequences and is in ",format," format."))})
     }
     return(genome.dt)
}




#' @title Read the proteome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the [Biostrings]{readDNAStringSet} function.
#' @return A data.table storing the gene id in the first column and the corresponding sequence as string in the second column.
#' @export
read_proteome <- function(file, format, ...)
{        
	if(!is.element(format,c("fasta","fastq")))
		stop("Please choose a file format that is supported by this function.")
	geneids <- seqs <- NULL
	tryCatch ({ 
	   proteome <- Biostrings::readAAStringSet(filepath = file, format = format, ...)
	   proteome_names <- as.vector(unlist(sapply(proteome@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
	   proteome.dt <- data.table::data.table(geneids = proteome_names,seqs =toupper(as.character(proteome)))
	   data.table::setkey(proteome.dt,geneids)}, error = function(e) { stop(paste0("File ",file, " could not be read properly. \n","Please make sure that ",file," contains only amino acid sequences and is in ",format," format."))})
     return(proteome.dt)
}





#' @title Save a genome in fasta format
#'
#' This function writes to a file. It assumes that the object to be written #' is a data.table object.
#' @param genome A data.table object storing DNA sequences in the first column and gene ids in the second column.
#' @param file.name A character string specifying the name of the fasta output file
#' @param as.string When set to TRUE sequences are in the form of strings instead of vectors of single characters.
#' @return NULL
#' @export
write_genome <- function(genome, file.name, nbchar = 80, open = "w", as.string = TRUE)
{
	seqs <- geneids <- NULL
    seqinr::write.fasta(sequences = as.list(genome[,seqs]),names = genome[,geneids],nbchar = nbchar, open = open, file.out = file.name,as.string = as.string)
}




#' @title Save a proteome in fasta format
#'
#' This function writes to a file. It assumes that the object to be written #' is a data.table object.
#' @param proteome A data.table object storing amino acid sequences in the first column and gene ids in the second column.
#' @param file.name A character string specifying the name of the fasta output file
#' @param as.string When set to TRUE sequences are in the form of strings instead of vectors of single characters.
#' @return NULL
#' @export
write_proteome <- function(proteome, file.name, nbchar = 80, open = "w", as.string = TRUE)
{
	seqs <- geneids <- NULL
    seqinr::write.fasta(sequences = as.list(proteome[,seqs]),names  = proteome[,geneids],nbchar = nbchar, open = open, file.out = file.name,as.string = as.string)
}							 