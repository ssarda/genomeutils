#' @title Gene Ontology Analysis
#' @description This function conducts GO analysis for a set of selevcted genes.
#' The background genes can also be changed.
#'
#' @param all_genes A character vector with all genes
#' @param select_genes A character vector with a given set of interesting genes  
#' @return A topGO GenTable object with fold enrichment, pvalue, adjusted p-value and FDR. 
#' @export
GO_analysis = function(all_genes, select_genes)
{
	GO_genes <- rep(1,length(all_genes))
	names(GO_genes) = all_genes
	GO_genes[names(GO_genes) %in% select_genes] = 0

	GOdata <- topGO::new("topGOdata", ontology = "BP", allGenes = GO_genes, geneSel = function(p) p < 1e-2, description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

	resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")

	outtable = topGO::GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
	return(outtable)
}



#' @title Fetch gene lengths from gene IDs
#' @description Using the Ensembl hg19 assembly, the function established a SQL connection with the Ensembl DB, and retrieves all the transcripts associated with each gene.
#' A longest transcript is determined by merging all the exons of individual transcripts. 
#' The concatenated length is computed and returned per gene ID
#' @param select_genes A character vector with a given set of interesting genes  
#' @return A vector of gene lengths 
#' @export
fetch_gene_length = function(select_genes)
{
	edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
	get_longest_tx_length = function(x) {max(ensembldb::lengthOf(edb, of = "tx", filter = list(ensembldb::GeneidFilter(x))))}
	gene_lengths = sapply(select_genes, get_longest_tx_length)
	return(gene_lengths)
}




#' @title Fetch GC content from gene IDs
#' @description A Genomic Ranges object is constructed from  UCSC Genome Browser, identified by the Ensembl hg19 assembly gene IDs.
#' The function then retrieves the sequences of individual genes identified by the corrdinates of the gene IDs from BSgenome hg19 assembly.
#' The alphabet frequencies are computed, and then the G+C content is reported as a fraction. This is returned per gene ID
#' @param select_genes A character vector with a given set of interesting genes  
#' @return A vector of gene GC contents 
#' @export
fetch_gc_content = function(select_genes)
{
	txdb.Hs = GenomicFeatures::makeTxDbFromUCSC(genome="hg19",tablename = "ensGene")
	genes_features <- GenomicFeatures::genes(txdb.Hs, filter=list(gene_id=rownames(select_genes)))

	seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, genes_features)

	alfs <- Biostrings::alphabetFrequency(seqs, collapse=F)

	get_gc = function(fs) {(sum(fs['G'], fs['C']) / sum(fs['A'], fs['T'], fs['G'], fs['C'])) * 100}
	
	gcs = apply(alfs, 1, get_gc)
	return(gcs)
}