# genomeutils
An R helper package for routine tasks in genome analysis


This package provides a set of helper tools that allow automating some common tasks one encounters during routine genome analysis.

## Fast loading and writing of files and genome objects
*Using the data.table package for performing these functions really fast.* 

Load from a file   
Write to a file   
Read .fasta/fastq to a Genome object      
Write a Genome object to file   
Read .fasta to a Proteome object   
Write a Proteome object to file      


## Fetch gene attributes and functional analysis
   
Fetch the gene length (coding exons) for a list of genes   
Fetch the GC content of the sequence for a list of genes      
Perform Gene Ontology Enrichment analysis for a set of interesting genes      
   

## Normalization and Factorization methods for gene matrices
   
Min max normalization   
Row median/deviation normalization   
Sample specific normalization  
Upper Quartile Normalization
Gene counts to expression in 'Counts per million': (CPM)  
Gene counts to expression in 'Transcripts per million': (TPM)  
Gene counts to expression in 'Relative Log Expression': (RLE) as used in DESeq     
Invertibility of matrix   
Principal Component Analysis + 2D plots      
Multi-Dimensional Scaling + 2D plots   
Singular Value Decomposition + 2D plots
    
       
## Pretty plotting functions

Modify heirarchical clustering to produce a plot colored by groups     
Produce an MA plot (to identify differentially expressed genes, for instance)   
Produce fancy heatmaps   
Produce a smooth histogram by modifying base R plotting parameters   
   
   
## Hypothesis testing for matched groups of data
*Includes ordering comparisons by significance*   

T-tests     
F-tests     
Significance of difference of means test      
Significance of difference of variances test  
Wilcoxon tests   
Kolmogorov–Smirnov test   


## Maximum likelihood and Bayesian inference

Maximum likelihood estimation of Gaussian distribution parameters + AIC   
Maximum likelihood estimation of Weibull distribution parameters + AIC   
Implementation of binomial generalized linear model + AIC + BIC      
Bayesian posterior estimation for a mixture of betas        
Bayesian posterior estimation for binomial beta distribution           
   
   
## Machine learning     

Support Vector Machine classifier   
Naive Bayes Classifier   
Random Forest Classifier   
Linear Discriminant Analysis   
Limma linear model for differential gene expression and computation of Residual sum of squares for downstream analysis      
   
   
## Variant Analysis (GWAS)   

Carries out genome-wide association analysis using parallelized code to perform it really fast.  
Input genotype and phenotype data and get significance measures by fitting a generalized linear model per SNP.      


