# -*- UTF-8 -*-
# !/usr/bin/env R4.2
# by zhiyuan

library(MatrixEQTL)
library(getopt)

spec <- matrix(c(
    "genelocfile", "geloc", 1, "character",
    "genofile", "geno", 1, "character",
    "expressionfile", "exp", 1, "character",
    "snpposinfo", "info", 1, "character",
    "covariate", "cov", 2, "character",
    "outpath", "out", 2, "character"
), byrow = TRUE, ncol = 4)
spec
opt <- getopt(spec)
opt

# define usage function
print_usage <- function(spec = NULL) {
    cat(getopt(spec, usage = TRUE))
    cat("Usage example: \n")
    cat("Usage example:
      1 Rscript xxx.R --genofile genotype.txt --phenofile phenotype.txt -- snpposinfo snpposinfo.txt --outpath /share/nas1/map
      2 Rscript GOClassificationMap.R --infile in.dat --rotation 0.5 --outpath /share/nas1 --method log

      Options:
      --genelocfile		-geloc	character	[forced]
      --genofile	-geno   character	[forced]
      --expressionfile	-exp  character	[forced]
      --snpposinfo 	-info	character	[forced]
      --covariate       -cov    character       [optional]
      --outpath		-out	character	[optional]
      \n")
    q(status = 1)
}

# Model selection
useModel <- modelLINEAR

# Genotype data
SNP_file_name <- opt$genofile
snps_location_file_name <- opt$snpposinfo
# Gene expression file name
expression_file_name <- opt$expressionfile
gene_location_file_name <- opt$genelocfile
# Covariates file name
# Set to character() for no covariates
covariates_file_name <- opt$covariate
# Output file name
output_file_name_cis <- paste(opt$outpath, "cis_qqnorm_result.txt", sep = "/")
output_file_name_tra <- paste(opt$outpath, "trans_qqnorm_result.txt", sep = "/")
# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 0.9
pvOutputThreshold_tra <- 0.9
# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric()
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist <- 100000
## Load genotype data

snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
# the TAB character
snps$fileOmitCharacters <- "NA"
# denote missing values;
snps$fileSkipRows <- 1
# one row of column labels
snps$fileSkipColumns <- 1
# one column of row labels
snps$fileSliceSize <- 20000
# read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)
## Load gene expression data

gene <- SlicedData$new()
gene$fileDelimiter <- "\t"
# the TAB character
gene$fileOmitCharacters <- "NA"
# denote missing values;
gene$fileSkipRows <- 1
# one row of column labels
gene$fileSkipColumns <- 1
# one column of row labels
gene$fileSliceSize <- 2000
# read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)
## Load covariates

cvrt <- SlicedData$new()
cvrt$fileDelimiter <- "\t"
# the TAB character
cvrt$fileOmitCharacters <- "NA"
# denote missing values;
cvrt$fileSkipRows <- 1
# one row of column labels
cvrt$fileSkipColumns <- 1
# one column of row labels
if (length(covariates_file_name) > 0) {
    cvrt$LoadFile(covariates_file_name)
}

## Run the analysis

snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)
