#!/usr/bin/env Rscript
# Example: Rscript ./UKBBurden/scripts/gene_burden_test.r -c ./covar/version2/ukbb_covar_v2_both.csv -g ./gene_burden/cnv_burden_v1/rare_del_burden/rare_del_burden_v1_c15.tsv -p ./phenotypes/quant-blood_counts-c100081/quant-lymphocyte_count_f30120-191546.txt -d quant -n 8

library(data.table)
library(foreach)
library(parallel)
library(speedglm)
library("optparse")
 
option_list = list(
    make_option(c("-c", "--covar"), type="character", default=NULL, 
              help="covariate file path", metavar="path"),
    make_option(c("-g", "--geno"), type="character", default=NULL, 
              help="genotype file path", metavar="path"),
    make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="phenotype file path", metavar="path"),
    make_option(c("-d", "--dtype"), type="character", default=NULL, 
              help="data type (quant or binary)", metavar="char"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default = %default]", metavar="char"),
    make_option(c("-n", "--nfork"), type="integer", default=1, 
              help="number of forks for parallel [default = %default]", metavar='int')
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## Helper function
# RankNorm() from package RNOmni
RankNorm <- function (u, k = 0.375)
{
  if (!is.vector(u)) {
    stop("A numeric vector is expected for u.")
  }
  if ((k < 0) || (k > 0.5)) {
    stop("Select the offset within the interval (0,0.5).")
  }
  if (sum(is.na(u)) > 0) {
    stop("Please exclude observations with missing measurements.")
  }
  n <- length(u)
  r <- rank(u)
  out <- qnorm((r - k)/(n - 2 * k + 1))
  return(out)
}

## Input
# covar <- fread('./covar/version2/ukbb_covar_v2_both.csv')
# geno <- fread('./gene_burden/cnv_burden_v1/rare_del_burden/rare_del_burden_v1_c15.tsv', header=T)
# pheno <- fread('./phenotypes/quant-blood_counts-c100081/quant-lymphocyte_count_f30120-191546.txt', header=T, col.names = c('eid', 'y'))
# dtype <- 'quant'

covar <- fread(opt$covar, header=T)
geno <-  fread(opt$geno, header=T)
pheno <- fread(opt$pheno, header=T, col.names = c('eid', 'y'))
dtype <- opt$dtype
out <- opt$out

setnames(geno, 1, 'Gene') # The first one must be 'Gene'

pheno <- merge(pheno, covar, by='eid')
pheno[, eid:=as.character(eid)]
setkey(pheno, eid)

keep <- intersect(pheno$eid, colnames(geno))
pheno <- pheno[keep, ]

## Subsampling for binary to 1:99 case:control
if (dtype == 'binary') {
    n_case <- pheno[y==1, .N]
    n_ctrl <- pheno[y==0, .N]
    if (n_ctrl/n_case > 99) {
        set.seed(1993)
        max_ctrl <- n_case * 99
        keep <- c(pheno[y==1, eid], pheno[y==0, sample(eid, max_ctrl)])
        pheno <- pheno[keep, ]
    }
}

## Normalize y for quant
if (dtype == 'quant') pheno[, y:=RankNorm(y)]


cat("Number of samples:", length(keep), '\n') # No. of samples used

geno_keep = as.matrix(geno[, -1])[,keep]
rownames(geno_keep) <- geno$Gene

# Run test in parallel
cl <- parallel::makeForkCluster(opt$nfork)
doParallel::registerDoParallel(cl)

# loop over genes
res <- foreach(g=geno$Gene, .combine='rbind') %dopar% {
    pheno[, x:=geno_keep[g, ]]
    hom_ref_count = pheno[x==0, .N]
    het_alt_count = pheno[x==1, .N]
    hom_alt_count = pheno[x==2, .N]
    
    if (dtype == 'quant') {
        fit <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]))
    }
    
    if (dtype == 'binary') {
        fit <- glm(data = pheno[, -1], formula = y ~ ., family = 'binomial')
    }
    
    coef_tb <- coef(summary(fit))
    if ('x' %in% rownames(coef_tb)) {
        ret <- c(g, coef_tb['x', ], hom_ref_count, het_alt_count, hom_alt_count)
    } else {
        ret <- c(g, rep(NA, 4), hom_ref_count, het_alt_count, hom_alt_count)
    }
    
    ret
}

res = setDT(as.data.frame(res))
if (dtype == 'binary') colnames(res) <- c('Gene', 'EST', 'SE', 'T', 'P', 'HOM_REF', 'HET_ALT', 'HOM_ALT')
if (dtype == 'quant') colnames(res)  <- c('Gene', 'EST', 'SE', 'Z', 'P', 'HOM_REF', 'HET_ALT', 'HOM_ALT')

fwrite(res, out, sep='\t')

parallel::stopCluster(cl)