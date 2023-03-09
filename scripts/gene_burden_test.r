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
    make_option(c("-s", "--spa"), type="logical", default=FALSE, 
              help="use SPA test for binary traits [default = %default]", metavar="TRUE/FALSE"),
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
    if (n_ctrl/n_case > 100) {
        set.seed(1993)
        max_ctrl <- n_case * 100
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

# Binary - SPA test
if (dtype == 'binary' & opt$spa) {
    library(SPAtest)
    # split geno into chunks for parallel computing
    x <- 1:nrow(geno_keep)
    chunks <- split(x, cut(seq_along(x), opt$nfork, labels = FALSE))
    res <- foreach(ix=1:opt$nfork, .combine='rbind') %dopar% {
        tmp <- SPAtest::ScoreTest_SPA(genos = geno_keep[chunks[[ix]],], pheno = pheno$y,
                                      cov = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]),
                                      beta.out = TRUE, beta.Cutoff = 1e-5)
        tmp <- as.data.table(tmp)
    }
}

# Binary - logistic regression
if (dtype == 'binary' & !opt$spa) {
    # loop over genes
    res <- foreach(g=geno$Gene, .combine='rbind') %dopar% {
        pheno[, x:=geno_keep[g, ]]

        fit <- speedglm.wfit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]), family = binomial('logit'))

        coef_tb <- coef(summary(fit))
        if ('x' %in% rownames(coef_tb)) {
            ret <- c(g, coef_tb['x', ])
        } else {
            ret <- c(g, rep(NA, 4))
        }
        ret
    }
}

# Continous traits - linear regression
if (dtype == 'quant') {
    # loop over genes
    res <- foreach(g=geno$Gene, .combine='rbind') %dopar% {
        pheno[, x:=geno_keep[g, ]]

        fit <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]))

        coef_tb <- coef(summary(fit))
        if ('x' %in% rownames(coef_tb)) {
            ret <- c(g, coef_tb['x', ])
        } else {
            ret <- c(g, rep(NA, 4))
        }

        ret
    }
}

parallel::stopCluster(cl)


# Prepare output tables
res <- setDT(as.data.frame(res))
if (dtype == 'binary' & opt$spa) {
    colnames(res) <- c('P', 'P.NA', 'is.converge', 'beta', 'SEbeta')
    res <- cbind(Gene=rownames(geno_keep), res)
}
if (dtype == 'binary' & !opt$spa) colnames(res) <- c('Gene', 'EST', 'SE', 'T', 'P')
if (dtype == 'quant') colnames(res)  <- c('Gene', 'EST', 'SE', 'Z', 'P')

# add allele counts
if (dtype == 'binary') {
    control_hom_ref_count <- rowSums(geno_keep[,pheno$y==0]==0)
    control_het_alt_count <- rowSums(geno_keep[,pheno$y==0]==1)
    control_hom_alt_count <- rowSums(geno_keep[,pheno$y==0]==2)
    case_hom_ref_count <- rowSums(geno_keep[,pheno$y==1]==0)
    case_het_alt_count <- rowSums(geno_keep[,pheno$y==1]==1)
    case_hom_alt_count <- rowSums(geno_keep[,pheno$y==1]==2)
    ac <- data.table(case_hom_ref_count, case_het_alt_count, case_hom_alt_count, control_hom_ref_count, control_het_alt_count, control_hom_alt_count)
    colnames(ac) <- c('Case_HOM_REF', 'Case_HET_ALT', 'Case_HOM_ALT', 'Control_HOM_REF', 'Control_HET_ALT', 'Control_HOM_ALT')
}

if (dtype == 'quant') {
    hom_ref_count <- rowSums(geno_keep==0)
    het_alt_count <- rowSums(geno_keep==1)
    hom_alt_count <- rowSums(geno_keep==2)
    ac <- data.table(hom_ref_count, het_alt_count, hom_alt_count)
    colnames(ac) <- c('HOM_REF', 'HET_ALT', 'HOM_ALT')
}


res <- cbind(res, ac)
fwrite(res, out, sep='\t')