#!/usr/bin/env Rscript
# Test the combined effects of SNV and CNV

# Example:
# Rscript gene_burden_test_combined.r -c COVAR -g1 CNV -g2 SNV -p PHENO -d quant -n 8
# Rscript gene_burden_test_combined.r -c COVAR -g1 CNV -g2 SNV -p PHENO -d binary -s -n 8

library(data.table)
library(foreach)
library(parallel)
library(speedglm)
library(lmtest)
library(optparse)
 
option_list = list(
    make_option(c("-c", "--covar"), type="character", default=NULL, 
              help="covariate file path", metavar="path"),
    make_option(c("--geno1"), type="character", default=NULL, 
              help="genotype file path for CNV", metavar="path"),
    make_option(c("--geno2"), type="character", default=NULL, 
              help="genotype file path for SNV", metavar="path"),
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
covar <- fread(opt$covar, header=T)
geno1 <- fread(opt$geno1, header=T)
geno2 <- fread(opt$geno2, header=T)
pheno <- fread(opt$pheno, header=T, col.names = c('eid', 'y'))
dtype <- opt$dtype
out <- opt$out

pheno <- merge(pheno, covar, by='eid')
pheno[, eid:=as.character(eid)]
setkey(pheno, eid)
keep <- Reduce(intersect, list(pheno$eid, colnames(geno1), colnames(geno2)))
pheno <- pheno[keep, ]
cat("Number of samples:", length(keep), '\n') # No. of samples used

# Only use shared genes for combined tests
setnames(geno1, 1, 'Gene') # The first one must be 'Gene'
setnames(geno2, 1, 'Gene') # The first one must be 'Gene'
genes1 <- geno1$Gene
genes2 <- geno2$Gene
geno1 <- as.matrix(geno1[,-1], rownames.value=genes1)[,keep]
geno2 <- as.matrix(geno2[,-1], rownames.value=genes2)[,keep]
mac1 <- rowSums(geno1)
mac2 <- rowSums(geno2)
# At least 3 in SNV and CNV 
keep_genes <- intersect(genes1[mac1 >= 3], genes2[mac2 >= 3])
geno1 <- geno1[keep_genes,]
geno2 <- geno2[keep_genes,]
geno_joint <- ifelse(geno1 | geno2, 1, 0)
cat("Number of genes:", length(keep_genes), '\n')


## Subsampling for binary to 1:100 case:control
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

# Run test in parallel
cl <- parallel::makeForkCluster(opt$nfork)
doParallel::registerDoParallel(cl)

# Binary - SPA test
if (dtype == 'binary' & opt$spa) {
    require(SPAtest)
    # split geno into chunks for parallel computing
    x <- 1:nrow(geno_joint)
    chunks <- split(x, cut(seq_along(x), opt$nfork, labels = FALSE))
    res <- foreach(ix=1:opt$nfork, .combine='rbind') %dopar% {
        tmp <- SPAtest::ScoreTest_SPA(genos = geno_joint[chunks[[ix]],], pheno = pheno$y,
                                      cov = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]),
                                      beta.out = TRUE, beta.Cutoff = 1e-5)
        tmp <- as.data.table(tmp)
    }
}

# Binary - logistic regression
if (dtype == 'binary' & !opt$spa) {
    # loop over genes
    res <- foreach(g=keep_genes, .combine='rbind') %dopar% {
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
    res <- foreach(g=keep_genes, .combine='rbind') %dopar% {
        if (TRUE) {
            pheno[, `:=`(x1=geno1[g, ], x2=geno2[g,])]
            # Fit different nested models
            fit_full <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid')]))
            fit_x1 <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid', 'x2')]))
            fit_x2 <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid', 'x1')]))
            fit_null <- speedlm.fit(y = pheno$y, X = as.matrix(pheno[, .SD, .SDcols=-c('y', 'eid', 'x1', 'x2')]))
            # LRT test
            lrt.x2 <- lrtest(fit_full, fit_x1)[2,c('Chisq', 'Pr(>Chisq)')] # test for x2
            lrt.x1 <- lrtest(fit_full, fit_x2)[2,c('Chisq', 'Pr(>Chisq)')] # test for x1
            lrt.x1x2 <- lrtest(fit_full, fit_null)[2,c('Chisq', 'Pr(>Chisq)')] # test for x1 + x2
            coef_tb <- coef(summary(fit_full))
            ret <- c(g, lrt.x1x2, lrt.x1, lrt.x2)
            if ('x1' %in% rownames(coef_tb)) {
                ret <- c(ret, coef_tb['x1', ])
            } else {
                ret <- c(ret, rep(NA, 4))
            }
            if ('x2' %in% rownames(coef_tb)) {
                ret <- c(ret, coef_tb['x2', ])
            } else {
                ret <- c(ret, rep(NA, 4))
            }
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
if (dtype == 'quant') colnames(res)  <- c('Gene', 'Chisq.joint', 'Chisq.p.joint',
                                       'Chisq.1', 'Chisq.p.1',
                                       'Chisq.2', 'Chisq.p.2',
                                       'EST.1', 'SE.1', 'Z.1', 'P.1',
                                       'EST.2', 'SE.2', 'Z.2', 'P.2')

# add allele counts
if (dtype == 'binary') {
    control_cnv_ref_count <- rowSums(geno1[,pheno$y==0]==0)
    control_cnv_alt_count <- rowSums(geno1[,pheno$y==0]>0)
    control_snv_ref_count <- rowSums(geno2[,pheno$y==0]==0)
    control_snv_alt_count <- rowSums(geno2[,pheno$y==0]>0)
    control_joint_ref_count <- rowSums(geno_joint[,pheno$y==0]==0)
    control_joint_alt_count <- rowSums(geno_joint[,pheno$y==0]>0)

    case_cnv_ref_count <- rowSums(geno1[,pheno$y==1]==0)
    case_cnv_alt_count <- rowSums(geno1[,pheno$y==1]>0)
    case_snv_ref_count <- rowSums(geno2[,pheno$y==1]==0)
    case_snv_alt_count <- rowSums(geno2[,pheno$y==1]>0)
    case_joint_ref_count <- rowSums(geno_joint[,pheno$y==1]==0)
    case_joint_alt_count <- rowSums(geno_joint[,pheno$y==1]>0)

    ac <- data.table(case_cnv_ref_count, case_cnv_alt_count,
                     case_snv_ref_count, case_snv_alt_count,
                     case_joint_ref_count, case_joint_alt_count,
                     control_cnv_ref_count, control_cnv_alt_count,
                     control_snv_ref_count, control_snv_alt_count,
                     control_joint_ref_count, control_joint_alt_count)
    colnames(ac) <- c('Case_CNV_REF', 'Case_CNV_ALT',
                      'Case_SNV_REF', 'Case_SNV_ALT',
                      'Case_JOINT_REF', 'Case_JOINT_ALT',
                      'Control_CNV_REF', 'Control_CNV_ALT',
                      'Control_SNV_REF', 'Control_SNV_ALT',
                      'Control_JOINT_REF', 'Control_JOINT_ALT')
}

if (dtype == 'quant') {
    cnv_ref_count <- rowSums(geno1==0)
    snv_ref_count <- rowSums(geno2==0)
    cnv_alt_count <- rowSums(geno1>0)
    snv_alt_count <- rowSums(geno2>0)
    joint_ref_count <- rowSums(geno_joint==0)
    joint_alt_count <- rowSums(geno_joint==1)
    ac <- data.table(cnv_ref_count, cnv_alt_count, snv_ref_count,
                     snv_alt_count, joint_ref_count, joint_alt_count)
    colnames(ac) <- c('CNV_REF', 'CNV_ALT', 'SNV_REF',
                      'SNV_ALT', 'JOINT_REF', 'JOINT_ALT')
}


res <- cbind(res, ac)
fwrite(res, out, sep='\t')