#!/usr/bin/env Rscript
# Example: gene_level_burden_test.r /nfs/research/birney/projects/association/gene_burden/genotypes/mask1_pLoF_canonical_rare_v1/ukb23156.mask1.c22.gene_burden.tsv /nfs/research/birney/users/tomas/gene_modules/cleaned_and_linked_phenotypes/cancer_ukb_f_id_40009_EFO_0000311-8.txt quant ukb23156.mask1.c22.cancer_ukb_f_id_40009_EFO_0000311-8.tsv
args = commandArgs(trailingOnly=TRUE)

library(data.table)

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
covar <- fread('/nfs/research/birney/projects/association/gene_burden/covariates//ukbb_200k_covar_v1.csv')
# geno <- fread('/nfs/research/birney/projects/association/gene_burden/genotypes/mask1_pLoF_canonical_rare_v1/ukb23156.mask1.c22.gene_burden.tsv', header=T)
geno <-  fread(args[1], header=T)
# pheno <- fread('/nfs/research/birney/users/tomas/gene_modules/cleaned_and_linked_phenotypes/cancer_ukb_f_id_40009_EFO_0000311-8.txt', header=T, col.names = c('eid', 'y'))
pheno <- fread(args[2], header=T, col.names = c('eid', 'y'))
# dtype <- 'quant' or 'binary'
dtype <- args[3]
out <- args[4]

pheno <- merge(pheno, covar, by='eid')
pheno[, eid:=as.character(eid)]
setkey(pheno, eid)

keep <- intersect(pheno$eid, colnames(geno))
print(length(keep)) # No. of samples used

pheno <- pheno[keep, ]
## Normalize y for quant
if (dtype == 'quant') pheno[, y:=RankNorm(y)]


geno_keep = as.matrix(geno[, -1])[,keep]
rownames(geno_keep) <- geno$Gene

res <- matrix(NA, nrow=nrow(geno), ncol=7, dimnames = list(rownames(geno_keep), c('est', 'stderr', 'tval', 'P', 'HOM_REF', 'HET_ALT', 'HOM_ALT')))
# loop over genes
for (g in geno$Gene) {
    pheno[, x:=geno_keep[g, ]]
    hom_ref_count = pheno[x==0, .N]
    het_alt_count = pheno[x==1, .N]
    hom_alt_count = pheno[x==2, .N]
    if (dtype == 'quant') fit <- lm(data = pheno[, -1], formula = y ~ .)
    if (dtype == 'binary') fit <- glm(data = pheno[, -1], formula = y ~ ., family = 'binomial')
    coef_tb <- coef(summary(fit))
    if ('x' %in% rownames(coef_tb)) res[g,] <- c(coef_tb['x', ], hom_ref_count, het_alt_count, hom_alt_count)
}

fwrite(as.data.table(res, keep.rownames=T), out, sep='\t')



