#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)

chrom <- args[1] # 22
mask <- args[2] # 1
file_dir <- args[3]
out <- args[4]

all.files <- list.files(file_dir, pattern=sprintf("c%s_.*mask%s.gene_burden.tsv", chrom, mask), full.names=T)


l <- lapply(all.files, fread, sep="\t", header=T)
dt <- rbindlist( l )
setkey(dt, Gene)

# Merge duplicated genes
dup_genes <- dt$Gene[duplicated(dt$Gene)]
if (length(dup_genes) > 0) {
    print(sprintf('Find %d dup genes!', length(dup_genes)))
    tmp = dt[dup_genes, lapply(.SD, max), by=Gene]
    dt = rbind(tmp, dt[!dup_genes])
}

fwrite(dt, out, sep='\t')

print('Done!')