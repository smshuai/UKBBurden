#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)

x <- fread(args[1], header=TRUE)

# There are 86 metadata columns in the file from VEP mask
# 87:ncol(x) are samples
# 0 - Hom_Ref; 1 - Het_Alt; 2 - Hom_Alt; 3 - Missing
# Missing to Hom_Ref
for (j in 87:ncol(x))
    set(x, which(x[[j]]==3), j, 0)

# Gene level (max of all variants)
xmat = x[, lapply(.SD, max), by=Gene, .SDcols=87:ncol(x)]

fwrite(xmat, args[2], sep='\t')