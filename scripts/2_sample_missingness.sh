#!/usr/bin/env bash

## QC pVCF files.
## Need VCFtools (0.1.16); picard

inVCF=$1
out=$2

vcftools --gzvcf ${inVCF} --out ${out} --missing-indv
