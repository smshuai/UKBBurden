#!/usr/bin/env bash

inVCF=$1
out=$2
BED=$3 # biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/xgen_plus_spikein.GRCh38.bed
REF=$4

# Norm and split variants
bcftools norm --threads 4 -m -any --check-ref w -R $BED -f $REF -Oz -o ${out}.norm.vcf.gz $inVCF
