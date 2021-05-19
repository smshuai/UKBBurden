#!/usr/bin/env bash

## QC pVCF files.
## Need VCFtools (0.1.16); picard

inVCF=$1
out=$2
BED=$3 # biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/xgen_plus_spikein.GRCh38.bed
REF=$4

# Norm and split variants
bcftools norm --threads 4 -m -any --check-ref w -R $BED -f $REF -Oz -o ${out}.norm.vcf.gz 

# Filter by depth, quality, variant missingness and targeted regions
vcftools --minDP 10 --minQ 20 --max-missing 0.9 --gzvcf ${inVCF} \
	--bed $BED --recode --recode-INFO-all --out ${out}.step1 --stdout | \
	bgzip -c > ${out}.step1.vcf.gz

# Filter by allele balance
picard FilterVcf -I ${out}.step1.vcf.gz -O ${out}.step2.vcf.gz --MIN_AB 0.2

vcftools --gzvcf ${inVCF} --out ${out} --missing-indv

cut=0.000000000000001 # 1e-15
vcftools --hardy --gzvcf ${inVCF} --out ${out}
