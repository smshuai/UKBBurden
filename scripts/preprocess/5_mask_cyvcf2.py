#!/usr/bin/env python
# coding: utf-8

# Conda activate cyvcf (4 threads, 10 GB mem)
# Usage: 5_mask_cyvcf2.py FE_plink.chr21.recode.vep.vcf.gz chr21

import sys
from cyvcf2 import VCF
import time
import logging
logging.basicConfig(format='%(message)s')

logging.warning('Starting. Args={}'.format(sys.argv))


vcf_reader = VCF(sys.argv[1], gts012=True, threads=4)
n_record = 0
n_mask1 = 0
n_mask2 = 0

info_format = vcf_reader.get_header_type('CSQ')['Description'].split('Format: ')[1].split('|')
info_format[-1] = 'LoF_info'
# print(info_format)
max_af = info_format.index('MAX_AF')
lof = info_format.index('LoF')
canonical = info_format.index('CANONICAL')
lrt_pred = info_format.index('LRT_pred')
consequence = info_format.index('Consequence')
impact = info_format.index('IMPACT')

def apply_mask(csq_records):
    global max_af, lof, canonical, impact, consequence
    for csq in csq_records:
        csq_info = csq.split('|')
        maf = 0 if csq_info[max_af] == '' else float(csq_info[max_af])
        # Rare, canonical
        if csq_info[canonical] == 'YES' and maf < 0.01:
            if csq_info[lof] == 'HC':
                return('m1', csq_info)
            elif csq_info[consequence] == 'missense_variant' and (csq_info[impact] in ('HIGH', 'MODERATE')):
                return('m2', csq_info)
    return(None, csq_info)


# Generate genotype matrix
fgt_m1 = open(sys.argv[2] + '.mask1.genotype.tsv', 'w')
fgt_m2 = open(sys.argv[2] + '.mask2.genotype.tsv', 'w')
# Write header of sample names.
fgt_m1.write("\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT'] + info_format + vcf_reader.samples) + '\n')
fgt_m2.write("\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT'] + info_format + vcf_reader.samples) + '\n')

    
# Input VCF with VEP annotated variants
t = time.process_time()
for record in vcf_reader:
    n_record += 1
    if n_record % 1000 == 0:
        delta_t = time.process_time() - t
        logging.warning('Processing {} records, {} secs passed.'.format(n_record, delta_t))
    if record.aaf > 0.01:
        # Not rare in the current cohort, skip
        continue
    mask, info = apply_mask(record.INFO['CSQ'].split(','))
    if mask == 'm1':
        n_mask1 += 1
        # The type of genotype. hom_ref = 0 het = 1 hom_alt = 2 uncalled = None
        fgt_m1.write("\t".join([record.CHROM, str(record.POS), record.ID, record.REF, record.ALT[0]] + info + record.gt_types.astype(str).tolist()) + '\n')
    elif mask == 'm2':
        n_mask2 += 1
        fgt_m2.write("\t".join([record.CHROM, str(record.POS), record.ID, record.REF, record.ALT[0]] + info + record.gt_types.astype(str).tolist()) + '\n')

fgt_m1.close()
fgt_m2.close()

logging.warning(f'{sys.argv[2]},stats,n_mask1,n_mask2,n_record')
logging.warning(f'{sys.argv[2]},stats,{n_mask1},{n_mask2},{n_record}')

