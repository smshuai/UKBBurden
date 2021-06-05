#!/usr/bin/env python
# coding: utf-8

import sys
import hail as hl
hl.init(quiet=True, spark_conf={"spark.ui.enabled": "False"})

if __name__ == '__main__':
    inVCF, outVCF = sys.argv[1:]
    mt = hl.import_vcf(inVCF, min_partitions=4, reference_genome='GRCh38', array_elements_required=False, force_bgz=True)
    n = mt.count()
    n_var0 = n[0]
    
    # Genotype quality control
    ## Set filter condition for AB
    mt = mt.annotate_entries(AB = (mt.AD[1] / hl.sum(mt.AD)))
    filter_condition_ab = ((mt.GT.is_hom_ref() & (mt.AB <= 0.1)) | (mt.GT.is_het() & (mt.AB >= 0.2) & (mt.AB <= 0.8)) | (mt.GT.is_hom_var() & (mt.AB >= 0.9)))
    filter_condition = (mt.GQ >= 20) & (mt.DP >= 10) & filter_condition_ab
    fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition))
    print(f'>{fraction_filtered * 100:.2f}% entries filtered out of downstream analysis.')
    mt = mt.filter_entries(filter_condition)
    ## Add frac of samples passing AB threshold per row
    mt = mt.annotate_rows(AB_count = hl.agg.count_where((mt.GT.is_het() & (mt.AB >= 0.2) & (mt.AB <= 0.8))))

    # Variant quality control
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((mt.AB_count > 0) & (mt.variant_qc.call_rate >= 0.9) & (mt.variant_qc.p_value_hwe >= 1e-15))

    n = mt.count()
    n_var1 = n[0]
    print(f'variantsPASS,{outVCF},{n_var1},{n_var0}')
    
    # Output
    mt = mt.annotate_rows(info = mt.info.annotate(AF=mt.variant_qc.AF[1]))
    hl.export_vcf(mt, outVCF)
    hl.stop()