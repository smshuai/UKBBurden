#!/usr/bin/env bash

## QC pVCF files.
## Need VCFtools (0.1.16); picard

inVCF=$1
out=$2

py_path=/nfs/research/birney/users/shimin/ukbb/wes/200k/UKBBurden/scripts
export PYSPARK_SUBMIT_ARGS='--driver-memory 60g --executor-memory 60g pyspark-shell'
python $py_path/3_hail_qc.py ${inVCF} ${out}
