#!/usr/bin/env bash

inVCF=$1
outVCF=$2
SIF="/nfs/research/birney/users/shimin/VEP/ensembl-vep_release_104.0.sif"
VEPDATA="/nfs/research/birney/users/shimin/VEP/vep_data/"

vep -i $inVCF -o $outVCF \
    --plugin dbNSFP,$VEPDATA/extdata/dbNSFP4.2a_grch38.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
    --plugin LoF,loftee_path:$VEPDATA/loftee,human_ancestor_fa:$VEPDATA/extdata/human_ancestor.fa.gz,gerp_bigwig:$VEPDATA/extdata/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:$VEPDATA/extdata/loftee.sql \
    --dir_plugins $VEPDATA/Plugins/ --dir_cache $VEPDATA --cache \
    --everything --offline --buffer_size 2000 --compress_output bgzip --vcf


# singularity exec -B $VEPDATA:/opt/vep/.vep/ $SIF \
#    vep -i $inVCF -o $outVCF \
#    --plugin dbNSFP,/opt/vep/.vep/extdata/dbNSFP4.2a_grch38.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
#    --plugin LoF,loftee_path:/opt/vep/.vep/loftee,human_ancestor_fa:/opt/vep/.vep/extdata/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/extdata/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/extdata/loftee.sql \
#    --dir_plugins /opt/vep/.vep/Plugins/ --dir_cache /opt/vep/.vep/ --cache \
#    --everything --offline --buffer_size 5000 --compress_output bgzip --vcf --fork 2
    
    
