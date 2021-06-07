# Program: calculate_principal_components.sh
# Author: Rene Welch 
# Created: 2020-09-01
# Updated: 2020-09-17
# Purpose: To compute top K principal components on a vcf.gz file

workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/extracted_data/SCCS/imputation/merge_impute/pca_tmp
outdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/results/2020_09_01_merge_impute_sccs/pca

plink2=/ua/rwelch/bin/plink2

vcf=$workdir/SCCS_merged_vdbp_cancer_rsq50.vcf.gz
outfile=$outdir/SCCS_merged_vdbp_cancer_rsq50

$plink2 --vcf $vcf \
  --threads 12 --memory 64000 \
  --pca 20 --out $outfile
