# Program: identity_by_descent.sh
# Author: Rene Welch
# Created: 2020-09-15
# Updated: 2021-05-03
# Purpose: Performs identidy by descent calculations

workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/extracted_data/SCCS/imputation/merge_impute/pca_tmp
vcf_file=SCCS_merged_vdbp_cancer_chr9_rsq50.vcf

PLINK=/s/pkg/linux64/plink/1.90/plink

$PLINK --vcf $workdir/$vcf_file \
  --genome full \
  --double-id \
  --out /z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/results/2021_03_19_identity_by_descent/plink/identity_by_descent