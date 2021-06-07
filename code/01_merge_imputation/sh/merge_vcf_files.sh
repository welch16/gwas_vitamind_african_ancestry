# Program: merge_vcf_files.sh
# Author: Rene Welch 
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Merge vcf files of the imputed data, after filtering variants with rsq score greater or equal than 30 or 50

workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/extracted_data/SCCS/imputation/merge_impute/qc_filtered

bcftools=/z/Comp/onglab/programs/bcftools/bin/bcftools

$bcftools concat --threads 12 --output $workdir/SCCS_merged_vdbp_cancer_rsq30.vcf.gz --output-type z \
  $workdir/SCCS_merged_vdbp_cancer_chr1_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr2_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr3_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr4_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr5_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr6_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr7_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr8_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr9_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr10_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr11_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr12_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr13_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr14_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr15_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr16_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr17_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr18_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr19_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr20_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr21_rsq30.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr22_rsq30.vcf

$bcftools concat --threads 12 --output $workdir/SCCS_merged_vdbp_cancer_rsq50.vcf.gz --output-type z \
  $workdir/SCCS_merged_vdbp_cancer_chr1_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr2_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr3_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr4_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr5_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr6_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr7_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr8_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr9_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr10_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr11_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr12_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr13_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr14_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr15_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr16_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr17_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr18_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr19_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr20_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr21_rsq50.vcf \
  $workdir/SCCS_merged_vdbp_cancer_chr22_rsq50.vcf

