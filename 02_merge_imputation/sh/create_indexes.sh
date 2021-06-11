# Program: create_indexes.sh
# Author: Rene Welch 
# Created: 2020-09-01
# Updated: 2020-09-06
# Purpose: To index chromosome specific vcf files

bcftools=/z/Comp/onglab/programs/bcftools/bin/bcftools
# workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD/extracted_data/SCCS/imputation/aux_impute
workdir=$1

for file in $workdir/*gz
do
  $bcftools index $file
done