# Program: create_merged_index.sh
# Author: Rene Welch 
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Create indexes for merged vcf files

/z/Comp/onglab/programs/bcftools/bin/bcftools index SCCS_merged_vdbp_cancer_rsq50.vcf.gz 
/z/Comp/onglab/programs/bcftools/bin/bcftools index SCCS_merged_vdbp_cancer_rsq30.vcf.gz 