# Program: 02_filter_and_flip_genotype.sh
# Author: Rene Welch 
# Created: 2019-08-30
# Updated: 2019-08-30
# Purpose: Filter the snps according to the criteria in the grant and makes the
#   snp alleles to be expressed in the '+' strand

set -u 

echo $(hostname)

basedr=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda

indr=$basedr/raw_data/Vit_D_MEGA_GSproject/PLINK_060718_0135
outdr=$basedr/extracted_data/genotype

echo "Input directory:" ${indr}
echo "Output directory:" ${outdr}

echo "PLINK version ~~~~~~"
plink --version
## this command only keeps 1 snp
## plink --file $indr/Vit_D_MEGA --snp 1:10952040-T-A --out $outdr/Vit_D_MEGA --make-bed

## saves the snps selected according to the criteria in the grant

echo "Filtering SNPs according to grant criteria ~~~~~~~"
plink --file $indr/Vit_D_MEGA --extract $outdr/2019_08_30_snps_to_subset.txt \
    --out $outdr/Vit_D_MEGA_grant --make-bed --snps-only

echo "Flipping bases for SNPs where the minor allele is in neg. strand"
plink --bfile $outdr/Vit_D_MEGA_grant --recode --out $outdr/Vit_D_MEGA_grant
plink --file $outdr/Vit_D_MEGA_grant --flip $outdr/2019_08_30_snps_to_flip.txt \
    --recode --out $outdr/Vit_D_MEGA_grant_wflipped_2019_08_30 --make-bed --snps-only
