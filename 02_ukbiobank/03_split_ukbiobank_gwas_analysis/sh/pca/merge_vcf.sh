
work_dir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD

in_dir=$work_dir/raw_data/UKBiobank/2020_02/genotype/ldprune
out_dir=$work_dir/results/2021_07_01_ukbiobank_stratified_gwas/plink2

echo $in_dir
echo $out_dir

outvcf=$in_dir/ukbiobank_prune_ld.vcf.gz

bcftools=/ua/rwelch/bin/bcftools-1.9/bin/bcftools

$bcftools concat --threads 12 --output $outvcf --output-type z \
  $in_dir/chr1.vcf $in_dir/chr2.vcf $in_dir/chr3.vcf $in_dir/chr4.vcf $in_dir/chr5.vcf \
  $in_dir/chr6.vcf $in_dir/chr7.vcf $in_dir/chr8.vcf $in_dir/chr9.vcf $in_dir/chr10.vcf \
  $in_dir/chr11.vcf $in_dir/chr12.vcf $in_dir/chr13.vcf $in_dir/chr14.vcf $in_dir/chr15.vcf \
  $in_dir/chr16.vcf $in_dir/chr17.vcf $in_dir/chr18.vcf $in_dir/chr19.vcf $in_dir/chr20.vcf \
  $in_dir/chr21.vcf $in_dir/chr22.vcf
