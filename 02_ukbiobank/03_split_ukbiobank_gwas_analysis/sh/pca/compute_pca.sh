
work_dir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD

in_dir=$work_dir/raw_data/UKBiobank/2020_02/genotype/ldprune
out_dir=$work_dir/results/2021_07_01_ukbiobank_stratified_gwas/plink2
# temp_dir=/scratch/rwelch/pca

echo $in_dir
echo $out_dir

plink2=/ua/rwelch/bin/plink2

$plink2 --vcf $in_dir/ukbiobank_prune_ld.vcf.gz  \
 --maf 0.01 --hwe 1e-5 --geno 0.05 \
  --pca 20 approx \
  --threads 12 \
  --out $out_dir/ukbiobank_pca
