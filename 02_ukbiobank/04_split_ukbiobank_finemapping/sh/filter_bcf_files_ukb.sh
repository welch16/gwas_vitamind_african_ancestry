
workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD

bash sh/qctools_filter_ukbiobank.sh \
  "$workdir"/raw_data/UKBiobank/2020_02/gwas/ukb_imp_chr2_colosubset_gwas.bgen \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr4_afrocarib.txt \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr4_afrocarib.bgen 

bash sh/qctools_filter_ukbiobank.sh \
  "$workdir"/raw_data/UKBiobank/2020_02/gwas/ukb_imp_chr2_colosubset_gwas.bgen \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr2_afrocarib.txt \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr2_afrocarib.bgen 

bash sh/qctools_filter_ukbiobank.sh \
  "$workdir"/raw_data/UKBiobank/2020_02/gwas/ukb_imp_chr2_colosubset_gwas.bgen \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr3_afrocarib.txt \
  "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/variants_chr3_afrocarib.bgen  

