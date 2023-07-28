
plink2=/ua/rwelch/bin/plink2

out=plink2/nosnp/vitd_colorectal_nosnp
samples=plink2/gwas_samples.txt
pheno=plink2/gwas_phenotype.txt
covar=plink2/gwas_covariates.txt

$plink2 --threads 12 --out $out \
  --memory 256000 --keep $samples \
  --pheno $pheno --covar covar \
  --glm no-snp 'mperm='1000 'standard-beta'