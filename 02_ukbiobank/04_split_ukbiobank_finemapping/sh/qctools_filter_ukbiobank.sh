
workdir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/GWAS_colorecal_vitD
QCTOOL="$workdir"/tools/QCTOOL

bgenfile=$1
snpids=$2
outfile=$3

samplefile="$workdir"/raw_data/UKBiobank/2020_02/genotype/ukb35256_imp_v3.sample

$QCTOOL -g $bgenfile \
  -incl-samples "$workdir"/results/2022_08_04_ukbiobank_finemapping/vcf/afrocaribbean_samplefile.txt \
  -incl-snpids $snpids \
  -og $outfile

