# Program: 01_compute_summary_statistics.sh
# Author: Rene Welch 
# Created: 2019-05-23
# Updated: 2019-08-30
# Purpose: Compute summary statistics per SNP for the rawdata

raw_file=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/raw_data/Vit_D_MEGA_GSproject/PLINK_060718_0135/Vit_D_MEGA
out_dir=/z/Comp/uwcccseq/rwelch/WarrenAndersen_Shaneda/extracted_data

# missing samples
mkdir -p $outdir/missing

plink --file $raw_file \
  --missing \
  --out $out_dir/missing

# allele frequency
mkdir -p $outdir/allele_freq

plink --file $raw_file \
  --freq \
  --out $out_dir/allele_freq
  
# Hardy-Weinberg 
mkdir -p $outdir/hwe

plink --file $raw_file \
  --hardy \
  --out $out_dir/hardy

