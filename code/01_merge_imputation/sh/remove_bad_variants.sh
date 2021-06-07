# Program: remove_bad_variants.sh
# Author: Rene Welch 
# Created: 2020-09-01
# Updated: 2020-09-11
# Purpose: To remove variants from a vcf file with not-identified alleles like '+' / '-' because this were causing problems when trying to impute the data

file=$1
tmp="$file"_tmp

echo $file

cat $file | grep -v "+" > $tmp
rm -f $file
mv $tmp $file
rm -fr $tmp