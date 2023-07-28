
bgenfile=$1
samplefile=$2

plink2 --bgen $bgenfile ref-first --sample $samplefile --export vcf 