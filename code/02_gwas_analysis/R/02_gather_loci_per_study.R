# Program: 02_gather_loci_per_study.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-08
# Purpose: Separates SNPs into different loci based on p.value threshold and
#   genomic distance. For downstream fine-mapping analysis

library(magrittr)
library(tidyverse)
library(UWCCC.GWAS)


snp_dist <- 250e3
pval_thr <- 5e-8

sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs")) %>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat))) %>%
  nest(chrom_data = c(pos, id, ref, alt, a_1, beta, se, t_stat, p)) %>%
  mutate(
    chrom_loci = map(chrom_data, dplyr::select, id, pos, p),
    chrom_loci = map(chrom_loci,
      UWCCC.GWAS::partition_by_snp, snp_dist, pval_thr))

sccs_vdbp %>%
  select(-chrom_data) %>%
  unnest(cols = c(chrom_loci)) %>%
  qs::qsave(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    str_c("sccs_vitdbp_loci_snpdist",
      scales::comma(snp_dist, scale = 1e-3, suffix = "kp"), "_pval_thr",
      scales::scientific(pval_thr),
    ".qs")))
