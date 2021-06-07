# Program: 05b_subset_vcf_files_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Writes queue to subset vcf files

library(magrittr)
library(tidyverse)

workdir <- here::here("extracted_data", "SCCS", "imputation", "merge_impute")
vcf_files <- list.files(workdir, pattern = "vcf")

vcf_files %<>%
  str_subset("csi", negate = TRUE)

queue <- tibble(vcf = vcf_files) %>%
  mutate(
    chr = str_split(vcf, "[\\_|\\.]"),
    chr = map_chr(chr, ~ str_subset(., "chr")))

queue_dir <- here::here("results", "2020_09_01_merge_impute_sccs", "variants")

snps <- tibble(
  rsq30 = list.files(queue_dir, pattern = "rsq30"),
  rsq50 = list.files(queue_dir, pattern = "rsq50")) %>%
  mutate(
    chr = str_split(rsq30, "\\_"),
    chr = map_chr(chr, 1))

queue %<>%
  inner_join(snps, by = "chr") %>%
  pivot_longer(starts_with("rsq"), names_to = "rsq", values_to = "snps") %>%
  mutate(
    out = str_replace(vcf, ".vcf.gz", str_c("_", rsq)),
    vcf = file.path(workdir, vcf),
    snps = file.path(queue_dir, snps),
    out = file.path(workdir, "qc_filtered", out)) %>%
  select(-chr, -rsq)

queue %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor", "queues",
      "2020_09_15_subset_snps.csv"), col_names = FALSE)