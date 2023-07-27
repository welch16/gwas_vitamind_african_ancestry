# Program: 01_build_filtering_list.R
# Author: Rene Welch
# Created: 2019-08-30
# Updated: 2019-08-30
# Purpose: Builds a filtering list of SNPs based on the grant criteria

library(magrittr)
library(tidyverse)

indr <- here::here("extracted_data")
dr <- here::here("results/2019_04_22_plink_files_review/rds")
outdr <- here::here("results/2019_08_30_filter_file_correction/rds")

### load minor allele frequency
freq <- readr::read_delim(file.path(indr, "allele_freq/Vit_D_MEGA.frq"),
  delim = " ", trim_ws = TRUE)
freq %<>% set_names(str_to_lower(names(freq)))
freq %<>% mutate(chr = as.character(chr))

### load hardy-weinberg equilibrium
hw <- readr::read_delim(
    file.path(indr, "hardy_weinberg/Vit_D_MEGA.hwe"),
      delim = " ", trim_ws = TRUE)
hw %<>% set_names(str_to_lower(names(hw))) %>% select(-x10)
hw %<>% mutate(chr = as.character(chr))


#### load missing variants
miss_v <- readr::read_delim(
    file.path(indr,"missing/Vit_D_MEGA.lmiss"), delim = " ", trim_ws = TRUE)
miss_v %<>% set_names(str_to_lower(names(miss_v)))
miss_v %<>% mutate(chr = as.character(chr))

plink_stats  <-
  inner_join(freq,
  select(hw, -chr, -a1, -a2), by = "snp") %>%
  inner_join(select(miss_v, -chr), by = "snp")

plink_stats %>%
  saveRDS(file.path(outdr, "2019_08_30_plink_summary_stats.rds"))

### Parameters in grant
maf_thr <- 0.05
f_miss_thr <- 0.95
pval_thr <- 1e-6

plink_stats %<>%
  filter(maf >= maf_thr & (1 - f_miss)  >= f_miss_thr & p > pval_thr)
plink_stats %>%
  saveRDS(file.path(outdr, "2019_08_30_plink_summary_stats_filtered.rds"))
