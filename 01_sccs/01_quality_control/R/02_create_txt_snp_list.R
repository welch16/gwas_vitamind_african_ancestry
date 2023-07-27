# Program: 02_create_txt_snp_list.R
# Author: Rene Welch
# Created: 2019-08-30
# Updated: 2019-08-30
# Purpose: Creates txt list of the snps to filter

library(magrittr)
library(tidyverse)

dr <- here::here("results/2019_08_30_filter_file_correction/rds")
outdr <- here::here("extracted_data", "genotype")

plink <- readRDS(file.path(dr, "2019_08_30_plink_summary_stats_filtered.rds"))
plink %>%
  select(snp) %>%
  write_delim(
    file.path(outdr, "2019_08_30_snps_to_subset.txt"), delim = "\t",
    col_names = FALSE)

flip <- readRDS(
  file.path(dr, "2019_08_30_filtered_variants_to_flip_w_summary_stats.rds"))
flip %>%
  select(snp.name) %>%
  write_delim(file.path(outdr, "2019_08_30_snps_to_flip.txt"),
    delim = "\t", col_names = FALSE)
