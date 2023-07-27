# Program: 04_convert_snp_stats_to_text.R
# Author: Rene Welch
# Created: 2019-08-30
# Updated: 2019-09-06
# Purpose: Converts the SNP statistics to text format

library(magrittr)
library(tidyverse)

indr <- here::here("results/2019_08_30_filter_file_correction/rds")

summary_stats <- readRDS(
  file.path(indr, "2019_08_30_plink_summary_stats.rds"))

summary_stats %>%
  write_csv(file.path(indr, "2019_08_30_plink_summary_stat.csv"))
