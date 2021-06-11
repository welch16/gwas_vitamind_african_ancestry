# Program: 02a_check_minimac4_results.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-01
# Purpose: Review minimac4 results

library(magrittr)
library(tidyverse)

data_dir <- here::here("extracted_data/SCCS/imputation/imputed_by_chrom")

info_file <- list.files(data_dir, pattern = "info")

old_queue <-
  read_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor",
      "queues", "2020_09_01_minimac4_impute.csv"),
    col_names = FALSE) %>%
  rename(prefix = X3)

old_queue %>%
  filter(
    ! str_c(prefix, ".info") %in%
      file.path(data_dir, info_file)) %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor",
      "queues", "2020_09_01_minimac4_impute_rerun.csv"),
    col_names = FALSE)