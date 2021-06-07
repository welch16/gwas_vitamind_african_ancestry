# Program: 06_create_ld_prune_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Create queue files for LD prunning with condor

library(magrittr)
library(tidyverse)

workdir <- here::here("extracted_data", "SCCS", "imputation",
  "merge_impute")

bedfiles <- list.files(file.path(workdir, "qc_filtered"),
  pattern = "bed")

queue <- tibble(file = bedfiles) %>%
  mutate(
    file = str_remove(file, ".bed"),
    outfile = file.path(workdir, "pca_tmp", file),
    file = file.path(workdir, "qc_filtered", file))

queue %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs",
      "condor", "queues", "2020_09_15_pca_ldprune.csv"),
    col_names = FALSE)

queue %>%
  mutate(
    queue_file = str_c(outfile, ".prune.in")) %>%
  select(file, queue_file, outfile) %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs",
      "condor", "queues", "2020_09_15_pca_ldprune_subset.csv"),
    col_names = FALSE)

  
