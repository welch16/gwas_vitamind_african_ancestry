# Program: 05c_bed_to_vcf_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Write queue to convert from bed to vcf with condor

library(magrittr)
library(tidyverse)

workdir <- here::here("extracted_data", "SCCS", "imputation",
  "merge_impute", "qc_filtered")

bedfiles <- list.files(workdir, pattern = "bed")

queue <- tibble(file = bedfiles) %>%
  mutate(
    file = file.path(workdir, file),
    file = str_remove(file, ".bed"))

queue %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs",
      "condor", "queues", "2020_09_15_convert_bed2vcf.csv"),
    col_names = FALSE)