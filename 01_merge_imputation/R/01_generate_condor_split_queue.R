# Program: 01_generate_condor_split_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-01
# Purpose: Generate queue for splitting samples using condor

library(magrittr)
library(tidyverse)

data_dir <- here::here("extracted_data/SCCS/Shaneda_SCCS_data")

files <- list.files(data_dir,
  full.names = TRUE, recursive = TRUE, pattern = "bed")

queue <- tibble(bfile = files)

queue %<>%
  mutate(bfile = str_remove(bfile, ".bed")) %>%
  crossing(chrom = seq_len(22)) %>%
  mutate(
    outfile = glue::glue(
      "{outdir}/{prefix}_chr{chrom}",
      chrom = chrom,
      prefix = basename(bfile),
      outdir = here::here(
        "extracted_data/SCCS/imputation/vcf_by_chrom/")))
    
queue %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor",
      "queues", "2020_09_01_vcf_by_chrom.csv"),
    col_names = FALSE)