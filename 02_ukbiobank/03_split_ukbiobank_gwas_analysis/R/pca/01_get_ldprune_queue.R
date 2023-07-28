library(magrittr)
library(tidyverse)

raw_dir <- here::here("raw_data/UKBiobank/2020_02/genotype")
bgen_files <- list.files(raw_dir, pattern = "bgen", full.names = TRUE)
bgen_files %<>%
  stringr::str_subset("bgi", negate = TRUE) %>%
  stringr::str_subset(regex("bgen$"))

out_dir <- file.path(raw_dir, "ldprune")

fs::dir_create(out_dir)

out <- tibble(bgen_files) %>%
  mutate(
    chrom = basename(bgen_files),
    chrom = str_remove(chrom, "ukb_imp_"),
    chrom = str_remove(chrom, "_v3.bgen"),
    outfile = file.path(out_dir, chrom)) %>%
  dplyr::select(-chrom)

out %>%
  readr::write_csv(
    here::here("results", "2021_07_01_ukbiobank_stratified_gwas",
      "condor", "pca", "2021_07_ldprune_queue.csv"), col_names = FALSE)
