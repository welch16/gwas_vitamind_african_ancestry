library(magrittr)
library(tidyverse)

folder <- "2021_07_01_ukbiobank_stratified_gwas"
geno_dir <- here::here("raw_data/UKBiobank/2020_02/genotype")
out_dir <- here::here("results", folder, "plink2", "stats")

fs::dir_create(out_dir)

bgen_files <- list.files(geno_dir, pattern = "\\.bgen", full.names = TRUE) %>%
  str_subset("bgi", negate = TRUE)

sample_dir <- here::here("results", folder, "plink2", "samples")

samples <- tibble::tibble(sample = list.files(sample_dir, full.names = TRUE))

queue <- tibble::tibble(file = bgen_files) %>%
  dplyr::mutate(
    chr = str_split(basename(file), "\\_"),
    chr = map_chr(chr, 3))

samples %<>%
  dplyr::mutate(
    pop = str_split(basename(sample), "\\_"),
    pop = map_chr(pop, 1))
    
queue %<>%
  crossing(samples) %>%
  dplyr::mutate(
    outfile = glue::glue("{outdir}/ukbiobank_{pop}_gwas_{chr}.txt",
      outdir = out_dir),
    outfile = as.character(outfile))

queue %>%
  dplyr::select(-chr, -pop) %>%
  readr::write_csv(
    here::here("results", folder, "condor", "gwas", "ukb_gwas_queue.csv"),
    col_names = FALSE)
