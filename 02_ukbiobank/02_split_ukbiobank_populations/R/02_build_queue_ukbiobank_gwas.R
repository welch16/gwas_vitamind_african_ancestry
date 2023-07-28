library(magrittr)
library(tidyverse)

input_dir <- here::here("extracted_data", "UKBiobank", "plink2")
output_dir <- here::here("extracted_data", "UKBiobank", "plink2", "analysis")

queue <- tibble::tibble(
  infile = list.files(input_dir, full.names = TRUE)) %>%
  dplyr::mutate(
    chr = basename(infile))

queue %<>%
  dplyr::filter(str_detect(infile, "pgen")) %>%
  dplyr::mutate(
    infile = stringr::str_remove(infile, ".pgen"),
    chr = stringr::str_remove(chr, "ukb_imp_"),
    chr = str_split(chr, "\\_"),
    chr = map_chr(chr, 2))
    
samples <- list.files(here::here("results",
  "2021_01_06_ukbiobank_split_populations", "plink2"), full.names = TRUE)

samples <- tibble::tibble(samples) %>%
  dplyr::mutate(
    prefix = basename(samples),
    prefix = str_remove(prefix, "_samples.tsv")) %>%
  dplyr::filter(!stringr::str_detect(prefix, "chr"))

queue %<>%
  tidyr::crossing(samples) %>%
  dplyr::mutate(
    outfile = glue::glue("ukb_{chr}_{prefix}_colorectal_vitd_gwas")) %>%
  dplyr::select(infile, samples, outfile, chr, prefix)

queue %>%
  readr::write_csv(
    here::here("results", "2021_01_06_ukbiobank_split_populations", "condor",
    "queues", "2021_01-ukbiobank_gwas_queue_pops.csv"), col_names = FALSE)
