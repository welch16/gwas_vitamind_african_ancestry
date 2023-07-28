
pacman::p_load(
  magrittr,
  tidyverse)

input_dir <- here::here("raw_data", "UKBiobank", "2020_02", "gwas")
output_dir <- here::here("extracted_data", "UKBiobank", "plink2")

queue <- tibble::tibble(
  infile = list.files(input_dir, full.names = TRUE)) %>%
  dplyr::mutate(
    chr = basename(infile))

queue %<>%
  dplyr::filter(str_detect(infile, "bgen")) %>%
  dplyr::filter(!str_detect(infile, "bgi")) %>%
  dplyr::mutate(
    chr = stringr::str_remove(chr, "ukb_imp_"),
    chr = str_split(chr, "\\_"),
    chr = map_chr(chr, 1),
    outfile = str_c("ukb_", chr, "_colorectal_vitD_gwas"),
    outfile = file.path(output_dir, outfile)) %>%
  dplyr::select(infile, outfile, chr)

queue %>%
  readr::write_csv(
    here::here("results", "2020_06_08_gwas_analysis", "condor",
    "2020_06-filter_queue.csv"), col_names = FALSE)