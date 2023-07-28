
pacman::p_load(
  magrittr,
  tidyverse)

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
    chr = map_chr(chr, 2),
    outfile = str_c("ukb_", chr, "_colorectal_vitD_gwas_out"),
    adj_outfile = str_c("ukb_", chr, "_colorectal_vitD_gwas_adjust"),
    adj_outfile = file.path(output_dir, adj_outfile)) %>%
  dplyr::select(infile, outfile, chr)

queue %>%
  readr::write_csv(
    here::here("results", "2020_06_08_gwas_analysis", "condor",
    "2020_06-gwas_queue.csv"), col_names = FALSE)
