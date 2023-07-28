library(magrittr)
library(tidyverse)

# this script saves txt files to filter the ukbiobank bgen files


# sample file
sample_file <- here::here("results/2021_07_01_ukbiobank_stratified_gwas",
  "plink2/samples/afrocaribbean_ukbiobank_samples.txt") %>%
  vroom::vroom()

sample_file %<>%
  dplyr::select(IID) %>%
  dplyr::rename(`0` = IID)

tibble::tibble("ID") %>%
  readr::write_csv(
    here::here("results", "2022_08_04_ukbiobank_finemapping",
      "vcf", "afrocaribbean_samplefile.txt"),
    col_names = FALSE)

sample_file %>%
  readr::write_csv(
    here::here("results", "2022_08_04_ukbiobank_finemapping",
      "vcf", "afrocaribbean_samplefile.sample"),
    col_names = TRUE, append = TRUE)

# variant file
stats_files <- here::here("results",
  "2021_07_01_ukbiobank_stratified_gwas/plink2/stats") %>%
  list.files(full.names = TRUE) %>%
  stringr::str_subset("log", negate = TRUE) %>%
  stringr::str_subset("afro")


get_chr <- function(file) {

  file %>%
    basename() %>%
    stringr::str_split("\\_") %>%
    purrr::map_chr(stringr::str_subset, "chr") %>%
    stringr::str_split("\\.") %>%
    purrr::map_chr(stringr::str_subset, "chr")

}

get_ext <- function(file) {

  file %>%
    basename() %>%
    stringr::str_split("\\.") %>%
    purrr::map_chr(~ .[length(.)])

}


stats <-
  tibble::tibble(file = stats_files) %>%
  dplyr::mutate(
    chr = get_chr(file),
    ext = get_ext(file)) %>%
  tidyr::pivot_wider(names_from = ext, values_from = file) %>%
  dplyr::select(-linear)


read_file <- function(file) {

  file %>%
    vroom::vroom() %>%
    dplyr::rename_with(snakecase::to_snake_case)
}

pick_variants <- function(afreq, hardy, smiss, vmiss) {

  `%<>%` <- magrittr::`%<>%`
  # read the files
  afreq %<>% read_file()
  hardy %<>% read_file()
  vmiss %<>% read_file()

  afreq %<>%
    dplyr::filter(pmin(alt_freqs, 1 - alt_freqs) >= 0.01) %>%
    dplyr::select(id, alt_freqs)

  hardy %<>%
    dplyr::filter(p > 1e-6) %>%
    dplyr::select(id, p) %>%
    dplyr::rename(hwe_pval = p)

  vmiss %<>%
    dplyr::filter(f_miss <= 0.05) %>%
    dplyr::select(id, f_miss)

  list(afreq, hardy, vmiss) %>%
    purrr::reduce(dplyr::inner_join, by = "id")

}

stats0 <- stats %>%
  dplyr::filter(chr %in% c("chr3", "chr4")) %>%
  dplyr::mutate(
    variants = purrr::pmap(list(afreq, hardy, smiss, vmiss),
      pick_variants))

stats0 %>%
  dplyr::select(-afreq, -hardy, -smiss, -vmiss) %>%
  dplyr::mutate(
    variants = purrr::map(variants, dplyr::filter, alt_freqs >= .05),
    variants = purrr::map(variants, dplyr::select, id),
    outfile = here::here("results", "2022_08_04_ukbiobank_finemapping",
      "vcf", glue::glue("variants_{chr}_afrocarib.txt")),
    outfile = purrr::walk2(
      variants, outfile, readr::write_delim, delim = " ",
        col_names = FALSE))


folder <- "2021_07_01_ukbiobank_stratified_gwas"
ukb_afroc_gwas <-
  here::here("results", folder, "plink2/ukb_gwas_filtered_afrocaribbean.qs")

ukb_afroc_gwas %<>% qs::qread()

add_stats <- function(stats, afreq, hardy, smiss, vmiss) {

  message(unique(stats$chrom))
  `%<>%` <- magrittr::`%<>%`
  # read the files
  afreq %<>% read_file()
  hardy %<>% read_file()
  vmiss %<>% read_file()

  afreq %<>%
    dplyr::filter(pmin(alt_freqs, 1 - alt_freqs) >= 0.05) %>%
    dplyr::select(id, alt_freqs)

  stats %<>%
    dplyr::inner_join(afreq, by = "id")

  hardy %<>%
    dplyr::filter(p > 1e-6) %>%
    dplyr::select(id, p) %>%
    dplyr::rename(hwe_pval = p)

  stats %<>%
    dplyr::inner_join(hardy, by = "id")

  vmiss %<>%
    dplyr::filter(f_miss <= 0.05) %>%
    dplyr::select(id, f_miss)

  stats %>%
    dplyr::inner_join(vmiss, by = "id")

}

ukb_afroc_gwas %<>%
  tidyr::nest(stats = -c(chr, pop))

ukb_afroc_gwas %<>%
  dplyr::inner_join(stats, by = "chr") %>%
  dplyr::mutate(
    stats = purrr::pmap(
      list(stats, afreq, hardy, smiss, vmiss), add_stats))


out <- ukb_afroc_gwas %>%
  dplyr::select(stats) %>%
  tidyr::unnest(cols = c(stats))

out %>%
  qs::qsave(
    here::here("results/2022_08_04_ukbiobank_finemapping/qs",
      "afroc_ukb_stats_filter.qs"))

sign_snps <- ukb_afroc_gwas %>%
  dplyr::filter(p <= 5e-8)


stats1 <- stats0 %>%
  dplyr::select(chr, variants) %>%
  tidyr::unnest(cols = c(variants))

sign_snps %>%
  dplyr::left_join(
    stats1, by = c("chr", "id")) %>%
  as.data.frame()

ukb_afroc_gwas
