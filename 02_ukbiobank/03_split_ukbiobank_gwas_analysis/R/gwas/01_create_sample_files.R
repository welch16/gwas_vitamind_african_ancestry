
library(magrittr)
library(tidyverse)

folder <- "2021_07_01_ukbiobank_stratified_gwas"

uk_data <- here::here("results", folder, "qs", "ukbiobank_grant_vars.qs") %>%
  qs::qread()

save_file <- function(data, filename) {

  data %>%
    dplyr::rename(id = eid) %>%
    dplyr::select(id) %>%
    dplyr::mutate(fid = id) %>%
    rlang::set_names("#FID", "IID") %>%
    readr::write_tsv(filename)

}

princomp <- here::here("results", folder, "qs", "ukbiobank_princomp.qs") %>%
  qs::qread()
  
princomp %<>%
  dplyr::mutate(
    eid = stringr::str_split(iid, "\\_") %>%
      purrr::map_chr(1) %>%
      as.integer()) %>%
  dplyr::select(-iid)
  
uk_data %<>%
  dplyr::inner_join(princomp, by = "eid")


uk_data %>%
  filter(!is.na(vitamin_d)) %>%
  count(pop)

uk_data %>%
  dplyr::filter(pop == "afrocaribbean") %>%
  dplyr::filter(!is.na(vitamin_d)) %>%
  save_file(
    here::here("results", folder, "plink2", "samples",
      "afrocaribbean_ukbiobank_samples.txt"))

uk_data %>%
  dplyr::filter(pop == "asian") %>%
  dplyr::filter(!is.na(vitamin_d)) %>%
  save_file(
    here::here("results", folder, "plink2", "samples",
      "asian_ukbiobank_samples.txt"))

uk_data %>%
  dplyr::filter(pop == "white") %>%
  dplyr::filter(!is.na(vitamin_d)) %>%
  save_file(
    here::here("results", folder, "plink2", "samples",
      "white_ukbiobank_samples.txt"))
