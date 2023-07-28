library(magrittr)
library(tidyverse)
library(UWCCC.GWAS)


data_dir <- here::here("extracted_data/UKBiobank/gwas_results/strat_pop")

gwas_results <- tibble::tibble(file = list.files(data_dir))
gwas_results %<>%
  dplyr::filter(! stringr::str_detect(file, stringr::regex("log$"))) %>%
  dplyr::filter(! stringr::str_detect(file, stringr::regex("id$")))

gwas_results %<>%
  dplyr::mutate(
    aux = stringr::str_split(file, "\\_"),
    chrom = purrr::map_chr(aux, 2),
    ethnicity = purrr::map_chr(aux, 3)) %>%
  dplyr::select(-aux)

folder <- "2021_01_06_ukbiobank_split_populations"

get_results <- function(files) {

  l <- purrr::map(files, vroom::vroom) %>%
    purrr::map(filter, abs(T_STAT) <= 38)

  data.table::rbindlist(l) %>%
    dplyr::rename_with(snakecase::to_snake_case)

}

gwas_results %<>%
  dplyr::mutate(file = file.path(data_dir, file))


# african
gwas_results %>%
  dplyr::filter(ethnicity == "african") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_african.qs"))


# afrocaribbean
gwas_results %>%
  dplyr::filter(ethnicity == "afrocaribbean") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_afrocaribbean.qs"))

# caribbean
gwas_results %>%
  dplyr::filter(ethnicity == "caribbean") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_caribbean.qs"))

# chinese
gwas_results %>%
  dplyr::filter(ethnicity == "chinese") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_chinese.qs"))
  
# european
gwas_results %>%
  dplyr::filter(ethnicity == "european") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_european.qs"))

# south asian
gwas_results %>%
  dplyr::filter(ethnicity == "south") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_southasian.qs"))

# asian
gwas_results %>%
  dplyr::filter(ethnicity == "asian") %>%
  dplyr::pull(file) %>%
  get_results() %>%
  tibble::as_tibble() %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukb_vitdbp_gwas_asian.qs"))
