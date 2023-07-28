library(magrittr)
library(tidyverse)

folder <- "2021_07_01_ukbiobank_stratified_gwas"

uk_data <- here::here("results", folder, "qs", "ukbiobank_grant_vars.qs") %>%
  qs::qread()

princomp <- here::here("results", folder, "qs", "ukbiobank_princomp.qs") %>%
  qs::qread()
  
princomp %<>%
  dplyr::mutate(
    eid = stringr::str_split(iid, "\\_") %>%
      purrr::map_chr(1) %>%
      as.integer()) %>%
  dplyr::select(-iid)
  
uk_data %<>%
  dplyr::inner_join(princomp, by = "eid") %>%
  filter(!is.na(vitamin_d))

uk_data %<>%
  dplyr::mutate(
    dplyr::across(starts_with("pc_"),
      list(~ . / sd(.)), .names = "{.col}"))

uk_data %<>%
    dplyr::rename(id = eid) %>%
    dplyr::mutate(fid = id) %>%
    dplyr::rename(`#FID` = fid, IID = id) %>%
    dplyr::select(`#FID`, IID, tidyselect::everything())

uk_data %>%
  dplyr::select(`#FID`,
    IID, age, bmi, tidyselect::starts_with("pc")) %>%
  readr::write_delim(here::here("results", folder, "plink2",
    "ukbiobank_covariates.txt"), delim = " ")

uk_data %>%
  dplyr::mutate(vitd = log(vitamin_d)) %>%
  dplyr::select(`#FID`, IID, vitd) %>%
  readr::write_delim(here::here("results", folder, "plink2",
    "ukbiobank_logvitd.txt"), delim = " ")
