
library(magrittr)
library(tidyverse)

princomp <- here::here("results/2021_07_01_ukbiobank_stratified_gwas",
  "plink2/ukbiobank_pca.eigenvec") %>%
  vroom::vroom()

princomp %<>%
  dplyr::rename_all(list(snakecase::to_snake_case))

# remove samples from email
to_remove <- here::here("results/2021_07_01_ukbiobank_stratified_gwas",
  "csv/sample_to_remove.csv") %>%
  readr::read_csv(col_names = FALSE) %>%
  dplyr::mutate(iid = stringr::str_c(X1, X1, sep = "_")) %>%
  dplyr::select(iid)

princomp %>%
  dplyr::anti_join(to_remove, by = "iid") %>%
  qs::qsave(here::here("results/2021_07_01_ukbiobank_stratified_gwas",
    "qs", "ukbiobank_princomp.qs"))
