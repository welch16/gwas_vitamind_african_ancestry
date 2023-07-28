
library(magrittr)
library(tidyverse)
library(qs)


# ukbiobank_vdbp <-
#   qs::qread(
#     here::here("results", "zz_methods", "qs", "ukbiobank_vitd_stats.qs"))

# ukbiobank_vdbp %<>%
#   dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))


# ukbiobank_vdbp %>%
#   readr::write_csv(here::here("results", "zz_methods", "csv",
#     "stats", "ukbiobank_hydroxy25_vitd_stat.csv"))

sccs_hydroxy <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_hydroxy_stats.qs"))

sccs_hydroxy %<>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))

sccs_hydroxy %>%
  readr::write_csv(here::here("results", "zz_methods", "csv",
    "stats", "sccs_hydroxy25_vitd_stat.csv"))

sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs"))

sccs_vdbp %<>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))

sccs_vdbp %>%
  readr::write_csv(here::here("results", "zz_methods", "csv",
    "stats", "sccs_vitdbp_stat.csv"))

