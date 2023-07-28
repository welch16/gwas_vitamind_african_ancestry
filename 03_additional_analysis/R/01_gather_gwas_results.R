
# we need to gather the 3 datasets:
# 
# UK Biobank: vdbp 
# SCCS: vdbp and hydroxyvitamin d

pacman::p_load(here, qs, magrittr, tidyverse, vroom, data.table)

ukbiobank_dir <- here::here("extracted_data/UKBiobank/gwas_results/log_vitd",
  "sex_bmi_age")

sccs_dir <- here::here("extracted_data/SCCS/gwas_results_bmi_age")

read_file <- function(file) {

  vroom::vroom(file, skip = 1,
    col_names = c("chrom", "pos", "id", "ref", "alt", "a_1", "test", "obs_ct",
      "beta", "se", "t_stat"),
  col_types = c("chrom" = "i", "pos" = "i", "id" = "c", "ref" = "c",
    "alt" = "c", "a_1" = "c", "test" = "c", "obs_ct" = "i", "beta" = "d",
    "se" = "d", "t_stat" = "d"),
    col_select = c("chrom", "pos", "id", "ref", "alt", "a_1", "beta", "se",
      "t_stat"))
}


ukbiobank_hydroxy <- list.files(ukbiobank_dir,
  pattern = "linear", full.names = TRUE) %>%
  str_subset("id", negate = TRUE) %>%
  purrr::map(read_file) %>%
  bind_rows()
  
ukbiobank_hydroxy %>%
  qs::qsave(
    here::here("results", "zz_methods", "qs", "ukbiobank_vitd_stats.qs"))
rm(ukbiobank_hydroxy)
gc()

sccs_hydroxy <- list.files(sccs_dir, pattern = "rsq50", full.names = TRUE) %>%
  str_subset(".linear") %>%
  str_subset("id", negate = TRUE) %>%
  str_subset("log") %>%
  str_subset("vit25") %>%
  purrr::map(read_file) %>%
  bind_rows()

sccs_hydroxy %>%
  qs::qsave(here::here("results", "zz_methods", "qs", "sccs_hydroxy_stats.qs"))
rm(sccs_hydroxy)
gc()


sccs_vdbp <- list.files(sccs_dir, pattern = "rsq50", full.names = TRUE) %>%
  str_subset(".linear") %>%
  str_subset("id", negate = TRUE) %>%
  str_subset("log") %>%
  str_subset("vitdbp") %>%
  purrr::map(read_file) %>%
  bind_rows()

sccs_vdbp %>%
  qs::qsave(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs"))
rm(sccs_vdbp)
gc()
