# Program: pre07_fix_fam_files.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-10-20
# Purpose: Fix `fam` files for plink / plink2 use


library(magrittr)
library(tidyverse)

datadir <- here::here("extracted_data", "SCCS")
vitd_data <- file.path(datadir, "vitamind_phenotypes.txt") %>%
  readr::read_tsv()

workdir <- here::here("extracted_data/SCCS/imputation/merge_impute/qc_filtered")
fam_files <- list.files(workdir, pattern = "\\.fam")

# read fam files
fams <- map(file.path(workdir, fam_files), read_tsv, col_names = FALSE)

vitd_data %<>%
  rename_all(list(snakecase::to_snake_case)) %>%
  select(id_number, sex) %>%
  mutate(
    sex = if_else(sex == "M", 1, 2),
    id_number = str_c(id_number, id_number, sep = "_"))

fams %<>%
  map(left_join, vitd_data, by = c(X2 = "id_number")) %>%
  map(mutate, X5 = sex) %>%
  map(select, -sex)

map2(
  fams,
  file.path(workdir, fam_files), write_tsv, col_names = FALSE)