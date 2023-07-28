
library(magrittr)
library(tidyverse)
library(vroom)

folder <- "2021_07_01_ukbiobank_stratified_gwas"
stats_dir <- here::here("results", folder, "plink2/stats")

all_files <- list.files(stats_dir, pattern = "glm")

stats <- tibble::tibble(
  file = all_files) %>%
  dplyr::mutate(
    stats = map(file.path(stats_dir, file), vroom),
    chr = str_split(file, regex("[\\_|\\.]")),
    chr = map_chr(chr, ~ str_subset(., "chr")),
    pop = case_when(
      str_detect(file, "asian") ~ "asian",
      str_detect(file, "afro") ~ "afrocaribbean",
      TRUE ~ "white"),
    stats = map(stats, dplyr::rename_all,
      snakecase::to_snake_case)) %>%
    select(-file)

count_errcode <- function(stats) {

  stats %>%
    filter(is.na(t_stat)) %>%
    count(errcode)

}

stats %<>%
  dplyr::mutate(
    errs = map(stats, count_errcode),
    stats = map(stats, filter, !is.na(t_stat)),
    stats = map(stats, select, -errcode))


stats %>%
  filter(pop == "white") %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukbiobank_strat_gwas_white.qs"))

stats %>%
  filter(pop == "asian") %>%
  qs::qsave(
    here::here("results", folder, "qs", "ukbiobank_strat_gwas_asian.qs"))
    
stats %>%
  filter(pop == "afrocaribbean") %>%
  qs::qsave(
    here::here("results", folder, "qs",
      "ukbiobank_strat_gwas_afrocaribbean.qs"))