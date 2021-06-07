# Program: 05_check_imputation_quality.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-15
# Purpose: Check data imputation quality scores

library(magrittr)
library(tidyverse)
library(vroom)
library(data.table)

data_dir <- here::here("extracted_data/SCCS/imputation/imputed_by_chrom")

info_file <- list.files(data_dir, pattern = "info")


data <- tibble::tibble(file = file.path(data_dir, info_file))
data %<>%
  mutate(
    data = map(file, vroom::vroom))

data %<>%
  mutate(
    file = basename(file),
    chr = str_split(file, "\\_"),
    chr = map_chr(chr, ~ .[length(.)]),
    chr = str_remove(chr, ".info"),
    group = str_split(file, "_chr"),
    group = map_chr(group, 1))

data %<>% select(-file, -data, data)
data %<>%
  mutate(
    data = map(data, rename_all, list(snakecase::to_snake_case)))

# quick test to see how to pick a consensus of snps to select 

summarize_by_chr <- function(my_chr, data) {

  chr_data <- data %>%
    filter(chr == my_chr) %>%
    unnest(cols = c(data))

  chr_data %<>%
    as.data.table()

  return(chr_data[,
    .(med_avg_call = median(avg_call),
      med_rsq = median(rsq),
      med_maf = median(maf)), by = .(snp)])
}

chr <- data %>%
  distinct(chr)

chr %<>%
  mutate(
    summary = map(chr, summarize_by_chr, data))

chr %>%
  unnest(cols = c(summary)) %>%
  qs::qsave(here::here("results", "2020_09_01_merge_impute_sccs", "qs",
    "imputation_summary_by_chrom.qs"))

chr <-
  qs::qread(here::here("results", "2020_09_01_merge_impute_sccs", "qs",
    "imputation_summary_by_chrom.qs"))

chr %<>%
  nest(summary = c(snp, med_avg_call, med_rsq, med_maf))

# the other alternative would be to try  0.05
chr %>%
  mutate(
    outfile = glue::glue("{dir}/{chrom}_filtered_variants_rsq50.txt",
      dir = here::here("results", "2020_09_01_merge_impute_sccs",
        "variants"),
      chrom = chr),
    map(summary, filter, med_maf >= 0.01 & med_rsq >= .5) %>%
      map(select, snp) %>%
      map2(outfile, write_csv, col_names = FALSE))

chr %>%
  mutate(
    outfile = glue::glue("{dir}/{chrom}_filtered_variants_rsq30.txt",
      dir = here::here("results", "2020_09_01_merge_impute_sccs",
        "variants"),
      chrom = chr),
    map(summary, filter, med_maf >= 0.01 & med_rsq >= .3) %>%
      map(select, snp) %>%
      map2(outfile, write_csv, col_names = FALSE))
