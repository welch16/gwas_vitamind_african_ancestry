# Program: 03a_select_samples_per_cohort.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2021-04-12
# Purpose: Picks samples with highest imputation quality when repeated in
#   multiple cohorts

library(magrittr)
library(tidyverse)
library(vroom)

work_dir <- here::here("extracted_data", "SCCS")

files <- list.files(work_dir, recursive = TRUE, pattern = "imiss")

clean_samples <- function(samples) {

  samples %>%
    mutate(x = str_squish(x)) %>%
    separate("x", into = c("fid", "iid", "miss_pheno", "n_miss", "n_geno",
      "f_miss"), sep = " ")

}


data <- tibble(
  file = file.path(work_dir, files),
  cohort = map_chr(str_split(files, "\\/"), 2)) %>%
  mutate(
    data = map(file, vroom, delim = "\t", skip = 1,
      col_names = "x"),
    data = map(data, clean_samples)) %>%
  select(-file)

data %<>% unnest(cols = c(data))

unique_samples <- data %>%
  mutate(
    n_miss = as.numeric(n_miss),
    f_miss = as.numeric(f_miss)) %>%
  group_by(fid, iid) %>%
  summarize(
    call_rate = max(1 - f_miss),
    n_miss = min(n_miss),
    cohort = cohort[ which.max(1 - f_miss)], .groups = "drop")
    
unique_samples %<>%
  nest(samples = c(fid, iid, call_rate, n_miss))
unique_samples %<>%
  mutate(
    outfile = glue::glue(
      "{dir}/{cohort}_bestcallrate_samples.txt",
      dir = here::here("results", "2020_09_01_merge_impute_sccs",
        "samples"),
      cohort = cohort),
    samples = map(samples, mutate, ss = str_c(fid, iid, sep = "_")),
    samples = map(samples, select, ss) %>%
      map2(outfile, write_csv, col_names = FALSE))
