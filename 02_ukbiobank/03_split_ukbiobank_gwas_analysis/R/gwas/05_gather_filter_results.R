library(magrittr)
library(vroom)
library(tidyverse)

folder <- "2021_07_01_ukbiobank_stratified_gwas"
work_dir <- here::here("results", folder, "plink2")

all_files <- list.files(work_dir, full.names = TRUE, recursive = TRUE)

chroms <- basename(all_files) %>%
  str_subset("chr") %>%
  str_split(regex("[\\.|\\_]")) %>%
  map_chr(str_subset, "chr") %>%
  unique()
pops <- c("afrocaribbean", "asian", "white")

results <- crossing(chr = chroms, pop = pops)

load_data_file <- function(chr, pop,
  hwe_thr, maf_thr, callrate_thr, all_files) {

  my_files <- all_files %>%
    str_subset("chr") %>%
    str_subset(str_c(chr, ".txt")) %>%
    str_subset(pop)

  stats <- my_files %>%
    str_subset("glm.linear") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    filter(!is.na(t_stat))

  maf <- my_files %>%
    str_subset("afreq") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    filter(alt_freqs >= maf_thr) %>%
    filter(alt_freqs <= 1 - maf_thr) %>%
    select(id)

  miss <- my_files %>%
    str_subset("vmiss") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    filter(f_miss <= 1 - callrate_thr) %>%
    select(id)

  hwe <- my_files %>%
    str_subset("hardy") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    filter(p > hwe_thr) %>%
    select(id)

  stats %>%
    inner_join(hwe, by = "id") %>%
    inner_join(miss, by = "id") %>%
    inner_join(maf, by = "id") %>%
    select(chrom, pos, id, ref, alt, a_1, t_stat, p)

}

results %<>%
  mutate(
    stats = map2(chr, pop, load_data_file,
      hwe_thr = 1e-5,
      maf_thr = 0.01,
      callrate_thr = 0.95,
      all_files))

results %>%
  filter(pop == "asian") %>%
  unnest(cols = c(stats)) %>%
  qs::qsave(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_asian.qs"))

results %>%
  filter(pop == "white") %>%
  unnest(cols = c(stats)) %>%
  qs::qsave(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_white.qs"))

results %>%
  filter(pop == "afrocaribbean") %>%
  unnest(cols = c(stats)) %>%
  qs::qsave(here::here("results", folder, "plink2",
    "ukb_gwas_filtered_afrocaribbean.qs"))
