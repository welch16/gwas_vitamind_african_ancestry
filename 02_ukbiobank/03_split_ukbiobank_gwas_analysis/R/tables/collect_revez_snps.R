library(magrittr)
library(tidyverse)
library(data.table)
library(readxl)
library(vroom)


revez <- here::here("results/2021_06_04_25hydroxy_gwas_analysis/xls",
  "Table 1_SCCS_Africanancestry.xlsx") %>%
  readxl::read_excel(sheet = 1, skip = 3) %>%
  dplyr::rename_with(snakecase::to_snake_case)

revez_snps <- dplyr::select(revez, snp)

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

load_data_file <- function(chr, pop, all_files, revez_snps) {

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
    dplyr::select(id, alt_freqs)

  miss <- my_files %>%
    str_subset("vmiss") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    dplyr::select(id, f_miss)

  hwe <- my_files %>%
    str_subset("hardy") %>%
    vroom::vroom() %>%
    rename_with(snakecase::to_snake_case) %>%
    dplyr::select(id, p) %>%
    dplyr::rename(hwe_p = p)

  revez_snps %>%
    dplyr::left_join(stats, by = c(snp = "id")) %>%
    dplyr::left_join(maf, by = c(snp = "id")) %>%
    dplyr::left_join(miss, by = c(snp = "id")) %>%
    dplyr::left_join(hwe, by = c(snp = "id"))

}

results %<>%
  mutate(
    stats = map2(chr, pop, load_data_file,
      all_files, revez_snps))

results %>%
  qs::qsave(here::here("results", folder, "qs", "ukb_strat_revez.qs"))
