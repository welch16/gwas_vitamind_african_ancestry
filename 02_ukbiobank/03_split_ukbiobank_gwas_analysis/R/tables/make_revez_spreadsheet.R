library(magrittr)
library(tidyverse)
library(data.table)
library(readxl)
library(vroom)

folder <- "2021_07_01_ukbiobank_stratified_gwas"


revez <- here::here("results/2021_06_04_25hydroxy_gwas_analysis/xls",
  "Table 1_SCCS_Africanancestry.xlsx") %>%
  readxl::read_excel(sheet = 1, skip = 3) %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  filter(snp != "SNP") %>%
  filter(!is.na(snp))


results <- here::here("results", folder, "qs", "ukb_strat_revez.qs") %>%
  qs::qread() %>%
  unnest(cols = c(stats))

results %<>%
  select(-test, -errcode, -ref, -alt, -a_1) %>%
  pivot_wider(
    names_glue = "{pop}_{.value}",
    names_from  = pop,
    names_repair = "minimal",
    values_from = one_of(c("obs_ct", "beta", "se", "t_stat",
      "p", "alt_freqs", "f_miss", "hwe_p"))) %>%
  dplyr::select(
    chr, snp, chrom, pos,
    starts_with("afro"),
    starts_with("asian"),
    starts_with("white"))

revez_snps <- revez %>%
  select(snp)

revez_snps %>%
  left_join(results, by = "snp") %>%
  readr::write_csv(
    here::here("results", folder, "csv", "revez_snps.csv"))
  
results <- here::here("results", folder, "qs", "ukb_strat_revez.qs") %>%
  qs::qread() %>%
  unnest(cols = c(stats))


revez_snps %>%
  left_join(
  results %>%
    filter(pop == "afrocaribbean") %>%
    select(snp, ref, alt, a_1), by = "snp") %>%
  readr::write_csv(
    here::here("results", folder, "csv", "revez_alleles_afrocaribbean.csv"))
  
revez_snps %>%
  left_join(
  results %>%
    filter(pop == "asian") %>%
    select(snp, ref, alt, a_1), by = "snp") %>%
  readr::write_csv(
    here::here("results", folder, "csv", "revez_alleles_asian.csv"))
  
revez_snps %>%
  left_join(
    results %>%
    filter(pop == "white") %>%
    select(snp, ref, alt, a_1), by = "snp") %>%
  readr::write_csv(
    here::here("results", folder, "csv", "revez_alleles_white.csv"))
  

