library(magrittr)
library(tidyverse)

source(here::here("src/gwas_funs_redo.R"))


qc_stats <- here::here("raw_data/UKBiobank/2020_02/summary_stats",
  "ukb_imp_chr4_colosubset_stats.txt")

# load summary statistics
folder <- "2022_06_29_prepare_imputation"
workdir <- here::here("results", folder, "qs", "gwas_clean")
all_files <- list.files(workdir, full.names = TRUE, pattern = "gwas_stats")

stats1 <- all_files %>%
  stringr::str_subset("covar_bmi_age") %>%
  stringr::str_subset("vit_2") %>%
  stringr::str_subset("chr4") %>%
  purrr::map(qs::qread) %>%
  bind_rows()

table_snp_keys <- c("4:72618334", "4:72627797", "4:72696465",
  "4:72819988", "4:72608115")

snps <- purrr::map(table_snp_keys,
  ~ dplyr::filter(stats1, stringr::str_detect(id, .))) %>%
  dplyr::bind_rows()

# Bioconductor part
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)
library(Matrix)

window <- 1

snps %<>%
  dplyr::mutate(id2 = id) %>%
  tidyr::nest(data = -c(id2)) %>%
  dplyr::mutate(
    gr = purrr::map(data, create_gr),
    gr = purrr::map(gr, IRanges::resize, window, fix = "center"))

vcf_files <- here::here("results", "2022_06_29_prepare_imputation",
  "impute_subset") %>%
  list.files(full.names = TRUE, pattern = ".dose.vcf.gz") %>%
  stringr::str_subset(".csi", negate = TRUE) %>%
  stringr::str_subset(".tbi", negate = TRUE)

snps %<>%
  dplyr::mutate(
    geno = purrr::map(gr, get_genotype, vcf_files, stats1))

# add covariates and response
folder <- "2022_06_29_prepare_imputation"
sampledir <- here::here("results", folder, "tsv", "samples")

samples <- tibble::tibble(
  file = list.files(sampledir, full.names = TRUE),
  platform = str_remove(basename(file), "_best_call_rates.tsv")) %>%
  dplyr::mutate(
    samples = purrr::map(file, readr::read_tsv, col_names = "iid")) %>%
  dplyr::select(-file)

data_dir <- here::here("extracted_data", "SCCS")
vitd_data <- file.path(data_dir, "vitamind_phenotypes.txt") %>%
  readr::read_tsv()

vitd_data %<>%
  dplyr::rename_all(list(snakecase::to_snake_case)) %>%
  dplyr::mutate(iid = str_c(id_number, id_number, sep = "_"))

samples %<>%
  dplyr::mutate(
    samples = purrr::map(samples, dplyr::left_join, vitd_data, by = "iid")) %>%
  tidyr::unnest(cols = c(samples))

# load pca
pcas <- here::here("results", folder, "ld_prune", "sccs_pca.eigenvec") %>%
  vroom::vroom() %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::mutate(
    across(where(is.numeric), list(~ . / sd(.)), .names = "{.col}"))


# one phenotype files
pheno <- samples %>%
  dplyr::select(iid, vit_25_ohd, vdbp) %>%
  dplyr::mutate(across(-iid, list(~ log1p(.)), .names = "{.col}")) %>%
  dplyr::select(iid, vit_25_ohd)

# 1 covariates file
covar <- samples %>%
  dplyr::select(iid, platform, bmi, enrollment_age, sex) %>%
  dplyr::mutate(
    age = case_when(
      enrollment_age < 50 ~ "A",
      enrollment_age <= 54 ~ "B",
      enrollment_age <= 59 ~ "C",
      enrollment_age <= 64 ~ "D",
      TRUE ~ "E"),
    sex = if_else(sex == "M", 1, 0)) %>%
  dplyr::select(-enrollment_age) %>%
  dplyr::inner_join(pcas, by = "iid")

all_data <- covar %>%
  dplyr::select(-platform) %>%
  dplyr::inner_join(pheno, by = "iid") %>%
  dplyr::select(iid, vit_25_ohd, tidyselect::everything())

vec_to_tibble <- function(vec, name) {

  namevar <- rlang::sym(name)
  tibble::tibble(iid = names(vec), x = vec) %>%
    dplyr::rename(!! namevar := x)

}

geno_tibb <- purrr::map2(snps$geno, snps$id2, vec_to_tibble) %>%
  purrr::reduce(dplyr::inner_join, by = "iid") %>%
  dplyr::rename(
    rs7041 = `4:72618334:A:C`,
    rs842998 = `4:72627797:G:C`,
    rs842873 = `4:72696465:C:G`,
    rs11731496 = `4:72819988:G:T`,
    rs819988 = `4:72608115:C:T`)

obs_ct <- 5116


data2 <- dplyr::inner_join(all_data, geno_tibb, by = "iid") %>%
  dplyr::select(-iid)

data1 <- data2 %>%
  dplyr::select(-rs819988)

model1 <- lm(vit_25_ohd ~ ., data = data1) %>%
  broom::tidy() %>%
  dplyr::filter(str_detect(term, "rs"))

model2 <- lm(vit_25_ohd ~ ., data = data2) %>%
  broom::tidy() %>%
  dplyr::filter(str_detect(term, "rs"))

snps %<>%
  dplyr::select(id2, data) %>%
  tidyr::unnest(cols = c(data))

all_stats <- here::here("results", folder, "gwas/stats/SCCS_chr4.acount") %>%
    vroom::vroom() %>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::mutate(maf = alt_cts / obs_ct)

stats <- purrr::map(table_snp_keys,
  ~ dplyr::filter(all_stats, stringr::str_detect(id, .))) %>%
  dplyr::bind_rows()

snp_keys <- tibble::tribble(
  ~ rsid, ~ id,
  "rs7041", "4:72618334:A:C",
  "rs842998", "4:72627797:G:C",
  "rs842873", "4:72696465:C:G",
  "rs11731496", "4:72819988:G:T",
  "rs819988", "4:72608115:C:T")

# part1 
snps %>%
  dplyr::inner_join(snp_keys, by = "id") %>%
  dplyr::inner_join(stats, by = "id", suffix = c("", ".x")) %>%
  dplyr::select(-ends_with("x")) %>%
  dplyr::select(-test, -errcode, -id2, -chrom, -pos) %>%
  dplyr::select(-ref, -alt, -a_1, -alt_cts) %>%
  dplyr::select(rsid, everything()) %>%
  dplyr::mutate(
    pve = pve(beta, se, maf, obs_ct),
    pve = scales::scientific(pve)) %>%
  readr::write_tsv(
    here::here("results", "2022_11_21_extra_tables",
      "tsv", "individual_snps.tsv"))

model1 %>%
  dplyr::rename(rsid = term) %>%
  dplyr::inner_join(snp_keys, by = "rsid") %>%
  dplyr::inner_join(snp_keys, by = "id", suffix = c("", ".x")) %>%
  dplyr::inner_join(stats, by = "id", suffix = c("", ".x")) %>%
  dplyr::select(-ends_with("x"), -chrom, -ref, -alt, -alt_cts) %>%
  dplyr::select(rsid, id, everything()) %>%
  dplyr::mutate(
    obs_ct = obs_ct / 2,
    pve = pve(estimate, std.error, maf, obs_ct),
    pve = scales::scientific(pve)) %>%
  readr::write_tsv(
    here::here("results", "2022_11_21_extra_tables",
      "tsv", "joint_model_our_snps.tsv"))

model2 %>%
  dplyr::rename(rsid = term) %>%
  dplyr::inner_join(snp_keys, by = "rsid") %>%
  dplyr::inner_join(snp_keys, by = "id", suffix = c("", ".x")) %>%
  dplyr::inner_join(stats, by = "id", suffix = c("", ".x")) %>%
  dplyr::select(-ends_with("x"), -chrom, -ref, -alt, -alt_cts) %>%
  dplyr::select(rsid, id, everything()) %>%
  dplyr::mutate(
    obs_ct = obs_ct / 2,
    pve = pve(estimate, std.error, maf, obs_ct),
    pve = scales::scientific(pve)) %>%
  readr::write_tsv(
    here::here("results", "2022_11_21_extra_tables",
      "tsv", "joint_model_our_snps_wmoy.tsv"))

