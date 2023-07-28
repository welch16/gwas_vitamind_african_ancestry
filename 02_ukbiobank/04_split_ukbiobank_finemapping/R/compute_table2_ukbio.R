
library(magrittr)
library(tidyverse)
library(gwasrapidd)

source(here::here("src/gwas_funs_redo.R"))

folder <- "2022_06_29_prepare_imputation"
workdir <- here::here("results", folder, "qs", "gwas_clean")
all_files <- list.files(workdir, full.names = TRUE, pattern = "gwas_stats")

stats1 <- all_files %>%
  stringr::str_subset("covar_bmi_age.vdbp") %>%
  purrr::map(qs::qread) %>%
  bind_rows()

# all of our significant variants are from chrom 4, so
# we are going to only read the data from those 

var_stats <- here::here("results", folder, "gwas/stats/SCCS_chr4.acount")

var_stats %<>%
  vroom::vroom() %>%
  dplyr::rename_with(snakecase::to_snake_case)

finemap_snps <- here::here("results", "2022_08_03_finemapping", "qs",
      "sccs_vdbp_fwd_step.qs") %>%
  qs::qread() %>%
  dplyr::rename(id = term)

var_stats %<>%
  dplyr::right_join(dplyr::select(finemap_snps, "id"), by = "id") %>%
  dplyr::mutate(maf = alt_cts / obs_ct)

finemap_snps %<>%
  dplyr::inner_join(
    dplyr::select(var_stats, id, maf, obs_ct), by = "id") %>%
  dplyr::mutate(
    pve = pve(estimate, std.error, maf, obs_ct))

# get covariates of the model
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
  dplyr::select(iid, vdbp)

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
  dplyr::select(iid, vdbp, tidyselect::everything())

# get genotype for the finemapped snps

snp_range <- finemap_snps %>%
  dplyr::pull(id) %>%
  stringr::str_split("\\:") %>%
  purrr::map_dbl(~ as.numeric(.[2])) %>%
  range()

library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)
library(Matrix)


gr <- GenomicRanges::GRanges(
  seqnames = 4,
  ranges = IRanges(start = snp_range[1], end = snp_range[2]))

vcf_files <- here::here("results", "2022_06_29_prepare_imputation",
  "impute_subset") %>%
  list.files(full.names = TRUE, pattern = ".dose.vcf.gz") %>%
  stringr::str_subset(".csi", negate = TRUE) %>%
  stringr::str_subset(".tbi", negate = TRUE)

gmat <- get_genotype(gr, vcf_files, finemap_snps)

gmat %<>%
  as.data.frame() %>%
  tibble::rownames_to_column("iid") %>%
  tibble::as_tibble()

all_data %<>%
  dplyr::inner_join(gmat, by = "iid") %>%
  dplyr::select(-iid)

joint_model <- lm(vdbp ~ ., data = all_data) %>%
  broom::tidy() %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::filter(! term %in% names(all_data)) %>%
  dplyr::filter(!stringr::str_detect(term, "age")) %>%
  dplyr::rename(id = term) %>%
  dplyr::mutate(id = stringr::str_remove_all(id, "`")) %>%
  dplyr::inner_join(
    finemap_snps %>%
      dplyr::select(id, rsids, maf, obs_ct), by = "id") %>%
  dplyr::mutate(
    pve = pve(estimate, std.error, maf, obs_ct))

joint_model %>%
  qs::qsave(
    here::here("results", "2022_08_03_finemapping", "qs",
      "sccs_vdbp_joint_model.qs"))
