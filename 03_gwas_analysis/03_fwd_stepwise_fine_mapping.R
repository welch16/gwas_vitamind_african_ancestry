# Program: 03_fwd_stepwise_fine_mapping_maf_corrected.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-22
# Purpose: Performs fine-mapping via forward stepwise selection

library(magrittr)
library(tidyverse)
library(snpStats)
library(vroom)

snp_dist <- 250e3
pval_thr <- 5e-8

future::plan(future::multiprocess, workers = 16)

sccs_vdbp_loci <-
  qs::qread(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    str_c("sccs_vitdbp_loci_snpdist",
      scales::comma(snp_dist, scale = 1e-3, suffix = "kp"), "_pval_thr",
      scales::scientific(pval_thr),
    ".qs")))

gwas_dir <- here::here("extracted_data", "SCCS", "imputation",
  "merge_impute", "qc_filtered")

files <- list.files(gwas_dir, pattern = "rsq50", full.names = TRUE)

get_plink_files <- function(all_files, chrom) {

  dir <- dirname(all_files) %>%
    unique()
  all_files <- basename(all_files) %>%
    str_subset(str_c("chr", chrom, "_"))
  file.path(dir, all_files)
}

# this gets genotype for every loci
sccs_vdbp_loci %<>%
  mutate(
    plink_files = map(chrom, ~ get_plink_files(files, .)),
    snpmat = map2(plink_files, locus,
     ~ snpStats::read.plink(
        bed = str_subset(.x, ".bed"),
        bim = str_subset(.x, ".bim"),
        fam = str_subset(.x, ".fam"),
        select.snps = .y$id)))

# need to get response, pca, and covariates
workdir <- here::here("results/2020_09_01_merge_impute_sccs/plink2")
files <- list.files(workdir, full.names = TRUE, pattern = "vitdbp") %>%
  str_subset("rsq50")

build_covariate <- function(plink_files, sccs_files) {

  fam <- plink_files %>%
    str_subset(".fam") %>%
    readr::read_tsv(col_names = FALSE)

  pheno <- sccs_files %>%
    str_subset("_log_pheno") %>%
    read_delim(" ") %>%
    select(-ends_with("FID"))

  covs <- sccs_files %>%
    str_subset("bmi_age") %>%
    str_subset("log", negate = TRUE) %>%
    read_delim(" ") %>%
    select(-ends_with("FID"))

  fam %>%
    rename(IID = X2, sex = X5) %>%
    select(IID, sex) %>%
    inner_join(pheno, by = "IID") %>%
    inner_join(covs, by = "IID") %>%
    select(IID, vitd, everything()) %>%
    mutate(sex = if_else(sex == 1, "m", "f"))

}

build_data <- function(snpmat, data_f) {

  genomat <- snpmat$genotype %>%
    as("numeric")

  genomat %<>%
    as.data.frame() %>%
    as_tibble(rownames = "IID")

  modeldata <- data_f %>%
    inner_join(genomat, by = "IID") %>%
    select(-IID) %>%
    select(vitd, everything()) %>%
    mutate(sex = if_else(sex == "m", 1, 0)) %>%
    na.omit()

  return(modeldata)
}

sccs_vdbp_loci %<>%
  mutate(
    sccs_files = map(chrom, ~ get_plink_files(files, .)),
    data_f = map2(plink_files, sccs_files, build_covariate),
    data_f = map2(snpmat, data_f, build_data))

fwd_step_regression <- function(modeldata, keep_vars) {

  `%<>%` <- magrittr::`%<>%`

  modeldata %<>%
      rename_all(list(~ str_replace_all(., ":", "\\_")))
  fwd_vars <- keep_vars
  stop_criteria <- FALSE
  all_vars <- modeldata %>%
    select(-vitd) %>%
    names()

  linear_model <- function(var, fwd_vars, modeldata) {

    my_data <- modeldata %>%
      select(one_of(c("vitd", fwd_vars, var)))
    model <- lm(vitd ~ ., data = my_data)
    broom::tidy(model) %>%
      filter(term == glue::glue("`{var}`"))
  }

  i <- 1
  while (length(fwd_vars) < length(all_vars) & !stop_criteria) {
    message("iteration:", i)
    rem_vars <- all_vars[!all_vars %in% fwd_vars]
    message("linear models")
    pvalues <- furrr::future_map(rem_vars, linear_model, fwd_vars, modeldata)
    pvalues %<>% dplyr::bind_rows()

    # select variable
    lowest <- pvalues %>%
      top_n(1, -p.value) %>%
      pull(term) %>%
      str_remove_all("`")
    pvalues %<>%
      filter(term != glue::glue("`{x}`", x = lowest))

    stop_criteria <- any(pvalues$p.value <= 1e-8)
    sum(pvalues$p.value <= 1e-8)
    fwd_vars <- c(fwd_vars, lowest)
    message(lowest)
    i <- i + 1
  }

  finaldata <- modeldata %>%
    select(one_of(c("vitd", fwd_vars)))
  finalmodel <- lm(vitd ~ ., data = finaldata)
  out <- broom::tidy(finalmodel)
  out %<>%
    filter(!term %in% keep_vars) %>%
    filter(term != "(Intercept)")
  message("--done")
  return(out)

}

keepvars <- c("sex", "age", "bmi")
keepvars <- c(keepvars, str_c("pc", seq_len(20), sep = "_"))

sccs_vdbp_loci %<>%
  mutate(
    fwd_pvals = map(data_f, safely(fwd_step_regression), keepvars))

sccs_vdbp_loci %>%
  select(locus, fwd_pvals) %>%
  qs::qsave(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    "fwd_step_reg_fine_mapping.qs"))
