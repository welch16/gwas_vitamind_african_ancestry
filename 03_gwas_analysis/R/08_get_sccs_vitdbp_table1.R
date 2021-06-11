# Program: 08_get_sccs_vitdbp_table1.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-05-03
# Purpose: Computes results for table 1

library(magrittr)
library(tidyverse)
library(gwasrapidd)
library(snpStats)
library(broom)
library(haploR)

prev_studies <-
  here::here("results/zz_methods/qs/previous_studies_results.qs") %>%
  qs::qread()
sccs_vitd_annots <-
  here::here("results/zz_methods/qs/sccs_vdbp_significant_with_genes.qs") %>%
  qs::qread()
sccs_vitd_finemap <-
  here::here("results/2021_02_19_update_analysis_vitdbp_manuscript/qs",
    "fwd_step_reg_fine_mapping.qs") %>%
  qs::qread()
 sccs_vitd_colsummary <-
  here::here(
    "results/2021_02_19_update_analysis_vitdbp_manuscript",
    "qs/sccs_vitdbp_col_summary.qs") %>%
  qs::qread()

sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs")) %>%
  dplyr::mutate(
    p = 2 * pnorm(-abs(t_stat)))

sccs_vitd_rsids <-
  here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
      "chrom_1_4_rsids.qs") %>%
  qs::qread() %>%
  dplyr::distinct(id, rsid) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, ":", "\\_"))

biomart_rsids <- qs::qread(
  here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript",
  "qs", "sccs_vitdbp_biomart_queries.qs")) %>%
  dplyr::select(rsid, term) %>%
  dplyr::rename(id = term)

# want tables with columns:
# rsid gene chrom pos ref/alt eaf beta std-err p.value

pve <- function(beta, se, maf, n) {

  num <- 2 * beta^2 * maf * (1 - maf)
  den <- num + 2 * n * se^2 * maf * (1 - maf)
  num / den
}

# piece 1: search the snps comming from fwd-stepwid
sccs_vdbp %<>%
  dplyr::mutate(id = stringr::str_replace_all(id, ":", "\\_"))

finemap_snps <- sccs_vitd_finemap %>%
  dplyr::filter(purrr::map_lgl(fwd_pvals, ~ is.null(.$error))) %>%
  dplyr::mutate(
    fwd_pvals = purrr::map(fwd_pvals, "result")) %>%
  tidyr::unnest(cols = c(fwd_pvals)) %>%
  dplyr::filter(p.value <= 5e-8) %>%
  dplyr::distinct(term) %>%
  dplyr::rename(id = term)

# piece 1: our data
base <- finemap_snps %>%
  dplyr::mutate(id = stringr::str_remove_all(id, "\\`")) %>%
  dplyr::left_join(sccs_vdbp, by = "id")
gene_annots <- sccs_vitd_annots %>%
  dplyr::select(id, genes) %>%
  tidyr::unnest(cols = c(genes)) %>%
  mutate(id = str_replace_all(id, ":", "\\_"))

sccs_vitd_colsummary %<>%
  dplyr::select(id, RAF, MAF) %>%
  mutate(id = str_replace_all(id, ":", "\\_"))

our_table <- base %>%
  left_join(gene_annots, by = "id") %>%
  inner_join(sccs_vitd_colsummary, by = "id") %>%
  inner_join(sccs_vitd_rsids, by = "id") %>%
  dplyr::mutate(
    allele = glue::glue("{ref}/{alt}"),
    allele = as.character(allele)) %>%
  dplyr::select(rsid, symbol, chrom, pos, allele, MAF, beta, se, p) %>%
  dplyr::rename(
    gene = symbol,
    `ref/alt` = allele) %>%
  dplyr::mutate(pve = pve(beta, se, MAF, 2531))


# piece 2: old data
prev_studies %<>%
  filter(author == "moy")

assocs <- prev_studies %>%
  pluck("assoc", 1)
variants <- prev_studies %>%
  pluck("variants", 1)

assocs %>% slotNames()

moy_study <- list(
  assocs %>% pluck("associations"),
  assocs %>% pluck("loci"),
  assocs %>% pluck("risk_alleles"),
  assocs %>% pluck("genes")) %>%
  reduce(partial(inner_join, by = "association_id")) %>%
  inner_join(variants@variants, by = "variant_id")

# rsid gene chrom pos ref/alt eaf beta std-err p.value
moy_study %<>%
  select(variant_id, gene_name, chromosome_name, chromosome_position,
    risk_allele, risk_frequency, beta_number, standard_error, pvalue) %>%
  rename(
    rsid = variant_id, gene = gene_name, chrom = chromosome_name,
    pos = chromosome_position, MAF = risk_frequency, beta = beta_number,
    std.error = standard_error, p.value = pvalue) %>%
  mutate(
    MAF = c(0.35, 0.44, 0.13),
    effect_allele = c("A", "C", "C"),
    risk_allele = str_c(effect_allele, risk_allele, sep = "/")) %>%
  rename(`ref/alt` = risk_allele) %>%
  add_row(
    rsid = "rs668443",
    gene = "GALNT2",
    chrom = "1",
    pos = 228276413,
    `ref/alt` = "A/T",
    MAF = 0.16,
    beta = 457.73,
    std.error = 109.81,
    p.value = 1.36e-5) %>%
  dplyr::mutate(pve = pve(beta, std.error, MAF, 1380))

moy_study_ours <- moy_study %>%
  dplyr::select(rsid)

query <- haploR::queryHaploreg(moy_study_ours$rsid, ldPop = "AFR",
    ldThresh = 0.3) %>%
  dplyr::rename(rsid = rsID)

query %<>%
  dplyr::select(rsid, query_snp_rsid, r2, is_query_snp)

moy_study_ours %<>%
  dplyr::rename(query_snp_rsid = rsid) %>%
  dplyr::inner_join(query, by = "query_snp_rsid") %>%
  dplyr::inner_join(sccs_vitd_rsids, by = "rsid") %>%
  dplyr::left_join(gene_annots, by = "id") %>%
  dplyr::left_join(sccs_vitd_colsummary, by = "id") %>%
  dplyr::inner_join(sccs_vdbp, by = "id") %>%
  dplyr::mutate(
    symbol = if_else(query_snp_rsid == "rs12144344", "ST6GALNAC3", symbol),
    allele = glue::glue("{allele1}/{allele2}",
      allele1 = if_else(ref == a_1, ref, alt),
      allele2 = if_else(ref == a_1, alt, ref)),
    allele = as.character(allele)) %>%
  dplyr::select(query_snp_rsid, rsid, r2, symbol,
    chrom, pos, allele, MAF, beta, se, p) %>%
  dplyr::rename(
    gene = symbol,
    `ref/alt` = allele) %>%
  dplyr::mutate(pve = pve(beta, se, MAF, 1380)) %>%
  dplyr::mutate(gene = if_else(is.na(gene), "GC", gene))

# piece 3: joint model

## 1. generate data_frame vitd, age, bmi, sex, pc_1,...,pc_20, table1 snps

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
datasets <- our_table %>%
  dplyr::bind_rows(moy_study_ours) %>%
  dplyr::select(rsid, chrom) %>%
  dplyr::inner_join(sccs_vitd_rsids, by = "rsid") %>%
  mutate(id = str_replace_all(id, "\\_", ":")) %>%
  dplyr::distinct(chrom, id) %>%
  nest(ids = c(id))

extra_rsids <- sccs_vitd_rsids %>%
  dplyr::filter(rsid == "rs705117" | rsid == "rs13142062") %>%
  dplyr::select(id) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, "\\_", ":"))

datasets$ids[[1]] %<>%
  dplyr::bind_rows(extra_rsids) %>%
  dplyr::distinct()

datasets %<>%
  mutate(
    plink_files = map(chrom, ~ get_plink_files(files, .)),
    snpmat = map2(plink_files, ids,
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

  summary <- snpStats::col.summary(snpmat$genotype)

  genomat <- snpmat$genotype %>%
    as("numeric")

  genomat %<>%
    as.data.frame() %>%
    as_tibble(rownames = "IID")

  to_change <- summary %>%
    dplyr::filter(RAF >= .5) %>%
    rownames()

  genomat %<>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::one_of(to_change), list( ~ 2 - .), .names = "{.col}"))

  modeldata <- data_f %>%
    inner_join(genomat, by = "IID") %>%
    select(-IID) %>%
    select(vitd, everything()) %>%
    mutate(sex = if_else(sex == "m", 1, 0)) %>%
    na.omit()

  return(modeldata)
}

## 2. perform linear models
datasets %<>%
  mutate(
    sccs_files = map(chrom, ~ get_plink_files(files, .)),
    data_f = map2(plink_files, sccs_files, build_covariate),
    data_f = map2(snpmat, data_f, build_data),
    id_summary = map(snpmat, ~ snpStats::col.summary(.$genotype)),
    id_summary = map(id_summary, as_tibble, rownames = "id"))

id_summary <- datasets %>%
  dplyr::select(chrom, id_summary) %>%
  tidyr::unnest(cols = c(id_summary)) %>%
  dplyr::mutate(id = stringr::str_replace_all(id, ":", "\\_"))
  

full_data <- datasets$data_f[[1]] %>%
  dplyr::inner_join(datasets$data_f[[2]],
  by = c("vitd", "sex", "age", "bmi", "pc_1", "pc_2", "pc_3", "pc_4", "pc_5",
    "pc_6", "pc_7", "pc_8", "pc_9", "pc_10", "pc_11", "pc_12",
    "pc_13", "pc_14", "pc_15", "pc_16", "pc_17", "pc_18", "pc_19", "pc_20"))
   
rename_rsids <- function(var, rsids) {

  tibble::tibble(var) %>%
    dplyr::mutate(
      id = stringr::str_replace_all(var, ":", "\\_")) %>%
    dplyr::left_join(rsids, by = "id") %>%
    dplyr::mutate(
      out = if_else(is.na(rsid), var, rsid)) %>%
    dplyr::pull(out)

}

full_data %<>%
  dplyr::rename_with(rename_rsids, rsids = sccs_vitd_rsids)
    
full_data %>%
  qs::qsave(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    "full_data_joint_models.qs"))
    
joint_model_ours <- lm(vitd ~ 1 + sex + age + bmi + pc_1 + pc_2 + pc_3 +
  pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10 + pc_11 + pc_12 +
  pc_13 + pc_14 + pc_15 + pc_16 + pc_17 + pc_18 + pc_19 + pc_20 +
  rs1979537 + rs13142062 + rs10805046 + rs10938008, data = full_data)

joint_model_ours_moy <- lm(vitd ~ 1 + sex + age + bmi + pc_1 + pc_2 + pc_3 +
  pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10 + pc_11 + pc_12 +
  pc_13 + pc_14 + pc_15 + pc_16 + pc_17 + pc_18 + pc_19 + pc_20 +
  rs1979537 + rs13142062 + rs10805046 + rs10938008 +
  rs705119 + rs12239582 + rs705117, data = full_data)

rs13142062_model <- lm(vitd ~ 1 + sex + age + bmi + pc_1 + pc_2 + pc_3 +
  pc_4 + pc_5 + pc_6 + pc_7 + pc_8 + pc_9 + pc_10 + pc_11 + pc_12 +
  pc_13 + pc_14 + pc_15 + pc_16 + pc_17 + pc_18 + pc_19 + pc_20 +
  rs13142062, data = full_data)


clean_model <- function(lm_model) {

  lm_model %>%
    broom::tidy() %>%
    dplyr::rename(rsid = term) %>%
    dplyr::filter(!stringr::str_detect(rsid, "pc")) %>%
    dplyr::filter(! rsid %in% c("sex", "age", "bmi")) %>%
    dplyr::filter(rsid != "(Intercept)")

}

joint_model_ours %<>% clean_model()
joint_model_ours_moy %<>% clean_model()


rs13142062_predictions <-
  list(
    g0 = predict(rs13142062_model,
      newdata = dplyr::mutate(full_data, rs13142062 = 0)),
    g1 = predict(rs13142062_model,
      newdata = dplyr::mutate(full_data, rs13142062 = 1)),
    g2 = predict(rs13142062_model,
      newdata = dplyr::mutate(full_data, rs13142062 = 2)))

rs13142062_model %<>% clean_model()

joint_model_ours %<>%
  dplyr::inner_join(sccs_vitd_rsids, by = "rsid")
joint_model_ours_moy %<>%
  dplyr::inner_join(sccs_vitd_rsids, by = "rsid")
rs13142062_model %<>%
  dplyr::inner_join(sccs_vitd_rsids, by = "rsid")


## 3. format table
joint_model_ours %<>%
  dplyr::left_join(gene_annots, by = "id") %>%
  dplyr::left_join(
  dplyr::select(id_summary, id, RAF, MAF), by = "id")
joint_model_ours_moy %<>%
  dplyr::left_join(gene_annots, by = "id") %>%
  dplyr::left_join(
    dplyr::select(id_summary, id, RAF, MAF), by = "id")
rs13142062_model %<>%
  dplyr::left_join(gene_annots, by = "id") %>%
  dplyr::left_join(
    dplyr::select(id_summary, id, RAF, MAF), by = "id")
  
  
joint_model_ours_moy %<>%
  dplyr::mutate(
    symbol = if_else(
      stringr::str_detect(id, stringr::regex("^1")),
        "ST6GALNAC3", symbol))
joint_model_ours %<>%
  dplyr::inner_join(
    dplyr::select(sccs_vdbp, id, chrom, pos, ref, alt, a_1), by = "id")
joint_model_ours_moy %<>%
  dplyr::inner_join(
    dplyr::select(sccs_vdbp, id, chrom, pos, ref, alt, a_1), by = "id")
rs13142062_model %<>%
  dplyr::inner_join(
    dplyr::select(sccs_vdbp, id, chrom, pos, ref, alt, a_1), by = "id")


joint_model_ours %<>%
  dplyr::rename(gene = symbol, beta = estimate)
joint_model_ours_moy %<>%
  dplyr::rename(gene = symbol, beta = estimate)
rs13142062_model %<>%  
  dplyr::rename(gene = symbol, beta = estimate)
  
joint_model_ours %<>%
  dplyr::mutate(
    `ref/alt` = as.character(glue::glue("{ref}/{alt}")),
    pve = pve(beta, std.error, MAF, 2531))
joint_model_ours_moy %<>%
  dplyr::mutate(
    `ref/alt` = as.character(glue::glue("{ref}/{alt}")),
    pve = pve(beta, std.error, MAF, 2531))
rs13142062_model %<>%
  dplyr::mutate(
    `ref/alt` = as.character(glue::glue("{ref}/{alt}")),
    pve = pve(beta, std.error, MAF, 2531))

    
joint_model_ours %<>%
  dplyr::mutate(gene = if_else(is.na(gene), "NPFFR2", gene))
joint_model_ours_moy %<>%
  dplyr::mutate(gene = if_else(is.na(gene), "NPFFR2", gene))





  
# save stuff
list(
  "sccs_vitd" = our_table,
  "moy_vitd" = moy_study,
  "moy_vitd_ourdata" = moy_study_ours,
  "sccs_vitd_joint" = joint_model_ours,
  "sccs_vitd_joint_w_moy" = joint_model_ours_moy) %>%
  qs::qsave(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript",
      "qs", "vitd_associations.qs"))


folder <- "2021_02_19_update_analysis_vitdbp_manuscript"
fs::dir_create(here::here("results", folder, "csv"))

our_table %>%
  write_csv(
    here::here("results", folder,
    "csv", "sccs_vitd_associations.csv"))

moy_study %>%
  write_csv(
    here::here("results", folder,
    "csv", "moy_vitd_associations.csv"))

moy_study_ours %>%
  write_csv(
  here::here("results", folder, "csv", "moy_vitd_assoc_sccs.csv"))

joint_model_ours %>%
  write_csv(
    here::here("results", folder, "csv", "sccs_vitd_associations_joint.csv"))

joint_model_ours_moy %>%
  write_csv(
    here::here("results", folder, "csv",
      "sccs_vitd_associations_joint_w_moy.csv"))

rs13142062_predictions %>%
  qs::qsave(here::here("results", folder, "qs", "rs13142062_predictions.qs"))
