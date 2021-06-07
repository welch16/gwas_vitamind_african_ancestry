# Program: 05_compute_row_summaries.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-08
# Purpose: Compute SNP summary statistics for table 1

pacman::p_load(magrittr, tidyverse, snpStats, broom)

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

sccs_vdbp_loci %<>%
  mutate(
    id_summary = map(snpmat, ~ snpStats::col.summary(.$genotype)),
    id_summary = map(id_summary, as_tibble, rownames = "id"))

sccs_vdbp_loci %>%
  dplyr::select(id_summary) %>%
  tidyr::unnest(cols = c(id_summary)) %>%
  dplyr::distinct() %>%
  qs::qsave(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    "sccs_vitdbp_col_summary.qs"))
