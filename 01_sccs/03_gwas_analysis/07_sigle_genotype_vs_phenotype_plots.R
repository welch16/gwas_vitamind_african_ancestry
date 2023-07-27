# Program: 07_single_genotype_vs_phenotype_plots.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-08
# Purpose: Plots genotype vs phenotype for significant SNPs

library(magrittr)
library(tidyverse)
library(snpStats)
library(ggbio)
library(ggrepel)
library(ggbeeswarm)
library(furrr)
library(future)


# get significant snps
snp_dist <- 250e3
pval_thr <- 5e-8

future::plan(future::multiprocess, workers = 20)

sccs_vdbp_loci <-
  qs::qread(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    str_c("sccs_vitdbp_loci_snpdist",
      scales::comma(snp_dist, scale = 1e-3, suffix = "kp"), "_pval_thr",
      scales::scientific(pval_thr),
    ".qs")))

sccs_vdbp_loci2 <- sccs_vdbp_loci %>%
  tidyr::unnest(cols = c(locus)) %>%
  dplyr::filter(p <= .1) %>%
  dplyr::distinct(chrom, pos, id)

rsids <- here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript",
  "qs", "chrom4_rsids.qs") %>%
  qs::qread()

sccs_vdbp_loci2 %<>%
  dplyr::left_join(rsids, by = c("chrom", "pos", "id"))

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

plink_files <- get_plink_files(files, unique(sccs_vdbp_loci2$chrom))
snp_mat <- snpStats::read.plink(
  bed = stringr::str_subset(plink_files, ".bed"),
  bim = stringr::str_subset(plink_files, ".bim"),
  fam = stringr::str_subset(plink_files, ".fam"),
  select.snps = unique(sccs_vdbp_loci2$id))

workdir <- here::here("results/2020_09_01_merge_impute_sccs/plink2")
files <- list.files(workdir, full.names = TRUE, pattern = "vitdbp") %>%
  stringr::str_subset("rsq50") %>%
  stringr::str_subset("chr4") %>%
  stringr::str_subset("log_pheno")

pheno <- readr::read_delim(files, " ") %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::select(-fid)

plot_genotype_pheno <- function(id, rsid, snp_mat, pheno) {

  geno_mat <- as(snp_mat$genotype, "character") %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "iid") %>%
    dplyr::select(iid, !! rlang::sym(id)) %>%
    dplyr::inner_join(pheno, by = "iid") %>%
    dplyr::rename(geno = !! rlang::sym(id))

  map <- snp_mat$map %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "posid") %>%
    dplyr::filter(posid == id)

  allele1 <- dplyr::pull(map, allele.1)
  allele2 <- dplyr::pull(map, allele.2)

  title <- if_else(is.na(rsid), id, rsid)

  relabel_geno <- function(lvs, allele1, allele2) {

    lvs %>%
      stringr::str_replace_all("A", allele1) %>%
      stringr::str_replace_all("B", allele2)

  }

  geno_mat %<>%
    dplyr::mutate(
      geno = factor(geno, levels = c("A/A", "A/B", "B/B")),
      geno = forcats::fct_relabel(geno, relabel_geno, allele1, allele2))

  geno_mat %>%
    ggplot(aes(geno, vitd)) +
    geom_boxplot() +
    theme_classic() +
    labs(
      x = NULL, y = expression(log(`Vitamin D Binding Protein`)),
      title = title
    )

}

figsdr <- here::here("results", 
  "2021_02_19_update_analysis_vitdbp_manuscript", "figs",
  "snp_boxplot")

fs::dir_create(figsdr)

sccs_vdbp_loci2 %<>%
  dplyr::mutate(
    geno_plot = furrr::future_map2(id, rsid, plot_genotype_pheno, 
      snp_mat, pheno),
    file = glue::glue("{dr}/boxplot_geno_vs_logvitd_{snp}.png",
      dr = figsdr, snp = if_else(is.na(rsid), id, rsid)))

sccs_vdbp_loci2 %>%
  dplyr::mutate(
    purrr::walk2(file, geno_plot, ggplot2::ggsave,
      width = 4.5, height = 3.5, unit = "in"))

