# Program: 10_calculations_regional_assoc_figs.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-04-12
# Purpose: Makes calculations for regional association figures

library(magrittr)
library(tidyverse)
library(snpStats)

# get t-statistics
sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs")) %>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))

vitd_rsids <-
  qs::qread(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs", "chrom4_rsids.qs"))

vitd_rsids %<>%
  dplyr::distinct(id, rsid)

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
    
sccs_vdbp %<>%
  inner_join(vitd_rsids, by = "id")
    
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

compute_ld <- function(locus, snpmat) {

  ld <- cor(as(snpmat$genotype, "numeric"))
  index_snp <- locus %>%
    dplyr::filter(index) %>%
    dplyr::pull(id)

  ld <- ld[, index_snp]
  locus %>%
    dplyr::mutate(ld = ld[locus$id])

}

vitd_assoc <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "vitd_associations.qs"))

sccs_vdbp_loci %<>%
  dplyr::mutate(
    locus = map(locus, ~ mutate(., index = pos %in% vitd_assoc[[1]]$pos))) %>%
  dplyr::filter(map_lgl(locus, ~ any(.$index))) %>%
  dplyr::mutate(
    locus = furrr::future_map2(locus, snpmat, compute_ld))
    
# loci 1
sccs_vdbp %>%
  dplyr::filter(rsid == "rs4282136")
  
id1 <- sccs_vdbp %>%
  dplyr::filter(rsid == "rs4282136") %>%
  pluck("id", 1)
  
id2 <- sccs_vdbp %>%
  dplyr::filter(rsid == "rs1979537") %>%
  pluck("id", 1)
 
id_locus <- sccs_vdbp_loci %>%
  dplyr::filter(map_lgl(locus, ~ any(id1 %in% .$id))) %>%
  dplyr::select(chrom, locus) %>%
  tidyr::unnest(cols = c(locus)) %>%
  dplyr::distinct()
    
id_locus %>% dplyr::filter(ld^2 >= .3) %>% dplyr::filter(!index) %>% nrow()  

# loci  2
sccs_vdbp %>%
  dplyr::filter(rsid == "rs13142062")

id1 <- sccs_vdbp %>%
  dplyr::filter(rsid == "rs13142062") %>%
  pluck("id", 1)

id_locus <- sccs_vdbp_loci %>%
  dplyr::filter(map_lgl(locus, ~ any(id1 %in% .$id))) %>%
  dplyr::select(chrom, locus) %>%
  tidyr::unnest(cols = c(locus)) %>%
  dplyr::distinct()
  
id_locus %>%
  dplyr::filter(ld^2 >= .7) %>%
  inner_join(vitd_rsids)                                                                                                                       
  
  
# other variants 
sccs_vdbp %>%
  dplyr::filter(pos == 72027476)


 sccs_vdbp_loci %>%
  dplyr::filter(map_lgl(locus, ~ any(.$pos == 72027476))) %>%
  dplyr::select(chrom, locus) %>%
  tidyr::unnest(cols = c(locus)) %>%
  dplyr::inner_join(vitd_rsids) %>%
  dplyr::distinct()
  
  
sccs_vdbp_loci %>%
  dplyr::filter(map_lgl(locus, ~ any(.$pos == 72888168))) %>%
  dplyr::select(chrom, locus) %>%
  tidyr::unnest(cols = c(locus)) %>%
  dplyr::inner_join(vitd_rsids) %>%
  dplyr::distinct()
    