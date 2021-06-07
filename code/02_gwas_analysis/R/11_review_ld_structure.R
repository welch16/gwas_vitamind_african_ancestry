# Program: 11_review_ld_structure.R
# Author: Rene Welch
# Created: 2020-05-03
# Updated: 2021-05-03
# Purpose: Reviews LD structure questions before submissions

library(magrittr)
library(tidyverse)
library(snpStats)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)

# get significant snps
vitd_assoc <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "vitd_associations.qs"))

vitd_rsids <-
  qs::qread(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs", "chrom4_rsids.qs"))

vitd_rsids %<>%
  dplyr::distinct(id, rsid)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
human <- org.Hs.eg.db
ensdb <- EnsDb.Hsapiens.v75

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
my_locus <- bind_rows(sccs_vdbp_loci$locus)
plink_files <- get_plink_files(files, "4")
snpmat <- snpStats::read.plink(
  bed = str_subset(plink_files, ".bed"),
  bim = str_subset(plink_files, ".bim"),
  fam = str_subset(plink_files, ".fam"),
  select.snps = unique(my_locus$id))

compute_ld <- function(snpmat) {

  ld <- cor(as(snpmat$genotype, "numeric"))
  ld
}

ldmat <- compute_ld(snpmat)


get_id <- function(rsid_q, vitd_dict) {

  vitd_dict %>%
    dplyr::filter(rsid == rsid_q) %>%
    dplyr::pull(id)

}

	

ldmat[get_id("rs13142062", vitd_rsids), get_id("rs705119", vitd_rsids)]

sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs"))

sccs_vdbp %<>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))


sccs_vdbp %>%
  dplyr::filter(id == get_id("rs705119", vitd_rsids))

haploR::queryHaploreg(c("rs705119", "rs13142062"),
  ldPop = "AFR", ldThresh = .3) %>%
  arrange(query_snp_rsid, desc(r2)) %>%
  dplyr::select(rsID, query_snp_rsid, r2)

# A tibble: 18 x 3
   rsID       query_snp_rsid    r2
   <chr>      <chr>          <dbl>
 1 rs13142062 rs13142062      1   
 2 rs222035   rs13142062      0.4 
 3 rs222047   rs13142062      0.56
 4 rs7041     rs13142062      0.44
 5 rs705117   rs13142062      0.32
 6 rs705119   rs13142062      0.53
 7 rs842999   rs13142062      0.64
 8 rs13142062 rs705119        0.53
 9 rs222004   rs705119        0.4 
10 rs222035   rs705119        0.65
11 rs222047   rs705119        0.69
12 rs7041     rs705119        0.7 
13 rs705117   rs705119        0.34
14 rs705118   rs705119        0.33
15 rs705119   rs705119        1   
16 rs71213586 rs705119        0.41
17 rs842877   rs705119        0.34
18 rs842999   rs705119        0.84

r$> haploR::queryHaploreg(c("rs705119", "rs13142062"), 
      ldPop = "AFR", ldThresh = .3) %>% 
      arrange(query_snp_rsid, desc(r2)) %>% 
      dplyr::select(rsID, query_snp_rsid, r2)                                                                     
      
# this returns:   
# A tibble: 18 x 3
#    rsID       query_snp_rsid    r2
#    <chr>      <chr>          <dbl>
#  1 rs13142062 rs13142062      1   
#  2 rs842999   rs13142062      0.64
#  3 rs222047   rs13142062      0.56
#  4 rs705119   rs13142062      0.53
#  5 rs7041     rs13142062      0.44
#  6 rs222035   rs13142062      0.4 
#  7 rs705117   rs13142062      0.32
#  8 rs705119   rs705119        1   
#  9 rs842999   rs705119        0.84
# 10 rs7041     rs705119        0.7 
# 11 rs222047   rs705119        0.69
# 12 rs222035   rs705119        0.65
# 13 rs13142062 rs705119        0.53
# 14 rs71213586 rs705119        0.41
# 15 rs222004   rs705119        0.4 
# 16 rs705117   rs705119        0.34
# 17 rs842877   rs705119        0.34
# 18 rs705118   rs705119        0.33