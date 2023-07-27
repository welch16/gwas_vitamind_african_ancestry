# Program: 09_get_rsids.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-08
# Purpose: Get rsids for SNPs
  
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(UWCCC.GWAS)
library(biomaRt)

all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snp_chr1 <- BSgenome::snpsBySeqname(all_snps, "1")
snp_chr4 <- BSgenome::snpsBySeqname(all_snps, "4")

pacman::p_load(magrittr, tidyverse)

# get significant snps
vitd_assoc <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs"))

# convert snp_chr into tibble
snps1 <- with(snp_chr1,
  tibble::tibble(chrom = 1, pos, rsid = RefSNP_id))
snps4 <- with(snp_chr4,
  tibble::tibble(chrom = 4, pos, rsid = RefSNP_id))

snps <- dplyr::bind_rows(snps1, snps4)

vitd_assoc %<>%
  dplyr::filter(chrom %in% c(1,4)) %>%
  dplyr::distinct(chrom, pos, id, ref, alt, a_1) %>%
  dplyr::left_join(snps, by = c("chrom", "pos"))
  
vitd_assoc %>%
  qs::qsave(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs", "chrom_1_4_rsids.qs"))
