# Program: 04_summarize_results.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-03-08
# Purpose: Summarizes fine-mapping results

library(magrittr)
library(tidyverse)
library(UWCCC.GWAS)

pacman::p_load(magrittr, tidyverse, UWCCC.GWAS)

vitd_rsids <-
  qs::qread(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs", "chrom4_rsids.qs"))

vitd_rsids %<>%
  dplyr::distinct(id, rsid) %>%
  dplyr::mutate(
    id = stringr::str_replace_all(id, ":", "\\_"))

sccs_vdbp_loci <-
  qs::qread(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    "fwd_step_reg_fine_mapping.qs"))

sccs_vdbp_loci %<>%
  dplyr::filter(map_lgl(fwd_pvals, ~ is.null(.$error))) %>%
  dplyr::mutate(
    fwd_pvals = map(fwd_pvals, "result")) %>%
  dplyr::select(-locus) %>%
  tidyr::unnest(cols = c(fwd_pvals))

sccs_vdbp_loci %<>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::arrange(p.value) %>%
  dplyr::distinct()

sccs_vdbp_loci %<>%
  dplyr::mutate(
    term = str_remove_all(term, "`"),
    ex = str_split(term, "\\_"),
    chrom = map_chr(ex, 1),
    pos = map_chr(ex, 2),
    ref = map_chr(ex, 3),
    alt = map_chr(ex, 4)) %>%
  dplyr::mutate(
    dplyr::across(one_of("chrom", "pos"), list(as.numeric),
      .names = "{.col}")) %>%
  dplyr::select(chrom, pos, ref, alt, everything(), -ex)

sccs_vdbp_loci %<>%
  dplyr::left_join(
    vitd_rsids, by = c(term = "id")) %>%
  dplyr::mutate(
    rsid = if_else(is.na(rsid), term, rsid))

sccs_vdbp_loci %>%
  qs::qsave(
    here::here("results", "2021_02_19_update_analysis_vitdbp_manuscript", "qs",
    "sccs_vitdbp_biomart_queries.qs"))

