pacman::p_load(magrittr, tidyverse, gwasrapidd, snpStats, broom, haploR)

vitd_assoc <-
  qs::qread(
    here::here("results", "zz_methods", "qs", "vitd_associations.qs"))

gc_snps <- vitd_assoc %>%
  purrr::map2(names(vitd_assoc), ~ mutate(.x, study = .y)) %>%
  purrr::map(filter, gene == "GC") %>%
  purrr::map(mutate, chrom = as.character(chrom)) %>%
  dplyr::bind_rows()


ld_query <- haploR::queryHaploreg(
  query = gc_snps$rsid,
  ldPop = "AFR",
  ldThresh = .6)

gc_snps %>%
  dplyr::inner_join(
    ld_query, by = c(rsid = "query_snp_rsid")) %>%
  dplyr::select(study, contains("rs"), r2) %>%
  dplyr::rename(ld_id = rsID) %>%
  dplyr::filter(rsid != ld_id)

ld_query <- haploR::queryHaploreg(
  query = gc_snps$rsid,
  ldPop = "ASN",
  ldThresh = .6)

gc_snps %>%
  dplyr::inner_join(
    ld_query, by = c(rsid = "query_snp_rsid")) %>%
  dplyr::select(study, contains("rs"), r2) %>%
  unique() %>%
  dplyr::rename(ld_id = rsID) %>%
  dplyr::filter(rsid != ld_id)



# A tibble: 8 x 4
#   study           rsid       ld_id       r2
#   <chr>           <chr>      <chr>    <dbl>
# 1 sccs_vitd       rs13142062 rs842999  0.64
# 2 sccs_vitd       rs13142062 rs842999  0.64
# 3 sccs_vitd_joint rs13142062 rs842999  0.64
# 4 sccs_vitd_joint rs13142062 rs842999  0.64
# 5 moy_vitd        rs7041     rs222035  0.84
# 6 moy_vitd        rs7041     rs222047  0.66
# 7 moy_vitd        rs7041     rs705119  0.7 
# 8 moy_vitd        rs7041     rs842999  0.68
