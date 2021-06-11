# Program: 03_get_snps_to_flip_list.R
# Author: Rene Welch
# Created: 2019-08-30
# Updated: 2019-08-30
# Purpose: Creates txt list of the snps to filter

library(magrittr)
library(tidyverse)

dr <- here::here("results/2019_08_30_filter_file_correction/rds")
outdr <- here::here("extracted_data", "genotype")

plink <- readRDS(file.path(dr, "2019_08_30_plink_summary_stats_filtered.rds"))

dr0 <- here::here("results/2019_04_22_plink_files_review/rds")

plink_map <- readRDS(file.path(dr0, "2019_04_24_plink_map_clean_with.rds"))

plink_map %<>%
  inner_join(plink, by = c(snp.name = "snp"), suffix = c(".parse", ".plink"))
plink_map %<>% filter(!is.na(maf))

nrow(plink_map) / nrow(plink) ### line 99.93%

plink_map %>% select(a1) %>% table()
plink_map %>% select(a2) %>% table()

# > plink_map %>% select(a1) %>% table()
# .
#      A      C      D      G      I      T
# 182178 157324    876 156317    499 183408
# > plink_map %>% select(a2) %>% table()
# .
#      A      C      D      G      I      T
# 155942 183611    499 182860    876 156814

plink_map %<>%
  mutate(
    base_hg19_neg = case_when(
      base_hg19 == "A" ~ "T",
      base_hg19 == "C" ~ "G",
      base_hg19 == "G" ~ "C",
      base_hg19 == "T" ~ "A",
      TRUE ~ NA_character_))

plink_map %<>%
  mutate(
    allele_sce = case_when(
      base_hg19 == a2 ~ "major_pos",
      base_hg19_neg == a2 ~ "major_neg",
      base_hg19 == a1 ~ "minor_pos",
      base_hg19_neg == a1 ~ "minor_neg",
      TRUE ~ "to_check"))

## prior analysis showed that we can fix some of these variants

plink_map_fix <- plink_map %>%
  filter(allele_sce == "to_check") %>%
  nest(-group) %>%
  mutate(
    n = map_dbl(data, nrow)) %>%
  filter(n  > 30) %>%
  arrange(-n) %>%
  select(-n)

to_fix <- list()

#### pos_based
gr <- "pos_based"
to_fix[[gr]] <- plink_map_fix %>%
  filter(group == gr) %>%
    unnest()

## the pattern for these snps was {chr}:{pos}-{ref}-{alt}
to_fix[[gr]] %<>%
  mutate(
    dashes = str_match_all(snp.name, "-"),
    ndashes = map_int(dashes, nrow),
    alleles = str_split(snp.name, "-") %>% map(~.[-1]),
    new_a1 = map_chr( alleles, 1),
    new_a2 = map_chr(alleles, ~ .[length(.)]),
    a1 = if_else(ndashes == 1, a1, new_a1),
    a2 = if_else(ndashes == 1, a2, new_a2)) %>%
  select(-contains("dashes"), -alleles, -new_a1, -new_a2)

plink_map_fix <- bind_rows(to_fix)
plink_map %<>% filter(allele_sce != "to_check")

### quick check
plink_map %>% select(allele_sce) %>% table()
plink_map %>% filter(allele_sce == "major_neg") %>% select(a2) %>% table()
plink_map %>% filter(allele_sce == "minor_neg") %>% select(a2) %>% table()

rev_comp <- function(nuc) {
  case_when(
    nuc == "A" ~ "T",
    nuc == "T" ~ "A",
    nuc == "C" ~ "G",
    nuc == "G" ~ "C")
}

plink_map %<>%
  mutate(
    a1 = if_else(str_detect(allele_sce, "neg"), rev_comp(a1), a1),
    a2 = if_else(str_detect(allele_sce, "neg"), rev_comp(a2), a2)) %>%
  bind_rows(plink_map_fix)
	
plink_map %>%
  filter(chr.parse != chr.plink) %>%
  filter(chr.plink != "0") %>%
  select(-nchrobs) %>%
  select(contains("chr")) %>%
  table()

clean_snps <- plink_map  %>%
  select(snp.name, chr.parse, pos, base_hg19, a1, a2, maf, p, f_miss) %>%
  rename(
    chr = chr.parse,
    minor = a1,
    major = a2,
    hw.pvalue = p)

clean_snps %>%
  saveRDS(file.path(dr, "2019_08_30_filtered_variants_with_summary_stats.rds"))

plink_map %>%
  filter(str_detect(allele_sce, "neg")) %>%
  select(snp.name, chr.parse, pos, base_hg19, a1, a2, maf, p, f_miss) %>%
  rename(
    chr = chr.parse,
    minor = a1,
    major = a2,
    hw.pvalue = p) %>%
  saveRDS(
    file.path(dr, "2019_08_30_filtered_variants_to_flip_w_summary_stats.rds"))
