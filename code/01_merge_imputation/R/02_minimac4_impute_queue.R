# Program: 02_minimac4_impute_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-01
# Purpose: Write queue for cohort imputation using minimac4 with condor

library(magrittr)
library(tidyverse)

panel_dir <- here::here("tools/minimac4_panel")
data_dir <- here::here("extracted_data/SCCS/imputation/vcf_by_chrom")
out_dir <- here::here("extracted_data/SCCS/imputation/imputed_by_chrom")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


panel <- tibble(panel = list.files(panel_dir, full.names = TRUE,
  pattern = "m3vcf"))

panel %<>%
  mutate(
    chrom = str_split(basename(panel), "\\."),
    chrom = map_chr(chrom, 1))

vcf <- tibble(vcf = list.files(data_dir, full.names = TRUE, pattern = "vcf"))
vcf %<>%
  mutate(
    chrom = str_split(basename(vcf), "\\_"),
    chrom = map_chr(chrom, ~ .[length(.)]),
    chrom = str_remove(chrom, ".vcf"),
    chrom = str_remove(chrom, "chr"))

queue <- inner_join(vcf, panel, by = "chrom") %>%
  mutate(
    outfile = file.path(out_dir, basename(vcf)),
    outfile = str_remove(outfile, ".vcf"))

queue %>%
  select(vcf, panel, outfile) %>%
  write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor",
      "queues", "2020_09_01_minimac4_impute.csv"),
    col_names = FALSE)