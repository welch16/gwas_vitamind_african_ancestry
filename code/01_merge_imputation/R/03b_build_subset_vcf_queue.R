# Program: 03b_build_subset_vcf_queue.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-09-11
# Purpose: Writes queue for subsetting samples with highest rsq scores that
#     were repeated in more than one cohort

library(magrittr)
library(tidyverse)

indir <- here::here("extracted_data", "SCCS", "imputation",
  "imputed_by_chrom")
outdir <- here::here("extracted_data", "SCCS", "imputation",
  "aux_impute")

infiles <- list.files(indir, pattern = "vcf")

samplefiles <- list.files(here::here("results", "2020_09_01_merge_impute_sccs",
  "samples"), full.names = TRUE)

infiles <- tibble::tibble(infiles = file.path(indir, infiles))
samples <- tibble::tibble(samplefiles)

get_cohort <- function(files) {

  case_when(
    str_detect(files, "AA_COPD") ~ "AA_COPD",
    str_detect(files, "BC_CRC") ~ "BC_CRC",
    str_detect(files, "BC_iCOGS") ~ "BC_iCOGS",
    str_detect(files, "BC_ROOT") ~ "BC_ROOT",
    str_detect(files, "CRC_GECCO") ~ "CRC_GECCO",
    str_detect(files, "CRC_Hawaii") ~ "CRC_Hawaii",
    str_detect(files, "Lung_GWAS") ~ "Lung_GWAS",
    str_detect(files, "Lung_exo") ~ "Lung_exo",
    str_detect(files, "pancre") ~ "pancre",
    str_detect(files, "prostate_con") ~ "prostate_con",
    str_detect(files, "Prostate_USC") ~ "Prostate_USC",
    TRUE ~ "vitd")

}

samples %<>% mutate(cohort = get_cohort(samplefiles))
infiles %<>%  mutate(cohort = get_cohort(infiles))

queue <- inner_join(infiles, samples, by = "cohort") %>% na.omit()
queue %<>%
  select(-cohort) %>%
  mutate(
    outfiles = file.path(outdir, basename(infiles)),
    outfiles = str_replace(outfiles, "dose", "sample_filt"),
    outfiles = str_remove(outfiles, ".gz"))

queue %>%
  readr::write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor/queues",
    "2020_09_03_subset_samples.csv"), col_names = FALSE)
