# Program: 07_clean_data_for_gwas.R
# Author: Rene Welch
# Created: 2020-09-01
# Updated: 2020-10-06
# Purpose: Create clean adjustment variables files for GWAS analysis

library(magrittr)
library(tidyverse)
library(ukbtools)

data_dir <- here::here("extracted_data", "SCCS")

# function to save stuff in plink2 format
save_vars <- function(data, file, sep = " ") {

  vars <- data %>%
    dplyr::select(-iid, -fid) %>%
    names()

  readr::write_lines(
    stringr::str_c(c("#FID", "IID", vars), collapse = sep), file)

  readr::write_delim(data, file, append = TRUE,
    col_names = FALSE, delim = sep)

}

# read and clean pca
pcadir <- here::here("results", "2020_09_01_merge_impute_sccs", "pca")

pca <- tibble(file = list.files(pcadir)) %>%
  filter(!str_detect(file, ".log")) %>%
  mutate(
    what = if_else(str_detect(file, "vec"), "components", "values"),
    rsq = if_else(str_detect(file, "rsq50"), "rsq50", "rsq30")) %>%
  pivot_wider(names_from = what, values_from = file) %>%
  mutate(
    components = map(file.path(pcadir, components), vroom::vroom),
    values = map(file.path(pcadir, values), read_csv, col_names = FALSE),
    components = map(components, rename_all, list(snakecase::to_snake_case)))

# covariates #################################################################

## followed this questionaire
## https://www.southerncommunitystudy.org/uploads/5/2/7/5/52750661/sccs-questionnaire.pdf
## to clean the data 

vitd_data <- file.path(data_dir, "vitamind_phenotypes.txt") %>%
  readr::read_tsv()

vitd_data %<>%
  dplyr::rename_all(list(snakecase::to_snake_case))


## age
vitd_data %<>%
  dplyr::mutate(
    age = case_when(
      enrollment_age < 50 ~ "A",
      enrollment_age <= 54 ~ "B",
      enrollment_age <= 59 ~ "C",
      enrollment_age <= 64 ~ "D",
      TRUE ~ "E"))

## education
vitd_data %<>%
  dplyr::mutate(
    education = case_when(
      education <= 2 ~ "A", # < high school
      education == 3 ~ "B", # high school,
      TRUE ~ "C")) # > high school %>

## income
vitd_data %<>%
  dplyr::mutate(
    hh_income = case_when(
      hh_income == 9999 ~ NA_character_,
      hh_income <= 1 ~ "A",
      hh_income <= 3 ~ "B",
      TRUE ~ "C"),
    hh_income = forcats::fct_explicit_na(hh_income, "NONE"))

## family history of colorectal cancer
summarize_family_history <- function(
  f_hist, m_hist, b_hist, s_hist) {

  # 0 no
  # 1 yes
  # 7777 - not aplicable
  # 9999 - dont know

  return(
    case_when(
      f_hist == 1 | m_hist == 1 | b_hist == 1 | s_hist == 1 ~ TRUE,
      f_hist %in% c(0, 7777) & m_hist %in% c(0, 7777) & 
      b_hist %in% c(0, 7777) & s_hist %in% c(0, 7777) ~ FALSE,
      TRUE ~ NA))
}

## colorectal and colon cancer family history
vitd_data %<>%
  dplyr::rename(
    si_colorectal_cancer = si_colorectalcancer) %>%
  dplyr::mutate(
    colorectal_fh = summarize_family_history(
      fa_colorectal_cancer,
      mo_colorectal_cancer,
      br_colorectal_cancer,
      si_colorectal_cancer),
    colon_fh = summarize_family_history(
      fa_colon_cancer,
      mo_colon_cancer,
      br_colon_cancer,
      si_colon_cancer),
    family_history = case_when(
      colorectal_fh | colon_fh ~ "yes",
      ! colorectal_fh & ! colon_fh ~ "no",
      TRUE ~ "NONE"))

## bmi obese
vitd_data %<>%
  dplyr::mutate(
    bmi_obese = if_else(bmi >= 30, "yes", "no"),
    bmi_obese = forcats::fct_explicit_na(bmi_obese, "NONE"))

## smoking status and alcohol
vitd_data %<>%
  dplyr::mutate(
    smoking_status = case_when(
      ever_smoker == 0 ~ "never",
      smoke_now == 1 ~ "current",
      smoke_now == 0 ~ "former"),
    alcohol = if_else(alcohol_per_day == 0, "no", "yes"),
    alcohol = forcats::fct_explicit_na(alcohol, "NONE"))

# insurance_coverage
vitd_data %<>%
  dplyr::mutate(
    insurance_coverage = case_when(
      insurance_coverage == 0 ~ "no",
      insurance_coverage == 1 ~ "yes",
      TRUE ~ "NONE"))


# enrollment source appears to be C for everything
vitd_data %<>%
  dplyr::rename(iid = id_number) %>%
  dplyr::mutate(fid = iid)

vitd_data %>% count(enrollment_source)
vitd_data %>% count(enrollment_state)

vitd_data %<>%
  dplyr::select(fid, iid, age, sex,
    education, insurance_coverage,
    hh_income, family_history,
    bmi_obese, smoking_status, alcohol,
    ends_with("met_hr"),
    vit_25_ohd, vdbp) %>%
  dplyr::select(-starts_with("total")) %>%
  dplyr::select(-contains("sport"))

# get samples
sample_data <- tibble(file = file.path(data_dir, "imputation",
  "merge_impute", "qc_filtered") %>%
  list.files(full.names = TRUE, pattern = ".fam", recursive = TRUE)) %>%
  mutate(
    rsq = if_else(str_detect(file, "rsq50"), "rsq50", "rsq30"),
    samples = map(file, vroom::vroom, col_names = FALSE),
    samples = map(samples, dplyr::rename,
      fid = X1,
      iid = X2,
      pat = X3,
      mat = X4,
      sex2 = X5,
      pheno = X6),
    name = basename(file)) %>%
  select(name, everything(), -file)

# want to create the following files
# for both phenotypes vit_25_ohd and vdbp:
# want to generate 3 phenotype files
#   - without transformation, square root and logarithms
 
# principal components need to be scales by its standard deviation

sd_components <- function(pcavec) {

  pcavec %>%
    mutate(across(starts_with("pc"), list(~ . / sd(.)),
      .names = "{.col}"))

}

pca %<>%
  dplyr::mutate(
    components = purrr::map(components, sd_components))

sample_data %<>%
  dplyr::select(name, rsq, samples) %>%
  dplyr::mutate(ids = map(samples, dplyr::select, iid, fid)) %>%
  dplyr::inner_join(pca, by = "rsq")


# add both phenotypes
vit25 <- select(vitd_data, fid, iid, vit_25_ohd) %>%
  rename(vitd = vit_25_ohd) %>%
  mutate(iid = str_c(fid, iid, sep = "_"), fid = 0) %>%
  select(fid, iid, everything())
vitdbp <- select(vitd_data, fid, iid, vdbp) %>%
  rename(vitd = vdbp) %>%
  mutate(iid = str_c(fid, iid, sep = "_"), fid = 0) %>%
  select(fid, iid, everything())


vitd_data %<>%
  mutate(iid = str_c(fid, iid, sep = "_"), fid = 0) %>%
  select(fid, iid, everything())


get_ncols <- function(x) purrr::map_dbl(x, ncol)
get_nrows <- function(x) purrr::map_dbl(x, nrow)

sample_data %<>%
  mutate(
    vit25ohd = map(ids, dplyr::inner_join,
      vit25, by = c("fid", "iid")),
    vitdbp = map(ids, dplyr::inner_join, vitdbp,
      by = c("fid", "iid"))) %>%
  mutate(across(starts_with("vit"), list(~ map(., na.omit)),
    .names = "{.col}"))
    
sample_data %<>%
  select(-samples, -values, -ids)
   
transform_vitd <- function(pheno, fun) {

  purrr::map(pheno, dplyr::mutate, vitd = fun(vitd))

}

sample_data %<>%
  mutate(across(one_of("vit25ohd", "vitdbp"),
    list( ~ transform_vitd(., sqrt)),
    .names = "{.col}_sqrt")) %>%
  mutate(across(one_of("vit25ohd", "vitdbp"),
    list( ~ transform_vitd(., log)),
    .names = "{.col}_log"))

sample_data %<>%
  pivot_longer(starts_with("vit"),
    names_to = "vit",
    values_to = "pheno") %>%
  select(name, vit, everything())

sample_data %<>%
  mutate(
    covar = map(pheno, select, fid, iid) %>%
      map(inner_join,
        select(vitd_data, -sex, -vit_25_ohd, -vdbp),
        by = c("fid", "iid")),
    covar = map2(covar, components, inner_join, by = c("iid"))) %>%
  select(-components)

# save files and queue
out_dir <- here::here("results", "2020_09_01_merge_impute_sccs", "plink2")

sample_data %<>%
  mutate(
    pheno = map(pheno, select, fid, iid, everything()))

sample_data %<>%
  mutate(
    name = stringr::str_remove(name, ".fam"),
    pheno_file = glue::glue("{dir}/{cohort}_{pheno}_pheno.txt",
      dir = out_dir, cohort = name, pheno = vit),
    cov_file = glue::glue("{dir}/{cohort}_{pheno}_covar.txt",
      dir = out_dir, cohort = name, pheno = vit),
    pheno = map2(pheno, pheno_file, save_vars),
    covar = map2(covar, cov_file, save_vars))

queue <- sample_data %>%
  select_if(is.character) %>%
  mutate(
    plink_file = file.path(data_dir,
      "imputation", "merge_impute", "qc_filtered",
      name),
    name = str_c(name, vit, sep = "_"))

queue %<>% select(ends_with("file"), name)
queue %<>% select(plink_file, pheno_file, cov_file, name)
queue %<>% distinct()

queue %<>% mutate_all(list(as.character))

queue %>%
  readr::write_csv(
    here::here("results", "2020_09_01_merge_impute_sccs", "condor", "queues",
    "2020_09_18_sccs_queue.csv"), col_names = FALSE)