
pacman::p_load(
  magrittr,
  tidyverse,
  ukbtools)

colo_cohort <- qs::qread(
  here::here("results", "2020_06_08_gwas_analysis",
    "qs", "ukbiobank_colorectal_cancer_2020_03_02.qs"))

principal_components <- qs::qread(
  here::here("results", "2020_05_29_principal_component_analysis", "qs",
    "pca_results_plink2.qs")) %>%
  pluck("vectors") %>%
  rename(eid = id)

colo_cohort %<>%
  dplyr::mutate(eid = as.character(eid), fid = eid) %>%
  dplyr::left_join(principal_components, by = "eid") %>%
  dplyr::select(eid, fid, everything())

out_dir <- here::here("results", "2020_06_08_gwas_analysis", "plink2")

colo_cohort %<>%
  dplyr::mutate_if(is.factor,
    list( ~ fct_explicit_na(., "NONE"))) %>%
  dplyr::mutate_if(is.factor,
    list( ~ fct_recode(., NONE = "missing"))) %>%
  dplyr::mutate_if(is.character,
    list( ~ fct_explicit_na(., "NONE")))

save_vars <- function(data, file, sep = " ") {

  vars <- data %>%
    dplyr::select(-eid, -fid) %>%
    names()

  readr::write_lines(
    stringr::str_c(c("#FID", "IID", vars), collapse = sep), file)

  readr::write_delim(data, file, append = TRUE,
    col_names = FALSE, delim = sep)

}

# change levels
colo_cohort %<>%
  dplyr::mutate(
    colorectal = if_else(colorectal == "yes", 2, 1),
    age_bin = fct_recode(age_bin,
      A = "< 50",
      B = "(50,54]",
      C = "(54,59]",
      D = "(59,64]",
      E = ">= 65"),
    household_income = fct_recode(household_income,
      A = " < 18K",
      B = "18K - 31K",
      C = "31K - 52K",
      D = "52K - 100K",
      E = " > 100K"),
    ethnicity = fct_recode(ethnicity,
      A = "african",
      B = "any other asian background",
      B = "any other black background",
      B = "any other mixed background",
      B = "any other white background",
      B = "other",
      C = "bangladeshi",
      D = "british",
      E = "caribbean",
      F = "chinese",
      G = "indian",
      H = "irish",
      I = "pakistani",
      J = "white",
      K = "white and asian",
      L = "white and black african",
      M = "white and black caribbean"))

samples <- colo_cohort %>%
  dplyr::select(fid, eid) %>%
  save_vars(file.path(out_dir, "gwas_samples.txt"))

colorectal_pheno <- colo_cohort %>%
  dplyr::select(fid, eid, colorectal) %>%
  save_vars(file.path(out_dir, "gwas_pheno_colorectal.txt"))

vitd_pheno <- colo_cohort %>%
  dplyr::select(fid, eid, vitamin_d) %>%
  save_vars(file.path(out_dir, "gwas_pheno_vitd.txt"))

sqrt_vitd_pheno <- colo_cohort %>%
  dplyr::mutate(sqrt_vitd = sqrt(vitamin_d)) %>%
  dplyr::select(fid, eid, sqrt_vitd) %>%
  save_vars(file.path(out_dir, "gwas_pheno_sqrt_vitd.txt"))

log_vitd_pheno <- colo_cohort %>%
  dplyr::mutate(log_vitd = log10(vitamin_d)) %>%
  dplyr::select(fid, eid, log_vitd) %>%
  save_vars(file.path(out_dir, "gwas_pheno_log_vitd.txt"))


covs <- colo_cohort %>%
  select(fid, eid, age_bin, bmi_obese, educ_agg, household_income,
    family_history, smoking_status, alcohol_intake, contains("phys_act")) %>%
  save_vars(file.path(out_dir, "gwas_cov_no_pca.txt"))

covs <- colo_cohort %>%
  mutate_at(vars(starts_with("pc")), list(as.character)) %>%
  select(fid, eid, starts_with("pc"), age_bin, bmi_obese,
    educ_agg, household_income, family_history, smoking_status,
    alcohol_intake, contains("phys_act")) %>%
  save_vars(file.path(out_dir, "gwas_cov_all_pca.txt"))

covs <- colo_cohort %>%
  mutate_at(vars(starts_with("pc")), list( ~ . / sd(.))) %>%
  mutate_at(vars(starts_with("pc")), list(as.character)) %>%
  select(fid, eid, starts_with("pc"), age_bin, bmi_obese,
    educ_agg, household_income, family_history, smoking_status,
    alcohol_intake, contains("phys_act")) %>%
  save_vars(file.path(out_dir, "gwas_cov_all_pca_std.txt"))

covs <- colo_cohort %>%
  mutate_at(vars(starts_with("pc")), list( ~ . / sd(.))) %>%
  mutate_at(vars(starts_with("pc")), list(as.character)) %>%
  select(fid, eid, starts_with("pc"), age_bin, bmi_obese,
    educ_agg, household_income, family_history, smoking_status,
    alcohol_intake, contains("phys_act")) %>%
  save_vars(file.path(out_dir, "gwas_cov_all_pca_std.txt"))

covs <- colo_cohort %>%
  select(fid, eid, starts_with("pc"), age, bmi) %>%
  mutate(across(starts_with("pc"),
    list( ~ . / sd(.)), .names = "{.col}")) %>%
  mutate_at(vars(starts_with("pc")), list(as.character))
  
covs %>% ggplot(aes(bmi)) + geom_histogram()
covs %>% ggplot(aes(age)) + geom_histogram()
dev.off()

covs %>%
  save_vars(file.path(out_dir, "gwas_age_bmi_pca_std.txt"))

covs <- colo_cohort %>%
  select(fid, eid, starts_with("pc"), age, bmi) %>%
  mutate(across(starts_with("pc"),
    list( ~ . / sd(.)), .names = "{.col}")) %>%
  mutate(across(one_of("age", "bmi"),
    list( ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)),
      .names = "{.col}")) %>%
  mutate_at(vars(starts_with("pc")), list(as.character))
  
covs %>%
  save_vars(file.path(out_dir, "gwas_center_age_bmi_pca_std.txt"))

# NOTE: removed the ethnicity variable due to the 
# warnings:
# Warning: Skipping --glm regression on phenotype 'colorectal' since variance
# inflation factor for covariate 'ethnicity=D' is too high. You may want to
# remove redundant covariates and try again.
# Warning: Skipping --glm regression on phenotype 'sqrt_vitd' since variance
# inflation factor for covariate 'ethnicity=D' is too high. You may want to
# remove redundant covariates and try again.
