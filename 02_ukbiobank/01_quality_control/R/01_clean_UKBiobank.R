
pacman::p_load(
  magrittr,
  tidyverse,
  ukbtools,
  future,
  furrr)

future::plan( future::multiprocess( workers = 12))

uk_data <- ukbtools::ukb_df("ukb40833",
  path = here::here("raw_data", "UKBiobank", "2020_02", "rformat"))

data_dir <- here::here("results",
  "2020_02_27_clean_UKbiobank_data_2nd_attempt", "rds")
out_dir <- here::here("results", "2020_03_24_genotype_QC_and_imputation",
  "samples")
keep_file <- file.path(out_dir, "colo_samples_to_keep.samples")

uk_data %<>% tibble::as_tibble()

## filter colorectal cancer patients or unclassificable

unique_values <- function(.data, key) {

  `%<>%` <- magrittr::`%<>%`
  count_all <- function(.data) {

    vars <- names(.data)
    purrr::map(vars,  ~ count(.data, !! rlang::sym(.)))

  }

  .data %<>% dplyr::select(tidyselect::contains(key))

  vars <- names(.data)
  values <- count_all(.data)
  values %<>% purrr::map(na.omit)
  nind <- purrr::map(values, dplyr::pull, n) %>%
    purrr::map_dbl(sum)
  ndiff <- purrr::map_int(values, nrow)
  tibble::tibble(vars, nind, ndiff)

}

nsamples_all <- nrow(uk_data)
icd10chapters %<>% tibble::as_tibble()
icd10codes %<>% tibble::as_tibble()

# two cohorts:
# variable 20001-0.0 asks whether cancer status was self-reported before 
#   the ukbiobank project started
# variable 40006-0.0 indicates the ICD10 code determining possible cancers

uk_data %>% unique_values("20001")
uk_data %>% select(contains("20001"))
uk_data %>% count(cancer_code_selfreported_f20001_0_0)

cases_self_reported <- uk_data %>%
  dplyr::filter(cancer_code_selfreported_f20001_0_0 %in%
    c("1020", "1021", "1022", "1023")) %>%
  dplyr::select(eid, cancer_code_selfreported_f20001_0_0)

cases_self_reported %>% count(cancer_code_selfreported_f20001_0_0)

## colon and rectum codes
rectum_colon_codes <- bind_rows(
  ukb_icd_keyword("colon", 10),
  ukb_icd_keyword("rectum", 10))

cases_icd10 <- uk_data %>%
  dplyr::filter(type_of_cancer_icd10_f40006_0_0 %in%
    pull(rectum_colon_codes, code)) %>%
  dplyr::filter(!str_detect(type_of_cancer_icd10_f40006_0_0, "D")) %>%
  dplyr::filter(type_of_cancer_icd10_f40006_0_0 != "C785") %>%
  dplyr::select(eid, type_of_cancer_icd10_f40006_0_0)



cases_icd10 %>% count(type_of_cancer_icd10_f40006_0_0) %>%
  rename(code = type_of_cancer_icd10_f40006_0_0) %>%
  left_join(rectum_colon_codes)


controls <- uk_data %>%
  dplyr::filter(is.na(type_of_cancer_icd10_f40006_0_0) &
    is.na(cancer_code_selfreported_f20001_0_0)) %>%
  dplyr::select(eid)

colo_cohort <-
  dplyr::bind_rows(
    cases_self_reported %>% select(eid) %>% mutate(colorectal = 1),
    cases_icd10 %>% select(eid) %>% mutate(colorectal = 1),
    controls %>% mutate(colorectal = 0)) %>%
  dplyr::left_join(uk_data, by = "eid") %>%
  dplyr::distinct()


## sex: removed one trans individual
colo_cohort %>% select(starts_with("sex")) %>% names

colo_cohort %<>%
  dplyr::rename(
    sex = sex_f31_0_0,
    genetic_sex = genetic_sex_f22001_0_0) %>%
  dplyr::filter(sex == genetic_sex) %>%
  dplyr::select(-genetic_sex)

## education
colo_cohort %>% select(contains("6138")) %>% names

colo_cohort %<>%
  dplyr::rename(educ = qualifications_f6138_0_0) %>%
  dplyr::select(-tidyselect::contains("f6138"))

colo_cohort %<>%
  dplyr::mutate(
    educ_agg = fct_collapse(educ,
      highschool = c("O levels/GCSEs or equivalent",
          "CSEs or equivalent"),
      college_technical = c("A levels/AS levels or equivalent",
        "NVQ or HND or HNC or equivalent",
        "College or University degree",
        "Other professional qualifications eg: nursing, teaching")),
    educ_agg = fct_explicit_na(educ_agg, "other"))

## household income
colo_cohort %>% select(contains("738")) %>% names

colo_cohort %<>%
  dplyr::rename(household_income = 
    average_total_household_income_before_tax_f738_0_0) %>%
  dplyr::select(-tidyselect::contains("f738"))

## family history

### illnesss of father
colo_cohort %>% select(contains("20107")) %>% names
colo_cohort %>% count(illnesses_of_father_f20107_0_0)

### illnesss of mother
colo_cohort %>% select(contains("20110")) %>% names
colo_cohort %>% count(illnesses_of_mother_f20110_0_0)


### illnesss of sibblings
colo_cohort %>% select(contains("20111")) %>% names
colo_cohort %>% count(illnesses_of_siblings_f20111_0_0)


family_history_code <- function(var) {

  dplyr::case_when(
    stringr::str_detect(var, "Bowel cancer") ~ "yes",
    stringr::str_detect(var, "not to answer") ~ "dont answer",
    stringr::str_detect(var, "not know") ~ "dont know",
    TRUE ~ "no")

}

colo_cohort %<>%
  dplyr::rename(
    father_history = illnesses_of_father_f20107_0_0,
    mother_history = illnesses_of_mother_f20110_0_0,
    sibblings_history = illnesses_of_siblings_f20111_0_0
  ) %>%
  dplyr::mutate_at(dplyr::vars(tidyselect::ends_with("history")),
    list( family_history_code)) %>%
  dplyr::mutate(
    family_history = case_when(
      father_history == "yes" | mother_history == "yes" |
        sibblings_history == "yes" ~ "yes",
      father_history == "dont know" | mother_history == "dont know" |
        sibblings_history == "dont know" ~ "dont know",
      father_history == "dont answer" | mother_history == "dont answer" |
        sibblings_history == "dont answer" ~ "dont answer",
      TRUE ~ "no")) %>%
  dplyr::select(
    -tidyselect::contains("f20107"),
    -tidyselect::contains("f20110"),
    -tidyselect::contains("f20111"))

## smoking status and packs per year
colo_cohort %>% select(contains("20116")) %>% names
colo_cohort %>% select(contains("20161")) %>% names

colo_cohort %<>%
  dplyr::rename(
    smoking_status = smoking_status_f20116_0_0) %>%
  dplyr::mutate(
    smoking_status = stringr::str_to_lower(smoking_status)) %>%
  dplyr::select(-tidyselect::contains("20116")) %>%
  dplyr::select(-tidyselect::contains("20161"))

## alcohol consumption
colo_cohort %>% unique_values("1558")
colo_cohort %>% count(alcohol_intake_frequency_f1558_0_0)

colo_cohort %<>%
  dplyr::rename(alcohol_intake = alcohol_intake_frequency_f1558_0_0) %>%
  dplyr::mutate(alcohol_intake = case_when(
      stringr::str_detect(alcohol_intake, "answer") ~ "dont answer",
      is.na(alcohol_intake) ~ "dont know",
      alcohol_intake == "Never" ~ "no",
      TRUE ~ "yes")) %>%
  dplyr::select(-tidyselect::contains("1558"))

## physical activity: these are categorical values
colo_cohort %>% unique_values("884")
colo_cohort %>% unique_values("904")

colo_cohort %>% select(contains("884")) %>% names
colo_cohort %>% select(contains("904")) %>% names

colo_cohort %<>%
  dplyr::rename(
    mod_phys_act =
      number_of_daysweek_of_moderate_physical_activity_10_minutes_f884_0_0,
    vig_phys_act =
      number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904_0_0) %>%
  dplyr::select(
    -tidyselect::contains("f884"),
    -tidyselect::contains("f904"))

## bmi
colo_cohort %<>%
  dplyr::rename(bmi = body_mass_index_bmi_f21001_0_0) %>%
  dplyr::mutate(bmi_obese =  if_else(bmi >= 30, "yes", "no")) %>%
  dplyr::select(-tidyselect::contains("f21001"))

## vitamin D
colo_cohort %>% unique_values("f30890")

colo_cohort %<>%
  dplyr::rename(vitamin_d = vitamin_d_f30890_0_0) %>%
  dplyr::select(-tidyselect::contains("f30890"))

## year of birth, date of death
colo_cohort %>% unique_values("34_")
colo_cohort %>% unique_values("40000")
colo_cohort %>% unique_values("40005")
colo_cohort %>% unique_values("21022")
colo_cohort %>% unique_values("53_0")

colo_cohort %<>%
  dplyr::rename(
    year_birth = year_of_birth_f34_0_0,
    age_recrui = age_at_recruitment_f21022_0_0,
    date_center = date_of_attending_assessment_centre_f53_0_0) %>%
  dplyr::mutate(
    age = lubridate::year(date_center) - year_birth) %>%
  dplyr::select(
    -tidyselect::contains("f34_"),
    -tidyselect::contains("f21022"),
    -tidyselect::contains("f53_"))

### race
colo_cohort %>% unique_values("21000")

colo_cohort %<>%
  dplyr::rename(ethnicity = ethnic_background_f21000_0_0) %>%
  dplyr::select(
    -tidyselect::contains("f21000"))

## bin the age
colo_cohort %<>%
  dplyr::mutate(
    age_bin = cut(age, breaks = c(0, 50, 54, 59, 64, Inf)),
    age_bin = forcats::fct_recode(age_bin,
      `< 50` = "(0,50]",
      `>= 65` = "(64,Inf]"),
    age_bin = forcats::fct_reorder(age_bin, age, .fun = min))

## remove remaining variables and save for 
colo_vars <- c("eid", "colorectal", "vitamin_d", "age", "age_bin",
  "bmi", "bmi_obese",
  "age_recrui", "sex", "educ", "educ_agg",
  "household_income", "father_history", "mother_history",
  "sibblings_history", "family_history", "smoking_status",
  "alcohol_intake", "mod_phys_act", "vig_phys_act", "year_birth",
  "ethnicity", "date_center")

colo_cohort %<>% dplyr::select(tidyselect::one_of(colo_vars))

colo_cohort %<>%
  dplyr::mutate_if(is.ordered, list(as.character)) %>%
  dplyr::mutate_if(is.character, list(stringr::str_to_lower))

colo_cohort %<>%
  dplyr::mutate(eid = as.character(eid)) %>%
  dplyr::distinct()

colo_cohort %>% nrow()
colo_cohort %>% na.omit() %>% nrow()


colo_cohort %>% as.list() %>%
  map(is.na) %>% map_int(sum)

colo_samples <- here::here(
  "results",
  "2020_03_24_genotype_QC_and_imputation",
  "samples/colo_samples_to_keep.samples") %>%
  readr::read_delim(" ")

samples <- colo_samples %>%
  slice(-1) %>%
  mutate(eid = as.character(ID_1)) %>%
  select(eid)

colo_cohort %<>%
  dplyr::inner_join(samples, by = "eid") %>%
  dplyr::filter(!is.na(vitamin_d))

# modification per variable

colo_cohort %<>%
  dplyr::mutate(
    colorectal = if_else(colorectal == 0, "no", "yes"),
    colorectal = factor(colorectal))


colo_cohort %<>%
  mutate(
    educ_agg = fct_recode(educ_agg,
      other = "none of the above",
      missing = "prefer not to answer"),
    household_income = fct_explicit_na(household_income,
      na_level = "missing"),
    household_income = fct_recode(household_income,
      `18K - 31K` = "18,000 to 30,999",
      `31K - 52K` = "31,000 to 51,999",
      `52K - 100K` = "52,000 to 100,000",
      ` > 100K` = "greater than 100,000",
      ` < 18K` = "less than 18,000",
      missing = "do not know",
      missing = "prefer not to answer"))

colo_cohort %<>%
  mutate_at(vars(ends_with("history")),
    list( ~ fct_recode(.,
      missing = "dont answer",
      missing = "dont know")))

colo_cohort %<>%
  mutate(
    smoking_status = fct_explicit_na(smoking_status, "missing"),
    smoking_status = fct_recode(smoking_status,
      missing = "prefer not to answer"),
    alcohol_intake = fct_recode(alcohol_intake,
      missing = "dont answer", missing = "dont know"))
    
clean_phys <- function(phys) {

  factor(phys) %>%
  fct_explicit_na("missing") %>%
    fct_recode(
      missing = "-3",
      missing = "-1",
      no = "0",
      no = "1",
      no = "2",
      no = "3",
      no = "4",
      yes = "5",
      yes = "6",
      yes = "7")
}

colo_cohort %<>%
  mutate_at(vars(ends_with("phys_act")),
    list(clean_phys))

colo_cohort %<>%
  mutate(
    ethnicity = fct_explicit_na(ethnicity, "missing"),
    ethnicity = fct_recode(ethnicity,
      missing = "do not know",
      missing = "prefer not to answer",
      other = "other ethnic group",
      other = "mixed",
      other = "black or black british",
      other = "asian or asian british"))

colo_cohort %<>%
  dplyr::select(-year_birth, -date_center)

colo_cohort %>%
  qs::qsave(here::here("results", "2020_06_08_gwas_analysis",
    "qs", "ukbiobank_colorectal_cancer_2020_03_02.qs"))


cases_icd10 %>%
  dplyr::mutate(eid = as.character(eid)) %>%
  dplyr::inner_join(samples, by = "eid") %>%
  dplyr::count(type_of_cancer_icd10_f40006_0_0) %>%
  dplyr::rename(code = type_of_cancer_icd10_f40006_0_0) %>%
  dplyr::left_join(rectum_colon_codes) %>%
  qs::qsave(
    here::here("results", "2020_06_08_gwas_analysis",
      "qs", "cases_meaning.qs"))






## define files and directories

## read colorectal cancer data
# colo <- read_rds(file.path(data_dir,
#   "ukbiobank_colorectal_cancer_2020_03_02.rds"))
# colo_samples <- colo %>% pull(eid)

# bgen_sample <- here::here("raw_data/UKBiobank/2020_02/genotype",
#   "ukb35256_imp_v3.sample") %>%
#   read_delim(delim = " ")

# second_header <- bgen_sample %>% slice(1)

# bgen_sample %<>% filter( ID_1 %in% colo_samples)
# bgen_sample <- bind_rows(second_header, bgen_sample)
# bgen_sample %>% write_delim(keep_file, delim = " ", col_names = TRUE)
