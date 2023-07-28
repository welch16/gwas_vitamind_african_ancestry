library(magrittr)
library(tidyverse)
library(ukbtools)
library(future)
library(furrr)

future::plan(future::multiprocess(workers = 12))

uk_data <- ukbtools::ukb_df("ukb40833",
  path = here::here("raw_data", "UKBiobank", "2020_02", "rformat"))

uk_data %<>%
  tibble::as_tibble()

## filter colorectal cancer patients or unclassificable

unique_values <- function(.data, key) {

  count_all <- function(.data) {

    vars <- names(.data)
    purrr::map(vars,  ~ count(.data, !! rlang::sym(.)))

  }

  .data %<>%
    dplyr::select(tidyselect::contains(key))

  vars <- names(.data)
  values <- count_all(.data)
  values %<>%
    purrr::map(na.omit)
  nind <- purrr::map(values, dplyr::pull, n) %>%
    purrr::map_dbl(sum)
  ndiff <- purrr::map_int(values, nrow)
  tibble::tibble(vars, nind, ndiff)

}

# remove samples
to_remove <- here::here("results/2021_07_01_ukbiobank_stratified_gwas",
  "csv/sample_to_remove.csv") %>%
  readr::read_csv(col_names = c("eid")) %>%
  dplyr::mutate(eid = as.integer(eid))

uk_data %<>%
  dplyr::anti_join(to_remove, by = "eid")

n_unique_values <- uk_data %>%
  as.list() %>%
  map_dbl(~ length(unique(.)))
vars_unique_values <- names(n_unique_values[n_unique_values <= 1])

# remove variables with unique values
uk_data %<>%
  dplyr::select(-tidyselect::one_of(vars_unique_values))

uk_data %>%
  qs::qsave(here::here("results", "2021_07_01_ukbiobank_stratified_gwas",
    "qs", "ukbiobank_full_data.qs"))

vars_na <- uk_data %>%
  as.list() %>%
  map_dbl(~ sum(is.na(.)) / length(.))
vars_na <- names(vars_na[vars_na >= 0.5])

uk_data %<>%
  dplyr::select(!tidyselect::one_of(vars_na))

uk_data %<>%
  dplyr::rename(
    sex = sex_f31_0_0,
    genetic_sex = genetic_sex_f22001_0_0) %>%
  dplyr::filter(sex == genetic_sex) %>%
  dplyr::select(-genetic_sex)

## education
uk_data %>% select(contains("6138")) %>% names

uk_data %<>%
  dplyr::rename(educ = qualifications_f6138_0_0) %>%
  dplyr::select(-tidyselect::contains("f6138"))

uk_data %<>%
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
uk_data %>% select(contains("738")) %>% names

uk_data %<>%
  dplyr::rename(household_income =
    average_total_household_income_before_tax_f738_0_0) %>%
  dplyr::select(-tidyselect::contains("f738"))

## family history

### illnesss of father
uk_data %>% select(contains("20107")) %>% names
uk_data %>% count(illnesses_of_father_f20107_0_0)

### illnesss of mother
uk_data %>% select(contains("20110")) %>% names
uk_data %>% count(illnesses_of_mother_f20110_0_0)


### illnesss of sibblings
uk_data %>% select(contains("20111")) %>% names
uk_data %>% count(illnesses_of_siblings_f20111_0_0)


family_history_code <- function(var) {

  dplyr::case_when(
    stringr::str_detect(var, "Bowel cancer") ~ "yes",
    stringr::str_detect(var, "not to answer") ~ "dont answer",
    stringr::str_detect(var, "not know") ~ "dont know",
    TRUE ~ "no")

}

uk_data %<>%
  dplyr::rename(
    father_history = illnesses_of_father_f20107_0_0,
    mother_history = illnesses_of_mother_f20110_0_0,
    sibblings_history = illnesses_of_siblings_f20111_0_0
  ) %>%
  dplyr::mutate_at(dplyr::vars(tidyselect::ends_with("history")),
    list(family_history_code)) %>%
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
uk_data %>% select(contains("20116")) %>% names
uk_data %>% select(contains("20161")) %>% names

uk_data %<>%
  dplyr::rename(
    smoking_status = smoking_status_f20116_0_0) %>%
  dplyr::mutate(
    smoking_status = stringr::str_to_lower(smoking_status)) %>%
  dplyr::select(-tidyselect::contains("20116"))

## alcohol consumption
uk_data %>% unique_values("1558")
uk_data %>% count(alcohol_intake_frequency_f1558_0_0)

uk_data %<>%
  dplyr::rename(alcohol_intake = alcohol_intake_frequency_f1558_0_0) %>%
  dplyr::mutate(alcohol_intake = case_when(
      stringr::str_detect(alcohol_intake, "answer") ~ "dont answer",
      is.na(alcohol_intake) ~ "dont know",
      alcohol_intake == "Never" ~ "no",
      TRUE ~ "yes")) %>%
  dplyr::select(-tidyselect::contains("1558"))

## physical activity: these are categorical values
uk_data %>% unique_values("884")
uk_data %>% unique_values("904")

uk_data %>% select(contains("884")) %>% names
uk_data %>% select(contains("904")) %>% names

uk_data %<>%
  dplyr::rename(
    mod_phys_act =
      number_of_daysweek_of_moderate_physical_activity_10_minutes_f884_0_0,
    vig_phys_act =
      number_of_daysweek_of_vigorous_physical_activity_10_minutes_f904_0_0) %>%
  dplyr::select(
    -tidyselect::contains("f884"),
    -tidyselect::contains("f904"))

## bmi
uk_data %<>%
  dplyr::rename(bmi = body_mass_index_bmi_f21001_0_0) %>%
  dplyr::mutate(bmi_obese =  if_else(bmi >= 30, "yes", "no")) %>%
  dplyr::select(-tidyselect::contains("f21001"))

## vitamin D
uk_data %>% unique_values("f30890")

uk_data %<>%
  dplyr::rename(vitamin_d = vitamin_d_f30890_0_0) %>%
  dplyr::select(-tidyselect::contains("f30890"))

## year of birth, date of death
uk_data %>% unique_values("34_")
uk_data %>% unique_values("21022")
uk_data %>% unique_values("53_0")

uk_data %<>%
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
uk_data %>% unique_values("21000")

uk_data %<>%
  dplyr::rename(ethnicity = ethnic_background_f21000_0_0) %>%
  dplyr::select(
    -tidyselect::contains("f21000"))

## bin the age
uk_data %<>%
  dplyr::mutate(
    age_bin = cut(age, breaks = c(0, 50, 54, 59, 64, Inf)),
    age_bin = forcats::fct_recode(age_bin,
      `< 50` = "(0,50]",
      `>= 65` = "(64,Inf]"),
    age_bin = forcats::fct_reorder(age_bin, age, .fun = min))

## remove remaining variables and save for
colo_vars <- c("eid", "vitamin_d", "age", "age_bin",
  "bmi", "bmi_obese",
  "age_recrui", "sex", "educ", "educ_agg",
  "household_income", "father_history", "mother_history",
  "sibblings_history", "family_history", "smoking_status",
  "alcohol_intake", "mod_phys_act", "vig_phys_act", "year_birth",
  "ethnicity", "date_center")

uk_data %<>% dplyr::select(tidyselect::one_of(colo_vars))

uk_data %<>%
  dplyr::mutate_if(is.ordered, list(as.character)) %>%
  dplyr::mutate_if(is.character, list(stringr::str_to_lower))

uk_data %<>%
  dplyr::mutate(
    pop = case_when(
      ethnicity %in% c("african", "any other black background",
    "caribbean", "black or black british") ~ "afrocaribbean",
      ethnicity %in% c("any other asian background",
    "bangladeshi", "chinese", "indian", "pakistani",
    "asian or asian british") ~ "asian",
      ethnicity %in% c("any other white background", "british",
    "irish", "white") ~ "white",
    TRUE ~ "other"))

uk_data %>%
  qs::qsave(here::here("results", "2021_07_01_ukbiobank_stratified_gwas",
    "qs", "ukbiobank_grant_vars.qs"))

