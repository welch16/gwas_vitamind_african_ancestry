
pacman::p_load(
  magrittr,
  tidyverse,
  broom,
  future,
  furrr)

colo_cohort <- qs::qread(
  here::here("results",
    "2020_06_08_gwas_analysis",
    "qs",
    "ukbiobank_colorectal_cancer_2020_03_02.qs"))

future::plan(multiprocess(workers = 8))

ind_var_models <- colo_cohort %>%
  dplyr::select(-eid, -colorectal, -age, -age_recrui, -educ) %>%
  names() %>% tibble::tibble(var = .)

build_data <- function(variable, colo ) {

  colo %<>%
    dplyr::select(colorectal, !! rlang::sym(variable)) %>%
    dplyr::rename( variable = !! rlang::sym(variable)) %>%
    dplyr::filter( !is.na(variable))

  if (is.numeric(colo$variable)) {

    colo %<>%
      dplyr::filter( variable > 0) %>%
      dplyr::filter(!is.na(variable))

  }else{

    colo %<>%
      dplyr::filter( !stringr::str_detect(variable, "missing"))

  }

  colo %<>%
    dplyr::mutate(variable = factor(variable))

  return(colo)

}

count_cases <- function(data, var) {

  if (is.numeric(var)) {

    out <- tibble::tibble(
            term = rep("", 0),
            no = rep(1, 0),
            yes = rep(1, 0),
            n = rep(1, 0))
  } else {

    out <- data %>%
      dplyr::count(colorectal, variable) %>%
      tidyr::spread(colorectal, n, fill = 0) %>%
      dplyr::rename(term = variable) %>%
      dplyr::mutate(
        n = yes + no,
        term = stringr::str_c(var, term, sep = ": "))

  }

  return(out)
}

compute_pvalues <- function(model, var) {

  out <- broom::tidy(model)
  return(
    out %>%
      dplyr::mutate(
        term = stringr::str_replace(term, "variable",
          stringr::str_c(var, ": "))))

}

compute_oddsratio_CI <- function(data, model, var) {

  browser()
  confint <- questionr::odds.ratio(model)
  confint %>%
    tibble::as_tibble(rownames = "term") %>%
    set_names(c("term", "oddsratio", "low", "high", "p.value")) %>%
    dplyr::select(-p.value) %>%
    dplyr::mutate(
      term = stringr::str_replace(term, "variable",
        stringr::str_c(var, ": ")))

}

merge_results <- function(npart, pvalues, or_ci) {

  return(
    purrr::reduce(
      list(npart, pvalues, or_ci),
      dplyr::inner_join, by = "term") %>%
    dplyr::select(-estimate,-std.error) %>%
    dplyr::mutate_if(is.numeric, list( ~ round(.,4))) %>%
    dplyr::mutate( high = if_else(high > 1e5, Inf, high)))

}

ind_var_models %<>%
  dplyr::mutate(
    data = furrr::future_map(var, build_data, colo_cohort),
    yes = purrr::map(data, filter, colorectal == "yes") %>%
      purrr::map_int(nrow),
    no = purrr::map(data, filter, colorectal == "no") %>%
      purrr::map_int(nrow),
    n = purrr::map_int(data, nrow),
    npart = purrr::map2(data, var, count_cases))
    


biglm::bigglm(colorectal ~ variable, data = ind_var_models$data[[3]],
  family = binomial(link = "logit"), chunksize = 100e3)


# ind_var_models %<>%
#   dplyr::mutate(
#     model = purrr::future_map(data, ~ biglm::bigglm(
#       colorectal ~  variable, data = .x,
#         family = "binomial")))
        
        
#         ,
#     pvalues = purrr::map2(model, var, compute_pvalues),
#     or_ci = purrr::pmap(list(data, model, var), compute_oddsratio_CI),
#     all = purrr::pmap(list(npart, pvalues, or_ci), merge_results))

# ind_var_models %>%
#   dplyr::select(all) %>%
#   tidyr::unnest( cols = c(all)) %>%
#   tidyr::separate(term, into = c("variable", "level"), sep = ": ") %>%
#   DT::datatable(rownames = FALSE, options = list(pageLength = 20))


