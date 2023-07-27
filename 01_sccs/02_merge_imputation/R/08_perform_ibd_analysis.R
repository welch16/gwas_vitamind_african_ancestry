# Program: 08_perform_ibd_analysis.R
# Author: Rene Welch
# Created: 2020-09-15
# Updated: 2021-05-03
# Purpose: Review identity by descent calculations to check whether there were
#   related individuals in the SCCS cohort

library(magrittr)
library(tidyverse)
library(vroom)
library(qs)
library(cowplot)
library(ComplexHeatmap)

folder <- "2021_03_19_identity_by_descent"
ibd <- here::here("results", folder, "plink",
  "identity_by_descent.genome") %>%
  data.table::fread()
  
ibd %<>%
  tibble::as_tibble() %>%
  dplyr::rename_with(snakecase::to_snake_case) %>%
  dplyr::select(-tidyselect::starts_with("fid"))
  
theme_set(cowplot::theme_minimal_grid())
fs::dir_create(figsdir <<- here::here("results", folder, "figs"))

# pi_hat = P(IBD == 2) + 1/2 P(IBD == 1)

pihat_hist <- ibd %>%
  dplyr::filter(pi_hat > 0.1) %>%
  ggplot(aes(pi_hat)) +
  geom_histogram(binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = 0.1875, colour = "red", linetype = 2) +
  scale_x_continuous(breaks = scales::breaks_width(.05)) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1, scale = 1e-3, sufix = "K")) +
  labs(
    y = "number of pairs",
    x = "estimated pairwise IBD (pi hat)",
    title = "IBD for sample pairs with PI HAT > 0.1")

ggsave(
  filename = file.path(figsdir, "pihat_histogram.png"),
  plot = pihat_hist,
  width = 7, height = 3.5, units = "in")

ff <- plinkQC::relatednessFilter(ibd,
  relatednessTh = .35,
  relatednessIID1 = "iid_1",
  relatednessIID2 = "iid_2",
  relatednessRelatedness = "pi_hat")

qs::qsave(ff, here::here("results", folder, "qs",
  "relatedness_filter.qs"))
