# Program: 01_manhattan_plot_with_snp_labels.R
# Author: Rene Welch
# Created: 2021-02-19
# Updated: 2021-05-03
# Purpose: Makes manhattan/qq plots, labelling peaks with annotated genes, it
#   also computes inflation factors

library(magrittr)
library(ggplot2)
library(scales)
library(cowplot)
library(qs)
library(rlang)
library(UWCCC.GWAS)
library(ggrepel, ggpubr)
library(tidyverse)

# general parameters

folder <- "2021_02_19_update_analysis_vitdbp_manuscript"
pvalue_thr <- 5e-8

figh <- 5
figw <- 12
units <- "in"


# load data and annotations
sccs_vdbp <-
  qs::qread(here::here("results", "zz_methods", "qs", "sccs_vdbp_stats.qs"))

sccs_vdbp %<>%
  dplyr::mutate(p = 2 * pnorm(-abs(t_stat)))

vitd_rsids <-
  qs::qread(here::here("results",
    "2021_02_19_update_analysis_vitdbp_manuscript", "qs", "chrom4_rsids.qs"))


sccs_vdbp_significant <- sccs_vdbp %>%
  dplyr::filter(p <= 5e-8)

sccs_vdbp_genes <-
  qs::qread(here::here("results", "zz_methods", "qs",
    "sccs_vdbp_significant_with_genes.qs")) %>%
  dplyr::select(id, genes) %>%
  tidyr::unnest(cols = c(genes))


sccs_vdbp_significant %<>%
  dplyr::inner_join(sccs_vdbp_genes, by = "id")

sccs_vdbp_significant0 <- sccs_vdbp_significant %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    chrom = unique(chrom),
    pos = pos[which.min(p)],
    p = min(p),
    .groups = "drop") %>%
  dplyr::top_n(2, wt = -p)
  
vdbp_manhattan_plot <- sccs_vdbp %>%
  manhattan_plot(
    sig = pvalue_thr, sig_color = "red", label_data = sccs_vdbp_significant0,
      label_var = symbol) +
  labs(
    x = NULL,
    color = "chr",
    y = expression(-log[10](p))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))

ggsave(
  filename = here::here("results", folder, "figs",
    "sccs_vitdbp_manhattan_plot_genes_label_p5e-8.png"),
  plot = vdbp_manhattan_plot,
  width = figw, height = figh, units = units)

median(qchisq(1 - sccs_vdbp$p, 1) / qchisq(.5, 1))

vdbp_qqplot <- sccs_vdbp %>%
  UWCCC.GWAS::qqplot_pvalue() +
   theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_blank())

ggsave(
  filename = here::here("results", folder, "figs",
    "sccs_vitdbp_qqplot.png"),
    plot = vdbp_qqplot, width = figh, height = figh, units = units)
