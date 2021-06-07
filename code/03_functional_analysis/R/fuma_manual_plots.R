# Program: fuma_manual_plots.R
# Author: Rene Welch
# Created: 2020-04-12
# Updated: 2021-04-12
# Purpose: Makes manual plots of FUMA results

library(magrittr)
library(tidyverse)
library(ComplexHeatmap)

folder <- "2021_04_12_fuma_manual_figs"
de_tissues <- here::here("results", folder, "FUMA_gene2func56796", 
  "gtex_v8_ts_general_DEG.txt") %>%
  readr::read_tsv()

de_tissues %<>%
  rename_with(snakecase::to_snake_case) %>%
  mutate(
    gene_set = fct_reorder(gene_set, p),
    tof = if_else(adj_p <= 0.05, "a", "b"))

library(tidytext)

de_plot <- de_tissues %>%
  filter(category != "DEG.twoside") %>%
  dplyr::mutate(
    category = factor(category, levels = c("DEG.up", "DEG.down")),
    category = fct_recode(category, "Up regulated" = "DEG.up",
      "Down regulated" = "DEG.down"),
    gene_set = reorder_within(gene_set, p, category),
    gene_set = fct_drop(gene_set)) %>%
  filter(p < 1) %>%
  ggplot(aes(x = gene_set, y = -log10(p), fill = tof)) +
  geom_col() +
  cowplot::theme_minimal_hgrid() +
  scale_x_reordered() +
  facet_wrap(vars(category), nrow  = 1, scales = "free_x") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(
    y = expression(-log[10](p.value)),
    x = "tissue gene set") +
  scale_fill_manual(values = c("a" = "lightblue", "b" = "grey")) +
  geom_hline(yintercept = -log10(.05), linetype = 2)

ggsave(
  filename = here::here("results", folder, "png", "upregulated_tissue.png"),
  plot = de_plot,
  units = "in", height = 3.5, width = 7)

log2expr <- here::here("results", folder, "FUMA_gene2func56193",
  "gtex_v8_ts_general_avg_log2TPM_exp.txt") %>%
  read_tsv() %>%
  select(-ensg) %>%
  as.data.frame() %>%
  column_to_rownames("symbol") %>%
  as.matrix()

ncols <- 11
cols <- circlize::colorRamp2(
  breaks = seq(0, max(log2expr), length.out = ncols),
  colors = rev(pals::brewer.spectral(ncols)))
  
Heatmap(log2expr, name = "expr", col = cols,
  show_row_dend = FALSE)
