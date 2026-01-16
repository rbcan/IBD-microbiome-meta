library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(targets)

source("IBD_strains_functions.R")

list(
  tar_target(
    metadata,
    filter(read.delim("data/map.0_Copy.txt"), !grepl("#", `X.SmplID`))
  ),
  tar_target(
    ts2,
    read_excel("data/Table_S2_overview.MGS.DA_strains_cor_genes.xlsx")
  ),
  tar_target(
    strain_res,
    read_in_strain_res("Z:/1/115b210e-aa45-4469-9eb1-d0d85d879fc0/data/GC_IBD8/Bin5SB/intra_phylo2/")
  ),
  tar_target(
    strain_res_df,
    left_join(strain_res, 
              unique(select(metadata, c("X.SmplID", "IndividualID", "Cohort", "IBD", "IBD2", "Time_rel", "IndividualID"))), 
              by = c("Sample" = "X.SmplID"))
  ),
  tar_target(
    pers_df,
    calculate_persistence(strain_res_df)
  ),
  tar_target(
    cor_df,
    get_df_for_correlations(pers_df, ts2)
  ),
  tar_target(
    fig4a,
    plot_persistence_per_ibd_association(cor_df)
  ),
  tar_target(
    fig4b,
    plot_ibd_association_correlations(cor_df, min_ind = 2)
  ),
  tar_target(
    mgs_persistent_ibd_enriched,
    get_species_persistent_enriched_in(cor_df, ts2, enriched_in = "IBD", min_ind = 2, persistence_threshold = 95)
  ),
  tar_target(
    mgs_persistent_control_enriched,
    get_species_persistent_enriched_in(cor_df, ts2, enriched_in = "control", min_ind = 2, persistence_threshold = 95)
  ),
  tar_target(
    fig4c,
    plot_strain_changes_definitions()
  ),
  tar_target(
    first5timepoints,
    get_first_x_timepoints_samples(metadata, n_timepoints = 5)
  ),
  tar_target(
   strain_res_df_5t,
   left_join(filter(strain_res, Sample %in% first5timepoints),
             unique(select(metadata, c("X.SmplID", "IndividualID", "Cohort", "IBD", "IBD2", "Time_rel", "IndividualID"))),
             by = c("Sample" = "X.SmplID"))
  ),
  tar_target(
    strain_changes_5t,
    get_strain_changes_df_5t(strain_res_df_5t, metadata, first5timepoints)
  ),
  tar_target(
    strain_changes_5t_df,
    summarise_strain_changes(strain_changes_5t, metadata, by_cohort = FALSE)
  ),
  tar_target(
    strain_changes_5t_df_coh,
    summarise_strain_changes(strain_changes_5t, metadata, by_cohort = TRUE)
  ),
  tar_target(
    strain_changes_test_coh_blocked,
    get_strain_changes_test_results_blocked_by_cohort(strain_changes_5t_df_coh, min_ind = 5)
  ),
  tar_target(
    fig4d,
    plot_strain_changes_boxplots_all(strain_changes_5t_df, strain_changes_test_coh_blocked)
  ),
  tar_target(
    figs15a,
    plot_ibd_association_correlations_ibd_healthy(pers_df, ts2, min_ind = 2)
  ),
  tar_target(
    figs15b,
    plot_ibd_association_correlations_per_cohort(pers_df, ts2, min_ind = 2)
  ),
  tar_target(
    figs16,
    plot_strain_changes_boxplots_per_cohort(strain_changes_5t_df_coh)
  ),
  tar_target(
    fig4_combined,
    combine_fig4(fig4a, fig4b, fig4c, fig4d)
  ),
  tar_target(
    fig4_png,
    ggsave("Fig4_08.12.2025.png", plot = fig4_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    fig4_pdf,
    ggsave("Fig4_08.12.2025.pdf", plot = fig4_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    fig4_svg,
    ggsave("Fig4_08.12.2025.svg", plot = fig4_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    figs15_combined,
    plot_grid(plotlist = list(figs15a, figs15b), labels = c("A", "B"), ncol = 1, rel_heights = c(4,7))
  ),
  tar_target(
    figs15_png,
    ggsave("FigS15_08.12.2025.png", plot = figs15_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    figs15_pdf,
    ggsave("FigS15_08.12.2025.pdf", plot = figs15_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    figs15_svg,
    ggsave("FigS15_08.12.2025.svg", plot = figs15_combined, width = 8, height = 12),
    format = "file"
  ),
  tar_target(
    figs16_png,
    ggsave("FigS16_08.12.2025.png", plot = figs16, width = 6, height = 6),
    format = "file"
  ),
  tar_target(
    figs16_pdf,
    ggsave("FigS16_08.12.2025.pdf", plot = figs16, width = 6, height = 6),
    format = "file"
  ),
  tar_target(
    figs16_svg,
    ggsave("FigS16_08.12.2025.svg", plot = figs16, width = 6, height = 6),
    format = "file"
  )
)
