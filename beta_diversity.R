# permanova_analysis.R
# Beta diversity and PERMANOVA analysis across countries at multiple taxonomy levels
# Author: Guangyi Li
# Date: 2025-06-17

library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pairwiseAdonis)
library(patchwork)
library(readr)

countries <- c("austria", "china", "india", "france", "usa", "italy", "japan")
taxonomy_levels <- c("genus", "species", "strain")
base_path <- "raw_data"  
result_base_dir <- "beta_diversity"

covariates_list <- list(
  after_fml_adj = c("Disease", "Age", "Gender", "BMI", "FML"),
  before_fml_adj = c("Disease", "Age", "Gender", "BMI")
)

run_permanova <- function(otu, meta, covariates) {
  common_samples <- intersect(colnames(otu), meta$Sample_ID)
  otu <- otu[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, meta$Sample_ID), ]
  otu <- otu[rowSums(otu > 0) >= 0.1 * ncol(otu), ]
  
  bray_dist <- vegdist(t(otu), method = "bray")
  formula <- as.formula(paste("bray_dist ~", paste(covariates, collapse = " + ")))
  set.seed(123)
  adonis2(formula, data = meta, permutations = 999, by = "margin")
}

create_pcoa_plot <- function(otu, meta, group_var = "Disease", country, tax_level) {
  bray_dist <- vegdist(t(otu), method = "bray")
  pcoa <- cmdscale(bray_dist, eig = TRUE)
  variance <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)
  
  plot_data <- data.frame(
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    Group = meta[[group_var]][match(colnames(otu), meta$Sample_ID)]
  )
  
  ggplot(plot_data, aes(PC1, PC2, color = Group)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(level = 0.8, linewidth = 0.5) +
    labs(
      x = paste0("PC1 (", variance[1], "%)"),
      y = paste0("PC2 (", variance[2], "%)"),
      title = paste(country, "-", tax_level),
      color = "Group"
    ) +
    theme_minimal(base_size = 10) +
    scale_color_manual(values = c("Control" = "#8ce99a", "CRC" = "#91a7ff")) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 10))
}

analyze_country_tax <- function(country, tax_level) {
  message("Processing ", country, " at level ", tax_level, " ...")
  otu_path <- file.path(base_path, tax_level, country, paste0(country, "_merged_sylph_output.tsv"))
  meta_path <- file.path(base_path, tax_level, country, paste0(country, "_metadata_full.csv"))
  result_dir <- file.path(result_base_dir, tax_level, country)
  if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

  otu <- read.table(otu_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_path, header = TRUE)

  colnames(otu) <- gsub("\\.", "-", colnames(otu))
  meta$Sample_ID <- gsub("\\.", "-", meta$Sample_ID)

  meta$Disease <- factor(meta$Disease)
  if ("Gender" %in% colnames(meta)) meta$Gender <- factor(meta$Gender)

  for (cov_type in names(covariates_list)) {
    covariates <- intersect(covariates_list[[cov_type]], colnames(meta))
    if (length(covariates) == 0) next

    message("  Running PERMANOVA with ", cov_type, " covariates ...")
    res <- tryCatch(run_permanova(otu, meta, covariates), error = function(e) {
      message("  PERMANOVA error: ", e$message)
      NULL
    })

    if (!is.null(res)) {
      sink(file.path(result_dir, paste0("permanova_results_", cov_type, ".txt")))
      cat("PERMANOVA Results for ", country, " - ", tax_level, " - ", cov_type, "\n\n")
      print(res)
      sink()

      if ("Disease" %in% covariates) {
        bray_dist <- vegdist(t(otu), method = "bray")
        pairwise <- pairwise.adonis(bray_dist, meta$Disease, p.adjust.m = "fdr")
        pairwise$sig <- cut(pairwise$p.adjusted, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                            labels = c("***", "**", "*", ".", ""), right = TRUE)
        write.csv(pairwise, file.path(result_dir, paste0("pairwise_comparisons_", cov_type, ".csv")), row.names = FALSE)
      }
    }
  }

  p <- create_pcoa_plot(otu, meta, "Disease", country, tax_level)
  ggsave(file.path(result_dir, "pcoa_plot.pdf"), p, width = 5, height = 5)
  message("Completed analysis for ", country, " at level ", tax_level)
  return(p)
}

# -------------------- Main --------------------
for (tax_level in taxonomy_levels) {
  message("\n========== Processing taxonomy level: ", tax_level, " ==========")
  all_plots <- lapply(countries, function(cn) analyze_country_tax(cn, tax_level))

  combined_pcoa <- patchwork::wrap_plots(all_plots, ncol = length(countries), guides = "collect") +
    patchwork::plot_annotation(
      title = paste("PCoA Plots - Bray-Curtis (", tax_level, " level)"),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
    ) & theme(legend.position = "bottom")

  ggsave(file.path(result_base_dir, tax_level, "combined_pcoa_plots.pdf"), combined_pcoa, width = 30, height = 5)
}

message("\nAll analyses completed successfully!")
message("Results saved to: ", normalizePath(result_base_dir))
