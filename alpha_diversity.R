# alpha_diversity_analysis.R
# Author: Guangyi Li
# Date: 2025-06-17
# Description: Alpha diversity analysis (Shannon & Richness) across taxonomy levels and countries.

library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)

taxonomy_levels <- c("genus", "species", "strain")
countries <- c("austria", "china", "india", "france", "usa", "italy", "japan")

base_path <- "raw_data" 
result_dir <- "alpha_diversity"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

covariates_list <- list(
  after_fml_adj = c("Disease", "Age", "Gender", "BMI", "fml"),
  before_fml_adj = c("Disease", "Age", "Gender", "BMI")
)

# Function to plot boxplot with p-value
plot_box_with_p <- function(data, yvar, ylab, filename, country) {
  max_y <- max(data[[yvar]], na.rm = TRUE)
  ggplot(data, aes(x = Disease, y = .data[[yvar]], fill = Disease)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    stat_compare_means(method = "wilcox.test", label.y = max_y * 1.05) +
    labs(title = paste(ylab, "by Disease Group -", country), x = "Disease", y = ylab) +
    theme_minimal(base_size = 12) +
    scale_fill_manual(values = c("Control" = "#8ce99a", "CRC" = "#91a7ff")) +
    theme(legend.position = "none") -> p
  ggsave(filename, plot = p, width = 6, height = 5)
}

# Main analysis function for one country & taxonomy level
analyze_country_taxonomy <- function(country, tax_level) {
  message("?? Processing: ", country, " | Taxonomy Level: ", tax_level)

  # Set paths
  otu_path <- file.path(base_path, tax_level, country, paste0(country, "_merged_sylph_output.tsv"))
  meta_path <- file.path(base_path, tax_level, country, paste0(country, "_metadata_full.csv"))
  country_dir <- file.path(result_dir, tax_level, country)
  dir.create(country_dir, recursive = TRUE, showWarnings = FALSE)

  # fml data
  if (!file.exists(otu_path) || !file.exists(meta_path)) {
    message("?? Missing data for ", country, " - ", tax_level)
    return(NULL)
  }

  otu <- read.table(otu_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_path)
  meta$Sample_ID <- gsub("\\.", "-", meta$Sample_ID)
  colnames(otu) <- gsub("\\.", "-", colnames(otu))

  # Compute alpha diversity
  shannon <- diversity(t(otu), index = "shannon")
  richness <- specnumber(t(otu))

  # Merge with metadata
  shannon_data <- merge(data.frame(Sample_ID = colnames(otu), Shannon = shannon), meta, by = "Sample_ID")
  richness_data <- merge(data.frame(Sample_ID = colnames(otu), Richness = richness), meta, by = "Sample_ID")

  for (df in list(shannon_data, richness_data)) {
    df$Disease <- factor(df$Disease)
    if ("Gender" %in% colnames(df)) df$Gender <- factor(df$Gender)
  }

  # Linear models
  for (cov_type in names(covariates_list)) {
    covs <- intersect(covariates_list[[cov_type]], colnames(shannon_data))
    if (length(covs) == 0) next

    f_shannon <- reformulate(covs, response = "Shannon")
    f_richness <- reformulate(covs, response = "Richness")

    lm_shannon <- try(lm(f_shannon, data = shannon_data), silent = TRUE)
    lm_richness <- try(lm(f_richness, data = richness_data), silent = TRUE)

    sink(file.path(country_dir, paste0("regression_", cov_type, ".txt")))
    cat("Country:", country, "\nTaxonomy Level:", tax_level, "\nModel:", cov_type, "\n\n")

    if (!inherits(lm_shannon, "try-error")) {
      cat("Shannon Diversity Model:\n")
      print(summary(lm_shannon))
    } else cat("Shannon model failed.\n")

    if (!inherits(lm_richness, "try-error")) {
      cat("\nRichness Model:\n")
      print(summary(lm_richness))
    } else cat("Richness model failed.\n")
    sink()

    # Diagnostic plots
    if (!inherits(lm_shannon, "try-error")) {
      pdf(file.path(country_dir, paste0("diagnostic_shannon_", cov_type, ".pdf")))
      par(mfrow = c(2, 2))
      plot(lm_shannon)
      dev.off()
    }

    if (!inherits(lm_richness, "try-error")) {
      pdf(file.path(country_dir, paste0("diagnostic_richness_", cov_type, ".pdf")))
      par(mfrow = c(2, 2))
      plot(lm_richness)
      dev.off()
    }
  }

  # Boxplots
  plot_box_with_p(shannon_data, "Shannon", "Shannon Diversity Index",
                  file.path(country_dir, "boxplot_shannon.pdf"), country)
  plot_box_with_p(richness_data, "Richness", "Richness Index",
                  file.path(country_dir, "boxplot_richness.pdf"), country)

  # Output data
  write.csv(shannon_data, file.path(country_dir, "shannon_metadata.csv"), row.names = FALSE)
  write.csv(richness_data, file.path(country_dir, "richness_metadata.csv"), row.names = FALSE)

  message("Done: ", country, " / ", tax_level)
}

# Loop over taxonomy levels and countries
for (tax_level in taxonomy_levels) {
  for (country in countries) {
    analyze_country_taxonomy(country, tax_level)
  }
}

message("Alpha diversity analysis complete!")
message("Results saved to: ", normalizePath(result_dir))
