# maaslin2_analysis.R
# Author: Guangyi Li
# Date: 2025-06-17
# Description: Perform MaAsLin2 analysis on microbial OTU data from multiple countries and sample splits; conducting linear modeling with and without FML adjustment, and save the results.

library(Maaslin2)
library(dplyr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- ifelse(length(args) >= 1, args[1], "raw_data_split/raw_data_strain")
#input_dir <- ifelse(length(args) >= 1, args[1], "raw_data_split/raw_data_species")
#input_dir <- ifelse(length(args) >= 1, args[1], "raw_data_split/raw_data_genus")

out_dir <- ifelse(length(args) >= 2, args[2], "maaslin_result_strain")
#out_dir <- ifelse(length(args) >= 2, args[2], "maaslin_result_species")
#out_dir <- ifelse(length(args) >= 2, args[2], "maaslin_result_genus")

countries <- c("austria", "china", "india", "france", "italy", "usa", "japan")
n_splits <- 100

run_maaslin <- function(otu, meta, out_path, include_fml = FALSE) {
  otu <- as.data.frame(otu)
  otu <- t(otu)
  colnames(otu) <- make.names(gsub(".*t__([^_]+_[^_]+).*", "\\1", colnames(otu)), unique = TRUE)
# colnames(otu_data) <- gsub(".*s__([^|]+).*", "\\1", colnames(otu_data))  
# colnames(otu_data) <- gsub(".*g__([^|]+).*", "\\1", colnames(otu_data))

  otu[] <- lapply(otu, function(x) as.numeric(as.character(x)))
  
  meta <- tibble::column_to_rownames(meta, "Sample_ID")
  meta$Disease <- factor(meta$Disease)
  meta$Gender <- factor(meta$Gender)
  meta$Age <- as.numeric(meta$Age)
  meta$BMI <- as.numeric(meta$BMI)
  meta$fml <- as.numeric(meta$fml)
  
  fixed <- if (include_fml) c("Disease", "Age", "Gender", "BMI", "fml") else c("Disease", "Age", "Gender", "BMI")
  
  Maaslin2(input_data = otu,
           input_metadata = meta,
           output = out_path,
           fixed_effects = fixed,
           normalization = "NONE",
           transform = "LOG",
           analysis_method = "LM",
           max_significance = 0.25,
           correction = "BH",
           standardize = TRUE,
           plot_heatmap = TRUE,
           heatmap_first_n = 50,
           plot_scatter = TRUE)
}

for (country in countries) {
  cat("Analyzing:", country, "\n")
  input_country <- file.path(input_dir, country)
  out_country <- file.path(out_dir, country)
  if (!dir.exists(out_country)) dir.create(out_country, recursive = TRUE)
  
  for (split_i in 1:n_splits) {
    split_path <- file.path(input_country, paste0("split_", split_i))
    meta_file <- file.path(split_path, "metadata_train.csv")
    otu_file <- file.path(split_path, "sylph_otu_train.csv")
    
    if (!file.exists(meta_file) || !file.exists(otu_file)) {
      cat("Missing files for", country, "split", split_i, "- skipping\n")
      next
    }
    
    meta <- read.csv(meta_file, check.names = FALSE)
    otu <- read.csv(otu_file, check.names = FALSE)
    
    run_maaslin(otu, meta, file.path(out_country, paste0("split_", split_i, "_before_fml_adj")), include_fml = FALSE)
    run_maaslin(otu, meta, file.path(out_country, paste0("split_", split_i, "_after_fml_adj")), include_fml = TRUE)
  }
}

cat("MaAsLin2 analysis completed for all countries.\n")

