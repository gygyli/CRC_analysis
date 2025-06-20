# training_test_set_split.R
# Author: Guangyi Li
# Date: 2025-06-17
# Description: Split microbiome data of each country by taxonomy level into balanced training and test sets, repeated 100 times with random sampling.

library(dplyr)
library(caret)

# Function to balance the test set
balance_test_set <- function(metadata, test_ratio, seed) {
  set.seed(seed)
  test_indices <- createDataPartition(metadata$Disease, p = test_ratio, list = FALSE)
  metadata_test <- metadata[test_indices, ]
  metadata_train <- metadata[-test_indices, ]
  
  crc_count <- sum(metadata_test$Disease == "CRC")
  control_count <- sum(metadata_test$Disease == "Control")
  
  if (crc_count > control_count) {
    extra_samples <- metadata_train %>%
      filter(Disease == "Control") %>%
      sample_n(min(crc_count - control_count, n()), replace = FALSE)
  } else if (control_count > crc_count) {
    extra_samples <- metadata_train %>%
      filter(Disease == "CRC") %>%
      sample_n(min(control_count - crc_count, n()), replace = FALSE)
  } else {
    extra_samples <- data.frame()
  }
  
  metadata_test_balanced <- bind_rows(metadata_test, extra_samples)
  
  total_samples <- nrow(metadata)
  target_test_size <- round(total_samples * test_ratio)
  if (target_test_size %% 2 != 0) {
    target_test_size <- ifelse(abs((target_test_size + 1) / total_samples - test_ratio) <
                                 abs((target_test_size - 1) / total_samples - test_ratio),
                               target_test_size + 1, target_test_size - 1)
  }
  
  if (nrow(metadata_test_balanced) > target_test_size) {
    metadata_test_balanced <- metadata_test_balanced %>%
      group_by(Disease) %>%
      sample_n(target_test_size / 2) %>%
      ungroup()
  }
  
  metadata_train_balanced <- anti_join(metadata, metadata_test_balanced, by = "Sample_ID")
  
  return(list(train = metadata_train_balanced, test = metadata_test_balanced))
}

# Function to save data to folders
save_data_to_folder <- function(metadata_train, metadata_test, sylph_otu_train, sylph_otu_test, split_index, save_base_dir) {
  split_dir <- file.path(save_base_dir, paste0("split_", split_index))
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(metadata_train, file.path(split_dir, "metadata_train.csv"), row.names = FALSE)
  write.csv(metadata_test, file.path(split_dir, "metadata_test.csv"), row.names = FALSE)
  write.csv(sylph_otu_train, file.path(split_dir, "sylph_otu_train.csv"), row.names = FALSE)
  write.csv(sylph_otu_test, file.path(split_dir, "sylph_otu_test.csv"), row.names = FALSE)
}

# Main processing function supporting different taxonomy levels
process_country_taxonomy_data <- function(taxonomy_level, country, test_ratio = 0.2, n_splits = 100) {
  message("Processing taxonomy level: ", taxonomy_level, " country: ", country)
  
  base_path <- "raw_data"
  metadata_path <- file.path(base_path, taxonomy_level, country, paste0(country, "_metadata_full.csv"))
  otu_path <- file.path(base_path, taxonomy_level, country, paste0(country, "_merged_sylph_output.tsv"))
  
  metadata <- read.csv(metadata_path, check.names = FALSE)
  metadata$Sample_ID <- as.character(metadata$Sample_ID)
  
  sylph_otu <- read.table(otu_path, sep = "\t", header = TRUE, check.names = FALSE)
  colnames(sylph_otu) <- gsub("_", "", colnames(sylph_otu))
  
  # Updated save directory with raw_data_ prefix
  save_base_dir <- file.path("raw_data_split", paste0("raw_data_", taxonomy_level), country)
  
  for (split_index in 1:n_splits) {
    random_seed <- 42 + split_index
    balanced_data <- balance_test_set(metadata, test_ratio, random_seed)
    metadata_train <- balanced_data$train
    metadata_test <- balanced_data$test
    
    sylph_otu_train <- sylph_otu %>% select(cladename, metadata_train$Sample_ID)
    sylph_otu_test <- sylph_otu %>% select(cladename, metadata_test$Sample_ID)
    
    save_data_to_folder(metadata_train, metadata_test, sylph_otu_train, sylph_otu_test, split_index, save_base_dir)
  }
}

# === Main execution ===
taxonomy_levels <- c("genus", "species", "strain")
countries <- c("austria", "china", "india", "france", "usa", "italy", "japan")

for (taxonomy_level in taxonomy_levels) {
  for (country in countries) {
    process_country_taxonomy_data(taxonomy_level, country, test_ratio = 0.2, n_splits = 100)
  }
}

