################################################################################
# Random Forest Classification for Cluster Prediction
# 
# This script trains a Random Forest classifier to predict cluster assignments
# and evaluates the minimum number of genes (variable features or TFs) required
# for accurate classification. The analysis includes:
# - Training/testing split and model evaluation
# - Feature importance ranking
# - Accuracy vs. number of genes analysis
# - Minimal gene set identification
#
################################################################################

# Load Required Libraries ------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(randomForest)
  library(ggplot2)
  library(caret)
  library(scales)
})

# Set Paths and Parameters -----------------------------------------------------
output_dir <- "Path_to_output_dir"
data_dir <- "Path_to_data_dir"

# Analysis parameters
set.seed(123)  # For reproducibility
n_trees <- 500
train_split <- 0.7
max_genes_to_test <- 50
accuracy_threshold <- 0.95

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

################################################################################
# SECTION 1: LOAD DATA AND PREPARE METADATA
################################################################################

# cat("Loading datasets...\n")

# Load pre-processed Seurat objects
# PN_subset <- readRDS(file.path(data_dir, "20241101_PNSubset.rds"))
# OM_cldn_5dpf_scSEQ_20211209 <- readRDS(file.path(data_dir, "20241101_OM_cldn_5dpf_scSEQ_20211209.rds"))

# Load transcription factor list (optional - for TF-based analysis)
# All_ZF_TF <- read.csv('Path_toTF_data', header = TRUE)

# Extract Barcodes -------------------------------------------------------------
filtered_barcodes <- Cells(PN_subset)
all_barcodes <- Cells(OM_cldn_5dpf_scSEQ_20211209)

cat("Total cells in full dataset:", length(all_barcodes), "\n")
cat("Cells in PN subset:", length(filtered_barcodes), "\n")

# Mark Filtered Cells in Original Dataset -------------------------------------
OM_cldn_5dpf_scSEQ_20211209$Filtered <- ifelse(
  Cells(OM_cldn_5dpf_scSEQ_20211209) %in% filtered_barcodes,
  "Filtered",
  "Unfiltered"
)

# Map Cluster Labels from PN_subset to Full Dataset ---------------------------
cat("Mapping cluster labels...\n")

# Initialize cluster column
OM_cldn_5dpf_scSEQ_20211209$FilteredCluster <- NA

# Extract cluster labels from PN_subset
filtered_cluster_labels <- FetchData(PN_subset, vars = "seurat_clusters") %>%
  rownames_to_column(var = "Barcode")

# Create barcode to cluster mapping
barcode_mapping <- filtered_cluster_labels
rownames(barcode_mapping) <- barcode_mapping$Barcode

# Map clusters to full dataset
OM_cldn_5dpf_scSEQ_20211209$FilteredCluster <- barcode_mapping[
  Cells(OM_cldn_5dpf_scSEQ_20211209),
  "seurat_clusters"
]

# Add prefix to distinguish PN_subset cells
OM_cldn_5dpf_scSEQ_20211209$FilteredCluster <- ifelse(
  !is.na(OM_cldn_5dpf_scSEQ_20211209$FilteredCluster),
  paste0("PN_", OM_cldn_5dpf_scSEQ_20211209$FilteredCluster),
  "Unfiltered"
)

cat("Cluster distribution:\n")
print(table(OM_cldn_5dpf_scSEQ_20211209$FilteredCluster))

################################################################################
# SECTION 2: PREPARE EXPRESSION MATRIX
################################################################################

cat("\nPreparing expression matrix...\n")

# Normalize data for consistency
OM_cldn_5dpf_scSEQ_20211209 <- NormalizeData(OM_cldn_5dpf_scSEQ_20211209)

# Extract expression data (genes x cells) and transpose to (cells x genes)
expr_data <- t(as.matrix(GetAssayData(OM_cldn_5dpf_scSEQ_20211209, slot = "data")))
expr_data <- as.data.frame(expr_data)
expr_data$Cluster <- OM_cldn_5dpf_scSEQ_20211209$FilteredCluster

################################################################################
# SECTION 3: FEATURE SELECTION
################################################################################

# Option 1: Use Top Variable Features -----------------------------------------
cat("Using top 2000 variable features...\n")
top_genes <- head(VariableFeatures(OM_cldn_5dpf_scSEQ_20211209), 2000)
expr_data <- expr_data[, c(top_genes, "Cluster"), drop = FALSE]

cat("Features in dataset:", ncol(expr_data) - 1, "\n")

# Option 2: Use Transcription Factors Only (Alternative) ----------------------
# Uncomment this section to use TFs instead of variable features
#
# cat("Filtering for transcription factors...\n")
# valid_tfs <- intersect(All_ZF_TF$gene_symbol, colnames(expr_data))
# cat("Number of TFs found in dataset:", length(valid_tfs), "\n")
# expr_data <- expr_data[, c(valid_tfs, "Cluster"), drop = FALSE]

################################################################################
# SECTION 4: TRAIN/TEST SPLIT
################################################################################

cat("\nSplitting data into training and testing sets...\n")

# Create stratified split (70% train, 30% test)
train_indices <- createDataPartition(
  expr_data$Cluster,
  p = train_split,
  list = FALSE
)

train_data <- expr_data[train_indices, ]
test_data <- expr_data[-train_indices, ]

# Prepare predictor matrices and response vectors
x_train <- train_data[, setdiff(colnames(train_data), "Cluster")]
y_train <- as.factor(train_data$Cluster)

x_test <- test_data[, setdiff(colnames(test_data), "Cluster")]
y_test <- test_data$Cluster

cat("Training samples:", nrow(train_data), "\n")
cat("Testing samples:", nrow(test_data), "\n")

################################################################################
# SECTION 5: TRAIN RANDOM FOREST MODEL
################################################################################

cat("\nTraining Random Forest model...\n")

rf_model <- randomForest(
  x = x_train,
  y = y_train,
  ntree = n_trees,
  importance = TRUE
)

print(rf_model)

################################################################################
# SECTION 6: MODEL EVALUATION
################################################################################

cat("\nEvaluating model performance...\n")

# Predict on test set
predictions <- predict(rf_model, newdata = x_test)

# Generate confusion matrix
confusion_matrix <- table(Actual = y_test, Predicted = predictions)

# Calculate overall accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Overall Test Accuracy:", round(accuracy * 100, 2), "%\n")

# Predict cluster assignments for entire dataset
all_predictions <- predict(rf_model, newdata = expr_data[, -ncol(expr_data)])
OM_cldn_5dpf_scSEQ_20211209$PredictedCluster <- all_predictions

# Evaluate predictions specifically for PN_subset cells
pn_cells <- Cells(PN_subset)
pn_predictions <- OM_cldn_5dpf_scSEQ_20211209$PredictedCluster[pn_cells]
true_clusters <- OM_cldn_5dpf_scSEQ_20211209$FilteredCluster[pn_cells]

pn_confusion_matrix <- table(Actual = true_clusters, Predicted = pn_predictions)
pn_accuracy <- sum(diag(pn_confusion_matrix)) / sum(pn_confusion_matrix)

cat("PN Subset Accuracy:", round(pn_accuracy * 100, 2), "%\n")

################################################################################
# SECTION 7: FEATURE IMPORTANCE ANALYSIS
################################################################################

cat("\nAnalyzing feature importance...\n")

# Extract and rank feature importance scores
importance_scores <- as.data.frame(importance(rf_model))
importance_scores$Gene <- rownames(importance_scores)

# Rank by MeanDecreaseAccuracy
importance_scores <- importance_scores[
  order(importance_scores$MeanDecreaseAccuracy, decreasing = TRUE),
]

cat("Top 10 most important features:\n")
print(head(importance_scores[, c("Gene", "MeanDecreaseAccuracy")], 10))

################################################################################
# SECTION 8: ACCURACY vs NUMBER OF GENES ANALYSIS
################################################################################

cat("\nTesting accuracy with varying number of top genes...\n")

accuracy_by_genes <- c()

for (i in 1:max_genes_to_test) {
  # Select top i genes based on importance
  selected_genes <- importance_scores$Gene[1:i]
  
  # Subset data for selected genes
  x_subset <- expr_data[, c(selected_genes, "Cluster"), drop = FALSE]
  x_train_subset <- x_subset[train_indices, -ncol(x_subset), drop = FALSE]
  y_train_subset <- x_subset[train_indices, "Cluster"]
  x_test_subset <- x_subset[-train_indices, -ncol(x_subset), drop = FALSE]
  y_test_subset <- x_subset[-train_indices, "Cluster"]
  
  # Train random forest model with subset
  rf_subset <- randomForest(
    x = x_train_subset,
    y = as.factor(y_train_subset),
    ntree = n_trees
  )
  
  # Predict on test set
  predictions_subset <- predict(rf_subset, x_test_subset)
  
  # Generate confusion matrix
  confusion_matrix_single <- table(
    Actual = y_test_subset,
    Predicted = predictions_subset
  )
  
  # Subset confusion matrix for PN clusters only
  pn_confusion_matrix_subset <- confusion_matrix_single[
    grepl("PN_", rownames(confusion_matrix_single)),
    grepl("PN_", colnames(confusion_matrix_single)),
    drop = FALSE
  ]
  
  # Calculate PN cluster-specific accuracy
  pn_accuracy_subset <- sum(diag(pn_confusion_matrix_subset)) / 
    sum(pn_confusion_matrix_subset)
  
  accuracy_by_genes <- c(accuracy_by_genes, pn_accuracy_subset)
  
  if (i %% 10 == 0 || i <= 10) {
    cat("Genes:", i, "| PN Accuracy:", round(pn_accuracy_subset * 100, 2), "%\n")
  }
}

################################################################################
# SECTION 9: IDENTIFY MINIMAL GENE SET
################################################################################

cat("\nIdentifying minimal gene set...\n")

# Find number of genes needed for 95% of maximum accuracy
max_accuracy <- max(accuracy_by_genes)
minimal_genes <- which(accuracy_by_genes >= accuracy_threshold * max_accuracy)[1]

cat("Maximum accuracy achieved:", round(max_accuracy * 100, 2), "%\n")
cat("Minimal genes for", accuracy_threshold * 100, "% of max accuracy:", 
    minimal_genes, "\n")

# Extract minimal gene set
minimal_gene_set <- importance_scores$Gene[1:minimal_genes]

cat("\nMinimal gene set:\n")
print(minimal_gene_set)

################################################################################
# SECTION 10: VISUALIZATION
################################################################################

cat("\nGenerating plots...\n")

# Prepare data for plotting
accuracy_data <- data.frame(
  NumberOfGenes = 1:length(accuracy_by_genes),
  Accuracy = accuracy_by_genes
)

# Create publication-quality plot
# Plot accuracy vs. number of genes
plot(1:50, accuracy_by_genes, type = "b", xlab = "Number of Genes (var)", ylab = "PN Cluster Accuracy", main = "PN Cluster Accuracy vs. Number of Genes")

# Convert data to a data frame
accuracy_data <- data.frame(
  NumberOfGenes = 1:length(accuracy_by_genes),
  Accuracy = accuracy_by_genes
)

# Create the plot with corrected margin settings
PN_accuracy_plot <- ggplot(accuracy_data, aes(x = NumberOfGenes, y = Accuracy)) +
  geom_line(color = "blue", size = 1) +  # Line styling
  geom_point(color = "darkred", size = 3) +  # Point styling
  labs(
    x = "Number of Genes", 
    y = "PN Cluster Accuracy", 
    title = "PN Cluster Accuracy vs. Number of Genes"
  ) +
  theme_minimal(base_size = 14) +  # Use a clean theme with a larger font size
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = ggplot2::margin(b = 15)),  # Correct margin
    axis.text = element_text(size = 12),  # Larger axis text
    axis.title = element_text(size = 14),  # Larger axis labels
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Lighter major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    plot.margin = ggplot2::margin(20, 20, 20, 20)  # Add padding around the plot
  ) +
  scale_x_continuous(breaks = seq(0, max(accuracy_data$NumberOfGenes), by = 5)) +  # Custom x-axis breaks
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))   # Format y-axis as percentages


# Save plot
ggsave(
  filename = file.path(
    output_dir,
    "PN_cluster_accuracy_prediction_by_number_of_genes.svg"
  ),
  plot = PN_accuracy_plot,
  device = "svg",
  width = 8,
  height = 6
)

cat("Plot saved to:", output_dir, "\n")

################################################################################
# SECTION 11: SAVE RESULTS
################################################################################

cat("\nSaving results...\n")

# Save feature importance scores
write.csv(
  importance_scores,
  file = file.path(output_dir, "feature_importance_scores.csv"),
  row.names = FALSE
)

# Save accuracy by number of genes
write.csv(
  accuracy_data,
  file = file.path(output_dir, "accuracy_by_number_of_genes.csv"),
  row.names = FALSE
)

# Save minimal gene set
write.csv(
  data.frame(Gene = minimal_gene_set),
  file = file.path(output_dir, "minimal_gene_set.csv"),
  row.names = FALSE
)

# Save model performance summary
model_summary <- data.frame(
  Metric = c(
    "Overall Test Accuracy",
    "PN Subset Accuracy",
    "Maximum Accuracy",
    "Minimal Genes Required",
    "Threshold Used"
  ),
  Value = c(
    accuracy,
    pn_accuracy,
    max_accuracy,
    minimal_genes,
    accuracy_threshold
  )
)

write.csv(
  model_summary,
  file = file.path(output_dir, "model_performance_summary.csv"),
  row.names = FALSE
)

# Optional: Save the trained model
# saveRDS(rf_model, file = file.path(output_dir, "random_forest_model.rds"))

# Optional: Save updated Seurat object with predictions
# saveRDS(
#   OM_cldn_5dpf_scSEQ_20211209,
#   file = file.path(data_dir, "OM_cldn_5dpf_scSEQ_20211209_with_predictions.rds")
# )

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("\n=== Random Forest Classification Complete ===\n")
cat("All results saved to:", output_dir, "\n")
cat("\nKey Findings:\n")
cat("- Overall accuracy:", round(accuracy * 100, 2), "%\n")
cat("- PN cluster accuracy:", round(pn_accuracy * 100, 2), "%\n")
cat("- Minimal genes required:", minimal_genes, "\n")
cat("- Max accuracy achieved:", round(max_accuracy * 100, 2), "%\n")
