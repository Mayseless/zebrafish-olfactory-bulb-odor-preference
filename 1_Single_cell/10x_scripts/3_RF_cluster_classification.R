# Load required packages
library(Seurat)
library(randomForest)
library(dplyr)
library(ggplot2)
library(caret)
library(reshape2)

# Define the function
RF_cluster_classification <- function(seurat_obj, gene_list, gene_list_name, output_dir = NULL) {
  
  # Set default assay to RNA
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Extract cluster identities
  clusters <- Idents(seurat_obj)
  
  # Extract gene names from the Seurat object
  genes_in_seurat <- rownames(seurat_obj)
  
  # Identify genes present in both the gene list and the dataset
  selected_genes <- intersect(gene_list, genes_in_seurat)
  
  # Extract normalized expression data for selected genes
  expr_data <- GetAssayData(seurat_obj, layer = "data")[selected_genes, ]
  expr_data <- as.data.frame(t(expr_data))
  
  # Add cluster identities
  expr_data$Cluster <- as.factor(clusters)
  
  # Split data into training and testing sets
  set.seed(123)
  train_indices <- sample(1:nrow(expr_data), size = 0.7 * nrow(expr_data))
  train_data <- expr_data[train_indices, ]
  test_data <- expr_data[-train_indices, ]
  
  # Prepare predictor matrix and response vector
  x_train <- train_data[, setdiff(colnames(train_data), 'Cluster')]
  y_train <- train_data$Cluster
  x_test <- test_data[, setdiff(colnames(test_data), 'Cluster')]
  y_test <- test_data$Cluster
  
  # Train the random forest model
  rf_model <- randomForest(x = x_train, y = y_train, importance = TRUE, ntree = 500)
  
  # Make predictions on the test set
  predictions <- predict(rf_model, newdata = x_test)
  
  # Calculate accuracy
  confusion_matrix <- table(Actual = y_test, Predicted = predictions)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  print(paste("Accuracy:", round(accuracy * 100, 2), "%"))
  
  # Extract top 50 important genes
  importance_scores <- importance(rf_model)
  top_genes <- importance_scores[order(importance_scores[, 'MeanDecreaseGini'], decreasing = TRUE), ]
  top_50_genes <- head(top_genes, 50)

  # Extract top 20 genes by MeanDecreaseAccuracy for plotting
  top_20_accuracy <- importance_scores[order(importance_scores[, 'MeanDecreaseAccuracy'], decreasing = TRUE), ][1:20, , drop = FALSE]
  top_20_accuracy <- as.data.frame(top_20_accuracy)
  top_20_accuracy$Gene <- rownames(top_20_accuracy)
  
  # Create a dot plot to replicate the style of varImpPlot
  accuracy_plot <- ggplot(top_20_accuracy, aes(x = MeanDecreaseAccuracy, y = reorder(Gene, MeanDecreaseAccuracy))) +
    geom_point(color = "darkblue", size = 3) +
    geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
    labs(title = "Top 20 Important Features by MeanDecreaseAccuracy",
         x = "MeanDecreaseAccuracy", y = "Gene") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
    
  importance_plot<-varImpPlot(rf_model, n.var = 20, main = "Top 20 Important Features")
  
  # Calculate per-cluster accuracy
  results_df <- data.frame(Actual = y_test, Predicted = predictions)
  per_cluster_accuracy <- results_df %>%
    group_by(Actual) %>%
    summarise(Total = n(), Correct = sum(Actual == Predicted), Accuracy = Correct / Total * 100) %>%
    arrange(desc(Accuracy))
  
  # Generate normalized confusion matrix for heatmap
  cm_df <- as.data.frame(confusion_matrix)
  colnames(cm_df) <- c("Actual", "Predicted", "Freq")
  cm_df <- cm_df %>%
    group_by(Actual) %>%
    mutate(Proportion = Freq / sum(Freq))
  
  # Generate misclassification heatmap data
  misclassified <- results_df %>%
    filter(Actual != Predicted)
  misclass_table <- table(misclassified$Actual, misclassified$Predicted)
  misclass_df <- as.data.frame(misclass_table)
  colnames(misclass_df) <- c("Actual", "Predicted", "Freq")
  
  # Plot normalized heatmap with dynamic title
  normalized_heatmap <- ggplot(data = cm_df, aes(x = Predicted, y = Actual, fill = Proportion)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkred") +
    labs(title = paste("Normalized Confusion Matrix Heatmap -", gene_list_name),
         x = "Predicted Cluster", y = "Actual Cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
  
  
  # Plot misclassification heatmap with dynamic title
  misclass_heatmap <- ggplot(data = misclass_df, aes(x = Predicted, y = Actual, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkblue") +
    labs(title = paste("Misclassification Heatmap -", gene_list_name),
         x = "Predicted Cluster", y = "Actual Cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
  
  # Save plots if output_dir is specified
  if (!is.null(output_dir)) {
    ggsave(filename = file.path(output_dir, paste0(gene_list_name, "_normalized_heatmap.svg")),
           plot = normalized_heatmap, device = "svg", width = 10, height = 6)
    ggsave(filename = file.path(output_dir, paste0(gene_list_name, "_misclass_heatmap.svg")),
           plot = misclass_heatmap, device = "svg", width = 10, height = 6)
    ggsave(filename = file.path(output_dir, paste0(gene_list_name, "_top20_MeanDecreaseAccuracy.svg")),
           plot = accuracy_plot, device = "svg", width = 8, height = 6)
  }
  
  # Return accuracy, top genes, and plots
  return(list(accuracy = accuracy, top_50_genes = top_50_genes, normalized_heatmap = normalized_heatmap, 
              misclass_heatmap = misclass_heatmap, accuracy_plot = accuracy_plot,importance_plot = importance_plot))
}

# save_dir = "Path_to_output_dir"
# 
# TF_reults <- RF_cluster_classification(PN_subset, All_ZF_TF$gene_symbol, "Transcription Factors", output_dir = save_dir)
# #ribo_genes<-read.csv2('Path_to_ribo_genes')
# ribo_results<-RF_cluster_classification(PN_subset, ribo_genes$Gene_Name, "Ribosomal genes", output_dir = save_dir)
# 
