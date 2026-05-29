################################################################################
# Single-Cell RNA-Seq Analysis: Zebrafish Telencephalon at 5 dpf
# 
# This script performs comprehensive single-cell RNA-seq analysis including:
# - Quality control and preprocessing
# - Cell type identification and annotation
# - Subset analysis of neuronal populations
# - Visualization with UMAP and dot plots
# - Regional annotation (Olfactory Bulb vs Pallium)
#
# Author: Oded Mayseless
# Date: 12.02.2026
# 
################################################################################

# Load Required Libraries ------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggraph)
library(ggdendro)
library(ggrepel)
library(scCustomize)
library(RColorBrewer)
library(pheatmap)
library(scales)
library(tidyr)
library(seriation)
library(viridis)
source("2_Silhouette_scoring.R")
source("3_RF_cluster_classification.R")

# Set Working Directory and Paths ---------------------------------------------
save_dir <- "Path_to_output_dir"
data_dir <- "Path_to_data_dir"

# Create output directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

################################################################################
# SECTION 1: DATA LOADING AND QUALITY CONTROL
################################################################################

# Load 10X Genomics Data -------------------------------------------------------
cat("Loading 10X data...\n")
OM_cldn_5dpf_scSEQ_20211209.data <- Read10X(data.dir = data_dir)

# Initialize Seurat Object -----------------------------------------------------
cat("Creating Seurat object...\n")
OM_cldn_5dpf_scSEQ_20211209 <- CreateSeuratObject(
  counts = OM_cldn_5dpf_scSEQ_20211209.data,
  project = "OM_cldn_5dpf_scSEQ_20211209",
  min.cells = 3,
  min.features = 200
)

# Calculate Mitochondrial Percentage -------------------------------------------
OM_cldn_5dpf_scSEQ_20211209[["percent.mt"]] <- PercentageFeatureSet(
  OM_cldn_5dpf_scSEQ_20211209,
  pattern = "^mt-"
)

# Visualize QC Metrics ---------------------------------------------------------
cat("Generating QC plots...\n")
QC_vlnplot <- VlnPlot(
  OM_cldn_5dpf_scSEQ_20211209,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
QC_vlnplot
ggsave(
  filename = file.path(save_dir, "QC_vlnplot.svg"),
  plot = QC_vlnplot,
  device = "svg",
  width = 8,
  height = 8
)

# Filter Cells Based on QC Metrics --------------------------------------------
cat("Filtering cells...\n")
OM_cldn_5dpf_scSEQ_20211209 <- subset(
  OM_cldn_5dpf_scSEQ_20211209,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10
)

cat("Cells remaining after QC:", ncol(OM_cldn_5dpf_scSEQ_20211209), "\n")

################################################################################
# SECTION 2: NORMALIZATION AND DIMENSIONALITY REDUCTION
################################################################################

# Data Normalization -----------------------------------------------------------
cat("Normalizing data...\n")
OM_cldn_5dpf_scSEQ_20211209 <- NormalizeData(
  OM_cldn_5dpf_scSEQ_20211209,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify Highly Variable Features -------------------------------------------
cat("Finding variable features...\n")
OM_cldn_5dpf_scSEQ_20211209 <- FindVariableFeatures(
  OM_cldn_5dpf_scSEQ_20211209,
  selection.method = "vst",
  nfeatures = 2000
)

# Scale Data -------------------------------------------------------------------
cat("Scaling data...\n")
all.genes <- rownames(OM_cldn_5dpf_scSEQ_20211209)
OM_cldn_5dpf_scSEQ_20211209 <- ScaleData(
  OM_cldn_5dpf_scSEQ_20211209,
  features = all.genes
)

# Run PCA ----------------------------------------------------------------------
cat("Running PCA...\n")
OM_cldn_5dpf_scSEQ_20211209 <- RunPCA(
  OM_cldn_5dpf_scSEQ_20211209,
  features = VariableFeatures(object = OM_cldn_5dpf_scSEQ_20211209)
)

################################################################################
# SECTION 3: CLUSTERING AND UMAP VISUALIZATION
################################################################################

# Clustering -------------------------------------------------------------------
cat("Clustering cells...\n")
OM_cldn_5dpf_scSEQ_20211209 <- FindNeighbors(
  OM_cldn_5dpf_scSEQ_20211209,
  dims = 1:30
)

OM_cldn_5dpf_scSEQ_20211209 <- FindClusters(
  OM_cldn_5dpf_scSEQ_20211209,
  resolution = 0.2
)

# Run UMAP ---------------------------------------------------------------------
cat("Running UMAP...\n")
OM_cldn_5dpf_scSEQ_20211209 <- RunUMAP(
  OM_cldn_5dpf_scSEQ_20211209,
  dims = 1:30
)

# Find Cluster Markers ---------------------------------------------------------
cat("Finding cluster markers...\n")
OM_cldn_5dpf_scSEQ_20211209.markers <- FindAllMarkers(
  OM_cldn_5dpf_scSEQ_20211209,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract Top 10 Markers per Cluster ------------------------------------------
All_cluster_markers_res02 <- OM_cldn_5dpf_scSEQ_20211209.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

################################################################################
# SECTION 4: CELL TYPE ANNOTATION
################################################################################

# Annotate Clusters ------------------------------------------------------------
cat("Annotating clusters...\n")
Idents(OM_cldn_5dpf_scSEQ_20211209) <- OM_cldn_5dpf_scSEQ_20211209$seurat_clusters
OM_cldn_5dpf_scSEQ_20211209$old.ident <- Idents(OM_cldn_5dpf_scSEQ_20211209)

new.cluster.ids <- c(
  "Replicating cells_1",        # Cluster 0
  "Replicating cells_2",        # Cluster 1
  "Neuronal_01",                # Cluster 2
  "Neuronal_02",                # Cluster 3
  "Neuronal_03",                # Cluster 4
  "Glia",                       # Cluster 5
  "Neuronal_04",                # Cluster 6
  "Non-Neuronal_1",             # Cluster 7
  "Neuronal_05",                # Cluster 8
  "Non-Neuronal_2",             # Cluster 9
  "Non-Neuronal_3",             # Cluster 10
  "Neuronal_06",                # Cluster 11
  "Non-Neuronal_4",             # Cluster 12
  "Non-Neuronal_5"              # Cluster 13
)

names(new.cluster.ids) <- levels(OM_cldn_5dpf_scSEQ_20211209)

OM_cldn_5dpf_scSEQ_20211209 <- RenameIdents(
  OM_cldn_5dpf_scSEQ_20211209,
  new.cluster.ids
)

OM_cldn_5dpf_scSEQ_20211209$celltype <- plyr::mapvalues(
  x = as.character(OM_cldn_5dpf_scSEQ_20211209$seurat_clusters),
  from = names(new.cluster.ids),
  to = new.cluster.ids
)

################################################################################
# SECTION 5: VISUALIZATION - FULL DATASET
################################################################################

# Highlighted UMAP Plot --------------------------------------------------------
cat("Creating highlighted UMAP plot...\n")
clusters_to_highlight <- c(
  "Neuronal_01", "Neuronal_02", "Neuronal_03",
  "Neuronal_04", "Neuronal_05", "Neuronal_06", "Glia"
)

num_clusters <- length(clusters_to_highlight)
highlight_colors <- brewer.pal(n = num_clusters, name = "Set1")
names(highlight_colors) <- clusters_to_highlight

# Calculate plot limits
um <- Embeddings(OM_cldn_5dpf_scSEQ_20211209, "umap")
xlim <- quantile(um[, 1], c(0.01, 0.99))
xpad <- 0.03 * diff(xlim)
ylim <- quantile(um[, 2], c(0.01, 0.99))
ypad <- 0.03 * diff(ylim)
xlim <- xlim + c(-xpad, xpad)
ylim <- ylim + c(-ypad, ypad)


umap_plot_highlighted <- Cluster_Highlight_Plot(
  seurat_object = OM_cldn_5dpf_scSEQ_20211209,
  cluster_name = clusters_to_highlight,
  label = FALSE,
  highlight_color = highlight_colors
) +
  theme_void() +
  theme(plot.title = element_text(hjust = 1)) +
  NoLegend() +
  coord_cartesian(xlim = xlim, ylim = ylim, clip = "on")  # <-- crop & clip
umap_plot_highlighted <- LabelClusters(
  umap_plot_highlighted,
  id   = "ident",  # change if your ID differs
  size = 4,                  # try 5–8
  repel = TRUE
)
umap_plot_highlighted

ggsave(
  filename = file.path(save_dir, "highlighted_umap_telencephalic.svg"),
  plot = umap_plot_highlighted,
  device = "svg",
  width = 8,
  height = 8
)

# Dot Plot - General Markers ---------------------------------------------------
marker_genes <- c(
  'gfap', 'elavl3', 'snap25a', 'slc17a6a', 'slc17a6b',
  'slc32a1', 'gad1b', 'gad2', 'tal1', 'uncx',
  'pcna', 'her4.2', 'dla'
)

telencephalic_dotPlot_with_names <- DotPlot(
  OM_cldn_5dpf_scSEQ_20211209,
  features = marker_genes,
  cluster.idents = FALSE,
  dot.min = 0.05,
  scale.by = "size"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  NoLegend()

ggsave(
  filename = file.path(save_dir, "combined_telencephalic_dotPlot_with_names.svg"),
  plot = telencephalic_dotPlot_with_names,
  device = "svg",
  width = 10,
  height = 8
)

# Dot Plot with Legend ---------------------------------------------------------
telencephalic_dotPlot_with_names_With_Legend <- DotPlot(
  OM_cldn_5dpf_scSEQ_20211209,
  features = marker_genes,
  cluster.idents = FALSE,
  dot.min = 0.05,
  scale.by = "size"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  filename = file.path(save_dir, "combined_telencephalic_dotPlot_with_names_withLegend.svg"),
  plot = telencephalic_dotPlot_with_names_With_Legend,
  device = "svg",
  width = 8,
  height = 8
)

################################################################################
# SECTION 6: NEURONAL SUBSET ANALYSIS
################################################################################

# Subset Postmitotic Neurons ---------------------------------------------------
cat("Subsetting neuronal populations...\n")
Idents(OM_cldn_5dpf_scSEQ_20211209) <- OM_cldn_5dpf_scSEQ_20211209$seurat_clusters

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- subset(
  OM_cldn_5dpf_scSEQ_20211209,
  idents = c(2, 3, 4, 5, 6, 8, 11)
)

# Normalize and Process Neuronal Subset ----------------------------------------
cat("Processing neuronal subset...\n")
OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- NormalizeData(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- FindVariableFeatures(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  selection.method = "vst",
  nfeatures = 2000
)

all.genes <- rownames(OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset)
OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- ScaleData(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  features = all.genes
)

# PCA and Clustering -----------------------------------------------------------
OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- RunPCA(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  npcs = 30,
  features = VariableFeatures(object = OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset)
)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- FindNeighbors(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  dims = 1:30
)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- FindClusters(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  resolution = 0.1
)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- RunUMAP(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  dims = 1:30
)

# Find Markers for Neuronal Subset ---------------------------------------------
OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset_markers_res0_1 <- FindAllMarkers(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_cluster_markers_subset_res0_1 <- OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset_markers_res0_1 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

top3_cluster_markers_subset_res0_1 <- OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset_markers_res0_1 %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# Annotate Neuronal Subset Clusters --------------------------------------------
new.cluster.ids <- c(
  "Excitatory neurons",           # Cluster 0
  "Immature excitatory neurons",  # Cluster 1
  "Glia",                         # Cluster 2
  "Immature inhibitory neurons",  # Cluster 3
  "Inhibitory neurons",           # Cluster 4
  "Excitatory neurons_02",        # Cluster 5
  "Neuroendocrine",               # Cluster 6
  "Immature excitatory neurons_02", # Cluster 7
  "Inhibitory neurons_02"         # Cluster 8
)

names(new.cluster.ids) <- levels(OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset <- RenameIdents(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  new.cluster.ids
)

OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset$celltype <- plyr::mapvalues(
  x = as.character(OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset$seurat_clusters),
  from = names(new.cluster.ids),
  to = new.cluster.ids
)

# Dot Plot - Neuronal Subset ---------------------------------------------------
combined_markers <- c(
  'slc17a6a', 'slc17a6b', 'slc32a1', 'gad1b', 'gad2',
  'tal1', 'uncx', 'gfap', 'snap25a', 'dla', 'her4.2',
  top3_cluster_markers_subset_res0_1$gene
)

neuronal_subset_dotplot_combined_markers <- DotPlot(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  features = unique(combined_markers),
  dot.min = 0.05,
  scale.by = "size",
  col.min = 0,
  scale = TRUE
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  NoLegend()

ggsave(
  filename = file.path(save_dir, "neuronal_subset_dotplot_combined_markers.svg"),
  plot = neuronal_subset_dotplot_combined_markers,
  device = "svg",
  width = 10,
  height = 8
)

# Version with Legend ----------------------------------------------------------
neuronal_subset_dotplot_combined_markers_withLegend <- DotPlot(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  features = unique(combined_markers),
  dot.min = 0.05,
  scale.by = "size",
  col.min = 0,
  scale = TRUE
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(
  filename = file.path(save_dir, "neuronal_subset_dotplot_combined_markers_withLegend.svg"),
  plot = neuronal_subset_dotplot_combined_markers_withLegend,
  device = "svg",
  width = 8,
  height = 8
)

# Highlighted UMAP - Neuronal Subset -------------------------------------------
clusters_to_highlight <- c(
  "Excitatory neurons",
  "Inhibitory neurons",
  "Excitatory neurons_02",
  "Inhibitory neurons_02"
)

num_clusters <- length(clusters_to_highlight)
highlight_colors <- brewer.pal(n = num_clusters, name = "Set1")
names(highlight_colors) <- clusters_to_highlight

um <- Embeddings(OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset, "umap")
xlim <- quantile(um[, 1], c(0.01, 0.99))
xpad <- 0.03 * diff(xlim)
ylim <- quantile(um[, 2], c(0.01, 0.99))
ypad <- 0.03 * diff(ylim)
xlim <- xlim + c(-xpad, xpad)
ylim <- ylim + c(-ypad, ypad)

umap_plot_highlighted <- Cluster_Highlight_Plot(
  seurat_object = OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  cluster_name = clusters_to_highlight,
  label = FALSE,
  highlight_color = highlight_colors
) +
  theme_void() +
  theme(plot.title = element_text(hjust = 1)) +
  NoLegend() +
  coord_cartesian(xlim = xlim, ylim = ylim, clip = "on")

umap_plot_highlighted <- LabelClusters(
  umap_plot_highlighted,
  id = "ident",
  size = 4,
  repel = TRUE
)

ggsave(
  filename = file.path(save_dir, "highlighted_umap_postmitotic.svg"),
  plot = umap_plot_highlighted,
  device = "svg",
  width = 8,
  height = 8
)

################################################################################
# SECTION 7: EXCITATORY/INHIBITORY NEURON SUBSET ANALYSIS
################################################################################

# Subset Glutamatergic and GABAergic Neurons -----------------------------------
cat("Subsetting excitatory and inhibitory neurons...\n")
PN_subset <- subset(
  OM_cldn_5dpf_scSEQ_20211209.Neuronal_Subset,
  idents = clusters_to_highlight
)

# Process PN Subset ------------------------------------------------------------
PN_subset <- NormalizeData(PN_subset)
PN_subset <- FindVariableFeatures(PN_subset, selection.method = 'vst')
PN_subset <- ScaleData(PN_subset, features = rownames(PN_subset))
PN_subset <- RunPCA(PN_subset)
PN_subset <- FindNeighbors(PN_subset, dims = 1:30, verbose = FALSE)

# find optimal resolution with silouhette scoring
scData.combined <- ChooseClusterResolutionDownsample(PN_subset, 30,figdir = save_dir)
PN_subset <- FindClusters(PN_subset, resolution = 3, verbose = FALSE)
PN_subset <- RunUMAP(PN_subset, dims = 1:30)

# Find Markers for PN Subset ---------------------------------------------------
cat("Finding markers for PN subset...\n")
PN_subset_markers <- FindAllMarkers(
  PN_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

PN_subset_markers_top10_markers <- PN_subset_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.table(
  PN_subset_markers_top10_markers,
  file = file.path(save_dir, "Top10Var_cluster_markers_OB_subset_res3.csv")
)

# Excitatory/Inhibitory Marker Dot Plot ----------------------------------------
PN_ex_in_dotPlot <- DotPlot(
  object = PN_subset,
  features = c('slc17a6a', 'slc17a6b', 'slc32a1', 'gad1b', 'gad2'),
  cluster.idents = TRUE,
  dot.min = 0.2,
  scale.by = "size"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )


################################################################################
# SECTION 8: RF VALIDATAION of OLFACTORY BULB CLUSTERING
################################################################################

All_ZF_TF <-read.csv('Path_to_TF_data', header = TRUE)
TF_reults <- RF_cluster_classification(PN_subset, All_ZF_TF$gene_symbol, "Transcription Factors", output_dir = save_dir)
ribo_genes<-read.csv2('Path_to_ribo_data')
ribo_results<-RF_cluster_classification(PN_subset, ribo_genes$Gene_Name, "Ribosomal genes", output_dir = save_dir)



################################################################################
# SECTION 9: REGIONAL ANNOTATION (OLFACTORY BULB vs PALLIUM)
################################################################################

# Store Original Identities ----------------------------------------------------
Idents(PN_subset) <- PN_subset$seurat_clusters
PN_subset$old.ident <- Idents(PN_subset)

# Define Cluster Groups --------------------------------------------------------
inh_clusters <- c(23, 13, 21, 15, 18, 1, 6, 12, 5, 20, 24)
exc_clusters <- c(11, 7, 4, 19, 17, 16, 8, 10, 14, 2, 9, 0, 3, 22, 25)
OB_clusters <- c(6, 18, 24, 7, 11, 19, 14, 25)

# Create Abbreviated Annotations -----------------------------------------------
new.cluster.ids <- sapply(
  levels(PN_subset),
  function(cl) {
    region <- ifelse(as.numeric(cl) %in% OB_clusters, "OB", "Pa")
    ntype <- ifelse(
      as.numeric(cl) %in% exc_clusters, "e",
      ifelse(as.numeric(cl) %in% inh_clusters, "i", "Ot")
    )
    paste0(ntype, "-", region, " (", cl, ")")
  }
)

names(new.cluster.ids) <- levels(PN_subset)

# Rename and Store Metadata ----------------------------------------------------
PN_subset <- RenameIdents(PN_subset, new.cluster.ids)

PN_subset$celltype <- plyr::mapvalues(
  x = as.character(PN_subset$seurat_clusters),
  from = names(new.cluster.ids),
  to = new.cluster.ids
)

# Annotated UMAP Plot ----------------------------------------------------------
umap_plot_annotated <- DimPlot(
  PN_subset,
  reduction = "umap",
  group.by = "ident",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.3
) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(save_dir, "umap_OB_Pallium_annotated_3.svg"),
  plot = umap_plot_annotated,
  device = "svg",
  width = 6,
  height = 6
)

################################################################################
# SECTION 10: COMPREHENSIVE DOT PLOT WITH HIERARCHICAL CLUSTERING
################################################################################

# Define Genes to Plot ---------------------------------------------------------
genes_to_plot <- c(
  "dlx5a", "dlx6a", "dlx2a", "dlx2b", "pbx1a", "foxp4", "rpp25b",
  "gad1b", "gad2", "ddc", "sp8a", "slc32a1", "emx1", "emx2", "emx3",
  "lhx1a", "lhx2b", "lhx9", "tbr1a", "tbr1b", "slc17a6a", "slc17a6b",
  "barhl2", "lhx5", "zic1", "zic4", "zic5", "zic2a", "dla", "dlb",
  "eomesa", "neurod1", "zbtb18", "bhlhe22", "th", "drd2b", "tac3b",
  "scrt1b", "cort", "sst1.1"
)

genes_to_plot <- intersect(genes_to_plot, rownames(PN_subset))

# Build Cluster to Group Mapping ----------------------------------------------
clu_levels <- levels(PN_subset$seurat_clusters)
clu_num <- as.numeric(clu_levels)

is_OB <- clu_num %in% OB_clusters
is_ex <- clu_num %in% exc_clusters
is_in <- clu_num %in% inh_clusters

group_per_clu <- ifelse(
  is_OB & is_in, "i-OB",
  ifelse(is_OB & is_ex, "e-OB",
         ifelse(!is_OB & is_in, "i-Pa",
                ifelse(!is_OB & is_ex, "e-Pa", "Other")))
)

# Create Abbreviated Labels ----------------------------------------------------
label_map <- setNames(
  paste0(
    ifelse(grepl("^i", group_per_clu), "i",
           ifelse(grepl("^e", group_per_clu), "e", "Ot")),
    "-",
    ifelse(grepl("OB$", sub(".* ", "", group_per_clu)), "OB", "Pa"),
    " (", clu_levels, ")"
  ),
  clu_levels
)

# Order Clusters by Groups -----------------------------------------------------
groups_order <- c("i-OB", "e-OB", "i-Pa", "e-Pa")
ord_idx <- order(
  factor(group_per_clu, levels = groups_order, ordered = TRUE),
  clu_num
)

ordered_clusters <- clu_levels[ord_idx]
ordered_groups <- group_per_clu[ord_idx]

# Calculate Group Boundaries ---------------------------------------------------
n_in_ob <- sum(ordered_groups == "i-OB")
n_ex_ob <- sum(ordered_groups == "e-OB")
n_in_pa <- sum(ordered_groups == "i-Pa")
n_ex_pa <- sum(ordered_groups == "e-Pa")

separator_rows <- c(
  if (n_in_ob > 0) n_in_ob else NULL,
  if (n_ex_ob > 0) (n_in_ob + n_ex_ob) else NULL,
  if (n_in_pa > 0) (n_in_ob + n_ex_ob + n_in_pa) else NULL
)

# Create Dot Plot --------------------------------------------------------------
cat("Creating comprehensive dot plot...\n")
dp <- DotPlot(
  PN_subset,
  features = genes_to_plot,
  group.by = "seurat_clusters",
  dot.min = 0.05,
  scale.by = "size",
  scale = TRUE,
  col.min = 0
)

# Build Gene x Cluster Matrix from Dot Plot Data ------------------------------
mat <- dp$data %>%
  mutate(
    id = factor(id, levels = ordered_clusters),
    features.plot = as.character(features.plot)
  ) %>%
  filter(!is.na(id)) %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(
    names_from = id,
    values_from = avg.exp.scaled,
    values_fill = NA_real_
  ) %>%
  tibble::column_to_rownames("features.plot") %>%
  as.matrix()

# Filter Genes with Sufficient Data -------------------------------------------
keep <- rowSums(!is.na(mat)) > 1
mat <- mat[keep, , drop = FALSE]

sdrow <- apply(mat, 1, sd, na.rm = TRUE)
mat <- mat[sdrow > 0, , drop = FALSE]

# Hierarchical Clustering with Optimal Leaf Ordering --------------------------
C <- stats::cor(t(mat), use = "pairwise.complete.obs")
D <- as.dist(1 - C)
hc <- hclust(D, method = "ward.D2")
o <- seriate(D, method = "OLO", control = list(hclust = hc))
gene_order <- rownames(mat)[get_order(o)]

# Apply Gene Order to Dot Plot -------------------------------------------------
all_genes <- unique(dp$data$features.plot)
rest <- setdiff(all_genes, gene_order)
dp$data$features.plot <- factor(
  dp$data$features.plot,
  levels = c(gene_order, rest)
)

# Apply Cluster Order and Labels -----------------------------------------------
dp$data$id <- factor(
  dp$data$id,
  levels = ordered_clusters,
  labels = label_map[ordered_clusters]
)

# Clip Values for Visualization ------------------------------------------------
dp$data$avg.exp <- pmin(pmax(dp$data$avg.exp, 0), 3)
dp$data$pct.exp <- pmin(pmax(dp$data$pct.exp, 0), 100)

# Create Publication-Ready Dot Plot -------------------------------------------
dotplot_pub <- dp +
  scale_size_continuous(limits = c(0, 100), range = c(1.4, 5.8)) +
  labs(
    x = NULL,
    y = NULL,
    color = "Avg log1p expression",
    size = "% expressing"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.border = element_blank()
  )

# Add Separators Between Groups ------------------------------------------------
if (length(separator_rows)) {
  for (k in separator_rows) {
    dotplot_pub <- dotplot_pub +
      geom_hline(yintercept = k + 0.5, size = 0.4, color = "white")
  }
}

# Display and Save -------------------------------------------------------------
print(dotplot_pub)

ggsave(
  file.path(save_dir, "avgExp_dotplot_perCluster_OBPA_sep.svg"),
  plot = dotplot_pub,
  width = 9.5,
  height = 8.5,
  device = "svg"
)



################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("\n=== Analysis Complete ===\n")
cat("All plots saved to:", save_dir, "\n")

# Session Information ----------------------------------------------------------
cat("\nSession Information:\n")
sessionInfo()