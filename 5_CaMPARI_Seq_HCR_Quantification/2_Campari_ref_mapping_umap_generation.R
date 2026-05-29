


library(ggplot2)
library(dplyr)
library(tidyr)
# ----------------------------------
# Step 1: Reference mapping for each condition
# ----------------------------------
combined_FS_data_TCA <- subset(OM_20240412_CaMPARI_FS, subset = sample == "TCA")
combined_FS_data_ATP<- subset(OM_20240412_CaMPARI_FS, subset = sample == "ATP")
combined_FS_data_Cyt<- subset(OM_20240412_CaMPARI_FS, subset = sample == "Cytidine")
combined_FS_data_E3<- subset(OM_20240412_CaMPARI_FS, subset = sample == "E3")

All_ZF_TF <-read.csv('C:/Oded_data/single_cell_seq/data/zf_tfs_with_motifs.csv', header = TRUE)

ref_mapping_anchors<-c(All_ZF_TF$gene_symbol, PN_subset_markers_top3_markers$gene)
CaMPARI_FS_ref_TCA <- ref_mapping(PN_subset, combined_FS_data_TCA, features_to_anchor = ref_mapping_anchors, ref_dataset_title = "PN_subset", query_dataset_title = "TCA campari FS")#, ref_label_col = "seurat_clusters")
CaMPARI_FS_ref_ATP <- ref_mapping(PN_subset, combined_FS_data_ATP, features_to_anchor = ref_mapping_anchors, ref_dataset_title = "PN_subset", query_dataset_title = "ATP campari FS")#, ref_label_col = "seurat_clusters")
CaMPARI_FS_ref_Cyt <- ref_mapping(PN_subset, combined_FS_data_Cyt, features_to_anchor = ref_mapping_anchors, ref_dataset_title = "PN_subset", query_dataset_title = "Cytidine campari FS")#, ref_label_col = "seurat_clusters")
CaMPARI_FS_ref_E3  <- ref_mapping(PN_subset, combined_FS_data_E3, features_to_anchor = ref_mapping_anchors, ref_dataset_title = "PN_subset", query_dataset_title = "E3 campari FS")#, ref_label_col = "seurat_clusters")


library(dplyr)
library(ggplot2)
library(tidyr)   # <-- add this

# -------------------------------------------------------------------
# 1. Extract cluster assignment and metadata
# -------------------------------------------------------------------
get_cluster_counts <- function(mapped_obj, sample_name) {
  df <- data.frame(
    cluster   = as.character(mapped_obj$predicted.celltype),
    replicate = as.character(mapped_obj$replicate),
    sample    = as.character(if (!is.null(mapped_obj$sample)) mapped_obj$sample else sample_name),
    stringsAsFactors = FALSE
  )
  df %>%
    dplyr::group_by(sample, replicate, cluster) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
}


counts_TCA <- get_cluster_counts(CaMPARI_FS_ref_TCA$mapped_query, "TCA")
counts_ATP <- get_cluster_counts(CaMPARI_FS_ref_ATP$mapped_query, "ATP")
counts_Cyt <- get_cluster_counts(CaMPARI_FS_ref_Cyt$mapped_query, "Cytidine")
counts_E3  <- get_cluster_counts(CaMPARI_FS_ref_E3$mapped_query,  "E3")

all_counts <- bind_rows(counts_TCA, counts_ATP, counts_Cyt, counts_E3)


suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr)})

# helper: extract confidence scores from mapped Seurat object
get_conf_df <- function(mapped_obj, sample_name) {
  md <- mapped_obj[[]]
  cand <- c("prediction.score.max","predicted.id.score","predicted.celltype.score",
            "prediction.score","pred_score","score","predicted.score.max")
  score_col <- cand[cand %in% colnames(md)][1]
  if (is.na(score_col)) stop("Could not find a prediction score column in mapped_obj metadata.")
  tibble(
    sample    = if ("sample" %in% colnames(md)) as.character(md$sample) else sample_name,
    replicate = if ("replicate" %in% colnames(md)) as.character(md$replicate) else NA_character_,
    cluster   = as.character(md$predicted.celltype),
    score     = as.numeric(md[[score_col]])
  )
}

conf_TCA <- get_conf_df(CaMPARI_FS_ref_TCA$mapped_query, "TCA")
conf_ATP <- get_conf_df(CaMPARI_FS_ref_ATP$mapped_query, "ATP")
conf_Cyt <- get_conf_df(CaMPARI_FS_ref_Cyt$mapped_query, "Cytidine")
conf_E3  <- get_conf_df(CaMPARI_FS_ref_E3$mapped_query,  "E3")

conf_df <- bind_rows(conf_TCA, conf_ATP, conf_Cyt, conf_E3) %>%
  filter(!is.na(cluster), !is.na(score)) %>%
  mutate(
    # nice numeric ordering if clusters are numbers
    cn = suppressWarnings(as.numeric(cluster))
  )
levs <- conf_df %>% arrange(cn, cluster) %>% pull(cluster) %>% unique()
conf_df <- conf_df %>% mutate(cluster = factor(cluster, levels = levs))

# colors (use your custom palette if defined)
if (!exists("custom_colors")) {
  custom_colors <- setNames(scales::hue_pal()(length(unique(conf_df$sample))),
                            sort(unique(conf_df$sample)))
}

# plot: violin per cluster, faceted by sample
p_conf <- ggplot(conf_df, aes(x = cluster, y = score, fill = sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.9, alpha = 0.9, color = NA) +
  stat_summary(fun = median, geom = "point", size = 1.4, color = "black") +
  facet_wrap(~ sample, ncol = 1) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(x = "Cluster", y = "Mapping confidence (prediction score)",
       title = "CaMPARI-seq mapping confidence per cluster per sample") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(face = "bold"))
print(p_conf)

# minimal printout for interpretation (copy/paste friendly)
conf_df %>%
  group_by(sample, cluster) %>%
  summarise(n = dplyr::n(),
            median = median(score, na.rm = TRUE),
            iqr = IQR(score, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(sample, as.numeric(as.character(cluster))) %>%
  print(n = Inf)

# -------------------------------------------------------------------
# 2. Relative enrichment Odor vs E3 (per replicate)
# -------------------------------------------------------------------
# ---------- Odor vs E3 enrichment per replicate ----------
compare_vs_E3 <- function(odor_counts, e3_counts, odor_name, pseudocount = 0.5) {
  reps <- intersect(unique(odor_counts$replicate), unique(e3_counts$replicate))
  
  res_list <- lapply(reps, function(rep_i) {
    df_odor <- dplyr::filter(odor_counts, replicate == rep_i)
    df_e3   <- dplyr::filter(e3_counts,   replicate == rep_i)
    
    merged <- dplyr::full_join(df_e3, df_odor, by = "cluster",
                               suffix = c("_E3", "_Odor")) |>
      dplyr::mutate(
        n_E3   = tidyr::replace_na(n_E3,   0L),
        n_Odor = tidyr::replace_na(n_Odor, 0L)
      )
    
    N1 <- sum(merged$n_E3)   # total E3 in this replicate
    N2 <- sum(merged$n_Odor) # total Odor in this replicate
    
    out <- merged |>
      dplyr::rowwise() |>
      dplyr::mutate(
        # Fisher on integer counts only; handle degenerate totals
        fisher_p = {
          if (N1 == 0L || N2 == 0L) {
            NA_real_
          } else {
            a <- as.integer(n_Odor)
            b <- as.integer(N2 - n_Odor)
            c <- as.integer(n_E3)
            d <- as.integer(N1 - n_E3)
            # Defensive clamp (shouldn't be needed if inputs are consistent)
            a <- max(a, 0L); b <- max(b, 0L); c <- max(c, 0L); d <- max(d, 0L)
            
            mtx <- matrix(c(a, b, c, d),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("Odor","E3"), c("in_k","not_in_k")))
            pval <- tryCatch(stats::fisher.test(mtx, alternative = "greater")$p.value,
                             error = function(e) NA_real_)
            pval
          }
        },
        # Smoothed log2 enrichment (pseudocount only affects ratio, not Fisher)
        log2_ratio = {
          p_odor <- (n_Odor + pseudocount) / (N2 + 2 * pseudocount)
          p_e3   <- (n_E3   + pseudocount) / (N1 + 2 * pseudocount)
          log2(p_odor / p_e3)
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        sample = odor_name,
        replicate = rep_i,
        cluster,
        n_E3, n_Odor,
        N_E3 = N1, N_Odor = N2,
        log2_ratio, fisher_p
      )
    
    out
  })
  
  dplyr::bind_rows(res_list) |>
    dplyr::group_by(sample, replicate) |>
    dplyr::mutate(fisher_q = p.adjust(fisher_p, method = "BH")) |>
    dplyr::ungroup()
}

res_TCA <- compare_vs_E3(counts_TCA, counts_E3, "TCA")
res_ATP <- compare_vs_E3(counts_ATP, counts_E3, "ATP")
res_Cyt <- compare_vs_E3(counts_Cyt, counts_E3, "Cytidine")
enrichment_per_rep <- bind_rows(res_TCA, res_ATP, res_Cyt)

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Your colors
custom_colors <- c(
  "TCA" = "#004000",
  "Cytidine" = "#00C000",
  "ATP" = "#FF00FF"
)

# ---------- Prep: ensure types ----------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Ensure odor order and names
enrichment_per_rep <- enrichment_per_rep %>%
  mutate(
    sample   = factor(as.character(sample), levels = c("TCA","Cytidine","ATP")),
    cluster  = as.character(cluster),
    replicate= as.character(replicate)
  )

# Pooled summary (mean across replicates) + SE + min FDR across reps for stars
enrichment_summary <- enrichment_per_rep %>%
  group_by(sample, cluster) %>%
  summarise(
    mean_log2 = mean(log2_ratio, na.rm = TRUE),
    sd_log2   = sd(log2_ratio,   na.rm = TRUE),
    n         = sum(!is.na(log2_ratio)),
    min_q     = suppressWarnings(min(fisher_q, na.rm = TRUE)),
    .groups   = "drop"
  ) %>%
  mutate(
    se_log2 = ifelse(n > 1 & is.finite(sd_log2), sd_log2 / sqrt(n), 0),
    sig_label = case_when(
      is.finite(min_q) & min_q < 0.001 ~ "***",
      is.finite(min_q) & min_q < 0.01  ~ "**",
      is.finite(min_q) & min_q < 0.05  ~ "*",
      TRUE ~ ""
    ),
    # nice numeric ordering if clusters are numeric-like
    cluster_num = suppressWarnings(as.numeric(cluster))
  ) %>%
  group_by(cluster) %>%
  mutate(cluster_ord = unique(replace_na(cluster_num, Inf))) %>%
  ungroup()

# Order clusters left→right by numeric if possible
lvl_order <- enrichment_summary %>%
  distinct(cluster, cluster_ord) %>%
  arrange(cluster_ord, cluster) %>%
  pull(cluster)

enrichment_summary <- enrichment_summary %>%
  mutate(cluster = factor(cluster, levels = lvl_order))

# Colors
custom_colors <- c("TCA"="#004000","Cytidine"="#00C000","ATP"="#FF00FF")

# ---- BAR PLOT: pooled across replicates ----
pd <- position_dodge2(width = 0.8, preserve = "single")

p_bar_pooled <- ggplot(enrichment_summary,
                       aes(x = cluster, y = mean_log2, fill = sample)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # stars (placed slightly above the error bar top)
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = custom_colors, name = "Odor", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Subtype enrichment vs E3 (pooled across replicates)",
    x = "Cluster",
    y = "log2 enrichment (Odor / E3)"
  )

print(p_bar_pooled)

# focus on selected clusters with same width behavior ----
keep_clusters <- c("7","11","14","19","25","6","18","24")
enrichment_summary_subset <- enrichment_summary %>%
  filter(cluster %in% keep_clusters) %>%
  mutate(cluster = factor(cluster, levels = keep_clusters))  # keep explicit order

p_bar_subset <- ggplot(enrichment_summary_subset,
                       aes(x = cluster, y = mean_log2, fill = sample)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = custom_colors, name = "Odor", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Selected clusters: enrichment vs E3 (pooled)",
    x = "Cluster",
    y = "log2 enrichment (Odor / E3)"
  )

print(p_bar_subset)
