# -------------------------------------------------------------------
# 3. Pairwise odor-vs-odor comparisons within replicate
#    Tests whether subtype proportions differ between odor conditions
#    independent of baseline activation (E3)
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# Helper: pairwise Fisher test between two odor conditions within replicate
compare_odors_pairwise <- function(counts_a, counts_b,
                                   odor_a, odor_b,
                                   pseudocount = 0.5,
                                   alternative = "two.sided") {
  
  reps <- intersect(unique(counts_a$replicate), unique(counts_b$replicate))
  
  if (length(reps) == 0) {
    warning(sprintf("No shared replicates between %s and %s", odor_a, odor_b))
    return(data.frame())
  }
  
  res_list <- lapply(reps, function(rep_i) {
    df_a <- dplyr::filter(counts_a, replicate == rep_i)
    df_b <- dplyr::filter(counts_b, replicate == rep_i)
    
    merged <- dplyr::full_join(
      df_a %>% dplyr::select(cluster, n),
      df_b %>% dplyr::select(cluster, n),
      by = "cluster",
      suffix = c("_A", "_B")
    ) %>%
      dplyr::mutate(
        n_A = tidyr::replace_na(n_A, 0L),
        n_B = tidyr::replace_na(n_B, 0L)
      )
    
    N_A <- sum(merged$n_A)
    N_B <- sum(merged$n_B)
    
    out <- merged %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        fisher_p = {
          if (N_A == 0L || N_B == 0L) {
            NA_real_
          } else {
            a <- as.integer(n_A)
            b <- as.integer(N_A - n_A)
            c <- as.integer(n_B)
            d <- as.integer(N_B - n_B)
            
            a <- max(a, 0L); b <- max(b, 0L); c <- max(c, 0L); d <- max(d, 0L)
            
            mtx <- matrix(c(a, b, c, d),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c(odor_a, odor_b), c("in_k", "not_in_k")))
            
            tryCatch(
              stats::fisher.test(mtx, alternative = alternative)$p.value,
              error = function(e) NA_real_
            )
          }
        },
        odds_ratio = {
          # Haldane-Anscombe correction for stability
          ((n_A + 0.5) / (N_A - n_A + 0.5)) / ((n_B + 0.5) / (N_B - n_B + 0.5))
        },
        log2_ratio = {
          p_A <- (n_A + pseudocount) / (N_A + 2 * pseudocount)
          p_B <- (n_B + pseudocount) / (N_B + 2 * pseudocount)
          log2(p_A / p_B)
        },
        direction = dplyr::case_when(
          log2_ratio > 0 ~ odor_a,
          log2_ratio < 0 ~ odor_b,
          TRUE ~ "equal"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(
        comparison = paste0(odor_a, "_vs_", odor_b),
        sample_a   = odor_a,
        sample_b   = odor_b,
        replicate  = rep_i,
        cluster,
        n_A, n_B,
        N_A, N_B,
        prop_A = n_A / N_A,
        prop_B = n_B / N_B,
        log2_ratio,
        odds_ratio,
        fisher_p,
        direction
      )
    
    out
  })
  
  dplyr::bind_rows(res_list) %>%
    dplyr::group_by(comparison, replicate) %>%
    dplyr::mutate(fisher_q = p.adjust(fisher_p, method = "BH")) %>%
    dplyr::ungroup()
}

# Run all odor-vs-odor comparisons
pair_TCA_vs_Cyt <- compare_odors_pairwise(counts_TCA, counts_Cyt, "TCA", "Cytidine")
pair_TCA_vs_ATP <- compare_odors_pairwise(counts_TCA, counts_ATP, "TCA", "ATP")
pair_Cyt_vs_ATP <- compare_odors_pairwise(counts_Cyt, counts_ATP, "Cytidine", "ATP")

pairwise_odor_enrichment <- bind_rows(
  pair_TCA_vs_Cyt,
  pair_TCA_vs_ATP,
  pair_Cyt_vs_ATP
)
# -------------------------------------------------------------------
# 4. Summary across replicates
# -------------------------------------------------------------------
pairwise_summary <- pairwise_odor_enrichment %>%
  group_by(comparison, cluster) %>%
  summarise(
    mean_log2 = mean(log2_ratio, na.rm = TRUE),
    sd_log2   = sd(log2_ratio, na.rm = TRUE),
    mean_or   = mean(odds_ratio, na.rm = TRUE),
    n_reps    = sum(!is.na(log2_ratio)),
    min_q     = suppressWarnings(min(fisher_q, na.rm = TRUE)),
    n_sig_005 = sum(fisher_q < 0.05, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    se_log2 = ifelse(n_reps > 1 & is.finite(sd_log2), sd_log2 / sqrt(n_reps), 0),
    sig_label = case_when(
      is.finite(min_q) & min_q < 0.001 ~ "***",
      is.finite(min_q) & min_q < 0.01  ~ "**",
      is.finite(min_q) & min_q < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Order clusters numerically when possible
cluster_levels_pair <- pairwise_summary %>%
  mutate(cluster_num = suppressWarnings(as.numeric(as.character(cluster)))) %>%
  distinct(cluster, cluster_num) %>%
  mutate(cluster_num = ifelse(is.na(cluster_num), Inf, cluster_num)) %>%
  arrange(cluster_num, cluster) %>%
  pull(cluster)

pairwise_summary <- pairwise_summary %>%
  mutate(cluster = factor(cluster, levels = cluster_levels_pair))

pairwise_odor_enrichment <- pairwise_odor_enrichment %>%
  mutate(cluster = factor(cluster, levels = cluster_levels_pair))

# -------------------------------------------------------------------
# 5. Print results table
# -------------------------------------------------------------------
cat("\n================ Pairwise odor-vs-odor subtype comparisons ================\n")
pairwise_summary %>%
  arrange(comparison, cluster) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print(n = Inf,width=Inf)

# Optional: show only clusters significant in at least one comparison
cat("\n================ Significant pairwise differences (min q < 0.05) ================\n")
pairwise_summary %>%
  filter(is.finite(min_q), min_q < 0.05) %>%
  arrange(comparison, min_q, desc(abs(mean_log2))) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print(n = Inf,width=Inf)

# -------------------------------------------------------------------
# 6. Plot: pooled pairwise comparison across all clusters
# -------------------------------------------------------------------
comparison_colors <- c(
  "TCA_vs_Cytidine" = "#7d8388",
  "TCA_vs_ATP"      = "#10eb55",
  "Cytidine_vs_ATP" = "#075707"
)

pd <- position_dodge2(width = 0.8, preserve = "single")

p_pairwise_bar <- ggplot(pairwise_summary,
                         aes(x = cluster, y = mean_log2, fill = comparison)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = comparison_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Pairwise odor-vs-odor subtype recruitment",
    subtitle = "Positive values indicate enrichment in the first odor of each comparison",
    x = "Cluster",
    y = "log2 enrichment (odor A / odor B)",
    fill = "Comparison"
  )

print(p_pairwise_bar)

# -------------------------------------------------------------------
# 7. Plot: selected clusters only
# -------------------------------------------------------------------
keep_clusters <- c("7","11","14","19","25","6","18","24")

pairwise_summary_subset <- pairwise_summary %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  mutate(cluster = factor(as.character(cluster), levels = keep_clusters))

p_pairwise_subset <- ggplot(pairwise_summary_subset,
                            aes(x = cluster, y = mean_log2, fill = comparison)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = comparison_colors, drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Selected clusters: pairwise odor-vs-odor recruitment",
    subtitle = "Positive values indicate enrichment in the first odor of each comparison",
    x = "Cluster",
    y = "log2 enrichment (odor A / odor B)",
    fill = "Comparison"
  )

print(p_pairwise_subset)


# -------------------------------------------------------------------
# Helper: add a sign-aware fill variable
# -------------------------------------------------------------------
pairwise_summary <- pairwise_summary %>%
  mutate(
    fill_var = case_when(
      comparison == "TCA_vs_Cytidine" & mean_log2 >= 0 ~ "TCA_pos",
      comparison == "TCA_vs_Cytidine" & mean_log2 <  0 ~ "TCA_neg",
      comparison == "Cytidine_vs_ATP" & mean_log2 >= 0 ~ "Cyt_pos",
      comparison == "Cytidine_vs_ATP" & mean_log2 <  0 ~ "Cyt_neg",
      comparison == "TCA_vs_ATP"      & mean_log2 >= 0 ~ "ATP_pos",
      comparison == "TCA_vs_ATP"      & mean_log2 <  0 ~ "ATP_neg"
    ),
    fill_var = factor(fill_var,
                      levels = c("TCA_pos","TCA_neg",
                                 "Cyt_pos","Cyt_neg",
                                 "ATP_pos","ATP_neg"))
  )

comparison_colors <- c(
  "TCA_pos" = "#10eb55",   # green      (TCA positive)
  "TCA_neg" = "#ff00ff",   # magenta    (TCA negative)
  "Cyt_pos" = "#075707",   # dark green (Cytidine positive)
  "Cyt_neg" = "#ff00ff",   # magenta    (Cytidine negative)
  "ATP_pos" = "#7d8388",   # grey       (ATP positive, unchanged)
  "ATP_neg" = "#7d8388"    # grey       (ATP negative, unchanged — adjust if needed)
)

# Legend labels so the fill_var codes don't show in the legend
comparison_labels <- c(
  "TCA_pos" = "TCA_vs_Cytidine",
  "TCA_neg" = "TCA_vs_Cytidine (neg)",
  "Cyt_pos" = "Cytidine_vs_ATP",
  "Cyt_neg" = "Cytidine_vs_ATP (neg)",
  "ATP_pos" = "TCA_vs_ATP",
  "ATP_neg" = "TCA_vs_ATP (neg)"
)

pd <- position_dodge2(width = 0.8, preserve = "single")

# -------------------------------------------------------------------
# 6. Plot: pooled pairwise comparison across all clusters
# -------------------------------------------------------------------
p_pairwise_bar <- ggplot(pairwise_summary,
                         aes(x = cluster, y = mean_log2, fill = fill_var)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = comparison_colors,
                    labels = comparison_labels,
                    drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title    = "Pairwise odor-vs-odor subtype recruitment",
    subtitle = "Positive values indicate enrichment in the first odor of each comparison",
    x        = "Cluster",
    y        = "log2 enrichment (odor A / odor B)",
    fill     = "Comparison"
  )

print(p_pairwise_bar)

# -------------------------------------------------------------------
# 7. Plot: selected clusters only
# -------------------------------------------------------------------
keep_clusters <- c("7","11","14","19","25","6","18","24")

pairwise_summary_subset <- pairwise_summary %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  mutate(cluster = factor(as.character(cluster), levels = keep_clusters))

p_pairwise_subset <- ggplot(pairwise_summary_subset,
                            aes(x = cluster, y = mean_log2, fill = fill_var)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                position = pd, width = 0.4, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = mean_log2 + se_log2 + 0.06),
            position = pd, size = 3, vjust = 0) +
  scale_fill_manual(values = comparison_colors,
                    labels = comparison_labels,
                    drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Selected clusters: pairwise odor-vs-odor recruitment",
    subtitle = "Positive values indicate enrichment in the first odor of each comparison",
    x        = "Cluster",
    y        = "log2 enrichment (odor A / odor B)",
    fill     = "Comparison"
  )

print(p_pairwise_subset)
ggsave(
  filename = "C:/Oded_data/single_cell_seq/campari_seq/p_pairwise_subset.pdf",
  plot     = p_pairwise_subset,
  width    = 8,
  height   = 6,
  units    = "in",
  useDingbats = FALSE   # keeps text/shapes as editable vectors in Illustrator / Inkscape
)
