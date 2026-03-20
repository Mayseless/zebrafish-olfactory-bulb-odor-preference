suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# =========================================================
# 1) Build a cell-level table from the already mapped objects
# =========================================================
get_cell_df <- function(mapped_obj, sample_name) {
  md <- mapped_obj[[]]
  
  out <- tibble(
    sample    = if ("sample" %in% colnames(md)) as.character(md$sample) else sample_name,
    replicate = if ("replicate" %in% colnames(md)) as.character(md$replicate) else "rep1",
    cluster   = as.character(md$predicted.celltype)
  ) %>%
    filter(!is.na(sample), !is.na(replicate), !is.na(cluster))
  
  # force the sample name if needed
  out$sample <- sample_name
  out
}

cells_TCA <- get_cell_df(CaMPARI_FS_ref_TCA$mapped_query, "TCA")
cells_Cyt <- get_cell_df(CaMPARI_FS_ref_Cyt$mapped_query, "Cytidine")
cells_ATP <- get_cell_df(CaMPARI_FS_ref_ATP$mapped_query, "ATP")
cells_E3  <- get_cell_df(CaMPARI_FS_ref_E3$mapped_query,  "E3")

all_cells <- bind_rows(cells_TCA, cells_Cyt, cells_ATP, cells_E3)

# keep only samples of interest
all_cells <- all_cells %>%
  filter(sample %in% c("TCA", "Cytidine", "ATP", "E3"))

# ---------------------------------------------------------
# Define valence groups
# ---------------------------------------------------------
positive_odors <- c("TCA", "Cytidine")
negative_odors <- c("ATP")

# =========================================================
# 2) Helpers: counts -> enrichment vs E3 -> valence contrast
# =========================================================
get_cluster_counts_from_cells <- function(cell_df) {
  cell_df %>%
    group_by(sample, replicate, cluster) %>%
    summarise(n = n(), .groups = "drop")
}

compare_vs_E3_from_counts <- function(all_counts, odor_name, pseudocount = 0.5) {
  odor_counts <- all_counts %>% filter(sample == odor_name)
  e3_counts   <- all_counts %>% filter(sample == "E3")
  
  reps <- intersect(unique(odor_counts$replicate), unique(e3_counts$replicate))
  
  res_list <- lapply(reps, function(rep_i) {
    df_odor <- odor_counts %>% filter(replicate == rep_i)
    df_e3   <- e3_counts   %>% filter(replicate == rep_i)
    
    merged <- full_join(df_e3, df_odor, by = "cluster", suffix = c("_E3", "_Odor")) %>%
      mutate(
        n_E3   = replace_na(n_E3,   0L),
        n_Odor = replace_na(n_Odor, 0L)
      )
    
    N1 <- sum(merged$n_E3)
    N2 <- sum(merged$n_Odor)
    
    merged %>%
      rowwise() %>%
      mutate(
        fisher_p = {
          if (N1 == 0L || N2 == 0L) {
            NA_real_
          } else {
            mtx <- matrix(
              c(as.integer(n_Odor), as.integer(N2 - n_Odor),
                as.integer(n_E3),   as.integer(N1 - n_E3)),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("Odor", "E3"), c("in_k", "not_in_k"))
            )
            tryCatch(fisher.test(mtx, alternative = "greater")$p.value,
                     error = function(e) NA_real_)
          }
        },
        log2_ratio = {
          p_odor <- (n_Odor + pseudocount) / (N2 + 2 * pseudocount)
          p_e3   <- (n_E3   + pseudocount) / (N1 + 2 * pseudocount)
          log2(p_odor / p_e3)
        }
      ) %>%
      ungroup() %>%
      transmute(
        sample = odor_name,
        replicate = rep_i,
        cluster,
        n_E3, n_Odor,
        N_E3 = N1, N_Odor = N2,
        log2_ratio, fisher_p
      )
  })
  
  bind_rows(res_list) %>%
    group_by(sample, replicate) %>%
    mutate(fisher_q = p.adjust(fisher_p, method = "BH")) %>%
    ungroup()
}

compute_valence_contrast <- function(all_counts,
                                     positive_odors = c("TCA", "Cytidine"),
                                     negative_odors = c("ATP"),
                                     pseudocount = 0.5) {
  
  all_res <- bind_rows(
    lapply(c(positive_odors, negative_odors), function(od) {
      compare_vs_E3_from_counts(all_counts, od, pseudocount = pseudocount)
    })
  )
  
  # cluster-wise positive vs negative contrast
  contrast_df <- all_res %>%
    mutate(valence = case_when(
      sample %in% positive_odors ~ "positive",
      sample %in% negative_odors ~ "negative",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(valence)) %>%
    group_by(cluster, valence) %>%
    summarise(
      mean_log2 = mean(log2_ratio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = valence, values_from = mean_log2, values_fill = 0) %>%
    mutate(
      contrast = positive - negative,
      abs_contrast = abs(contrast)
    ) %>%
    arrange(suppressWarnings(as.numeric(cluster)), cluster)
  
  list(
    all_res = all_res,
    contrast_df = contrast_df
  )
}

# =========================================================
# 3) Observed statistic
# =========================================================
observed_counts <- get_cluster_counts_from_cells(all_cells)
obs <- compute_valence_contrast(
  observed_counts,
  positive_odors = positive_odors,
  negative_odors = negative_odors
)

observed_contrast <- obs$contrast_df

print(observed_contrast, n=100)

# =========================================================
# 4) Permutation:
#    shuffle odor labels across odor-activated cells within replicate
#    while keeping:
#      - cluster identity fixed
#      - replicate fixed
#      - number of cells per odor within replicate fixed
#      - E3 unchanged
# =========================================================
permute_odors_within_replicate <- function(cell_df,
                                           odors_to_shuffle = c("TCA", "Cytidine", "ATP"),
                                           keep_fixed = "E3") {
  
  shuffled_part <- cell_df %>%
    filter(sample %in% odors_to_shuffle) %>%
    group_by(replicate) %>%
    group_modify(~{
      dat <- .x
      dat$sample <- sample(dat$sample, size = nrow(dat), replace = FALSE)
      dat
    }) %>%
    ungroup()
  
  fixed_part <- cell_df %>%
    filter(sample %in% keep_fixed)
  
  bind_rows(shuffled_part, fixed_part)
}

# =========================================================
# 5) Run permutations
# =========================================================
run_valence_permutation_test <- function(cell_df,
                                         positive_odors = c("TCA", "Cytidine"),
                                         negative_odors = c("ATP"),
                                         odors_to_shuffle = c("TCA", "Cytidine", "ATP"),
                                         keep_fixed = "E3",
                                         pseudocount = 0.5,
                                         n_perm = 1000,
                                         seed = 123) {
  set.seed(seed)
  
  # observed
  observed_counts <- get_cluster_counts_from_cells(cell_df)
  obs <- compute_valence_contrast(
    observed_counts,
    positive_odors = positive_odors,
    negative_odors = negative_odors,
    pseudocount = pseudocount
  )
  obs_df <- obs$contrast_df %>%
    select(cluster, observed_contrast = contrast, observed_abs_contrast = abs_contrast)
  
  # permutations
  perm_list <- vector("list", n_perm)
  
  for (i in seq_len(n_perm)) {
    perm_cells <- permute_odors_within_replicate(
      cell_df,
      odors_to_shuffle = odors_to_shuffle,
      keep_fixed = keep_fixed
    )
    
    perm_counts <- get_cluster_counts_from_cells(perm_cells)
    
    perm_contrast <- compute_valence_contrast(
      perm_counts,
      positive_odors = positive_odors,
      negative_odors = negative_odors,
      pseudocount = pseudocount
    )$contrast_df %>%
      select(cluster, contrast, abs_contrast) %>%
      mutate(perm = i)
    
    perm_list[[i]] <- perm_contrast
  }
  
  perm_df <- bind_rows(perm_list)
  
  # empirical p-values
  stats_df <- obs_df %>%
    left_join(
      perm_df %>%
        group_by(cluster) %>%
        summarise(
          p_emp_two_sided = (sum(abs_contrast >= first(obs_df$observed_abs_contrast[match(cluster, obs_df$cluster)])) + 1) / (n() + 1),
          p_emp_positive  = (sum(contrast >= first(obs_df$observed_contrast[match(cluster, obs_df$cluster)])) + 1) / (n() + 1),
          p_emp_negative  = (sum(contrast <= first(obs_df$observed_contrast[match(cluster, obs_df$cluster)])) + 1) / (n() + 1),
          null_mean = mean(contrast, na.rm = TRUE),
          null_sd   = sd(contrast, na.rm = TRUE),
          null_q95_abs = quantile(abs_contrast, 0.95, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "cluster"
    ) %>%
    mutate(
      q_emp_two_sided = p.adjust(p_emp_two_sided, method = "BH"),
      q_emp_positive  = p.adjust(p_emp_positive,  method = "BH"),
      q_emp_negative  = p.adjust(p_emp_negative,  method = "BH")
    ) %>%
    arrange(suppressWarnings(as.numeric(cluster)), cluster)
  
  list(
    observed = obs_df,
    permuted = perm_df,
    stats = stats_df
  )
}

perm_res <- run_valence_permutation_test(
  cell_df = all_cells,
  positive_odors = positive_odors,
  negative_odors = negative_odors,
  odors_to_shuffle = c("TCA", "Cytidine", "ATP"),
  keep_fixed = "E3",
  pseudocount = 0.5,
  n_perm = 5000,
  seed = 123
)

# =========================================================
# 6) Results table
# =========================================================
perm_stats <- perm_res$stats %>%
  mutate(across(where(is.numeric), ~round(., 4)))

print(perm_stats, n = 100)

# optionally focus on your highlighted clusters
keep_clusters <- c("7","11","14","19","25","6","18","24")
perm_stats_subset <- perm_stats %>%
  filter(cluster %in% keep_clusters)

print(perm_stats_subset, n = Inf, width = Inf)

# =========================================================
# 7) Plot observed contrast with permutation significance
# =========================================================
plot_df <- perm_res$stats %>%
  mutate(
    cluster_num = suppressWarnings(as.numeric(cluster)),
    cluster = factor(cluster, levels = unique(cluster[order(cluster_num, cluster)])),
    sig = case_when(
      q_emp_two_sided < 0.001 ~ "***",
      q_emp_two_sided < 0.01  ~ "**",
      q_emp_two_sided < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p_perm <- ggplot(plot_df, aes(x = cluster, y = observed_contrast)) +
  geom_col(fill = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(y = observed_contrast + 0.08 * sign(observed_contrast + 1e-8),
                label = sig),
            size = 3, vjust = ifelse(plot_df$observed_contrast >= 0, 0, 1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Valence contrast per cluster after odor-label permutation",
    x = "Cluster",
    y = "Observed contrast: mean positive enrichment - mean negative enrichment"
  )

print(p_perm)

# -------------------------------------------------------------------
# Add sign-based fill variable
# -------------------------------------------------------------------
plot_df <- plot_df %>%
  mutate(fill_sign = ifelse(observed_contrast >= 0, "positive", "negative"))

sign_colors <- c("positive" = "#10eb55", "negative" = "#ff00ff")

# -------------------------------------------------------------------
# Full plot: all clusters
# -------------------------------------------------------------------
p_perm <- ggplot(plot_df, aes(x = cluster, y = observed_contrast, fill = fill_sign)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(y = observed_contrast + 0.08 * sign(observed_contrast + 1e-8),
                label = sig),
            size = 3, vjust = ifelse(plot_df$observed_contrast >= 0, 0, 1)) +
  scale_fill_manual(values = sign_colors, guide = "none") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Valence contrast per cluster after odor-label permutation",
    x     = "Cluster",
    y     = "Observed contrast: mean positive enrichment - mean negative enrichment"
  )

print(p_perm)

# -------------------------------------------------------------------
# Subset plot: selected clusters only
# -------------------------------------------------------------------
keep_clusters <- c("7","11","14","19","25","6","18","24")

plot_df_subset <- plot_df %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  mutate(cluster = factor(as.character(cluster), levels = keep_clusters))

p_perm_subset <- ggplot(plot_df_subset, aes(x = cluster, y = observed_contrast, fill = fill_sign)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(y = observed_contrast + 0.08 * sign(observed_contrast + 1e-8),
                label = sig),
            size = 3, vjust = ifelse(plot_df_subset$observed_contrast >= 0, 0, 1)) +
  scale_fill_manual(values = sign_colors, guide = "none") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Selected clusters: valence contrast after odor-label permutation",
    x     = "Cluster",
    y     = "Observed contrast: mean positive enrichment - mean negative enrichment"
  )

print(p_perm_subset)

ggsave(
  filename    = "C:/Oded_data/single_cell_seq/campari_seq/p_perm_subset.pdf",
  plot        = p_perm_subset,
  width       = 8,
  height      = 6,
  units       = "in",
  useDingbats = FALSE
)

library(ggplot2)
library(dplyr)

keep_clusters <- c("7","11","14","19","25","6","18","24")

# ---------------------------------------------------------------
# Prepare observed subset (already has fill_sign from before)
# ---------------------------------------------------------------
plot_df_subset <- plot_df %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  mutate(cluster = factor(as.character(cluster), levels = keep_clusters))

# ---------------------------------------------------------------
# Prepare permutation null subset
# ---------------------------------------------------------------
perm_subset <- perm_res$permuted %>%
  filter(as.character(cluster) %in% keep_clusters) %>%
  mutate(cluster = factor(as.character(cluster), levels = keep_clusters))

sign_colors <- c("positive" = "#10eb55", "negative" = "#ff00ff")

# ===============================================================
# PLOT 2: Violin + colored bar
# ===============================================================
p_perm_violin <- ggplot() +
  # null distribution violin
  geom_violin(data = perm_subset,
              aes(x = cluster, y = contrast),
              fill = "grey80", color = "grey50",
              alpha = 0.7, width = 0.8, linewidth = 0.3) +
  # observed bar on top
  geom_col(data = plot_df_subset,
           aes(x = cluster, y = observed_contrast, fill = fill_sign),
           width = 0.4, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = plot_df_subset,
            aes(x = cluster,
                y = observed_contrast + 0.08 * sign(observed_contrast + 1e-8),
                label = sig),
            size = 3,
            vjust = ifelse(plot_df_subset$observed_contrast >= 0, 0, 1)) +
  scale_fill_manual(values = sign_colors, guide = "none") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title    = "Selected clusters: valence contrast vs. permutation null",
    subtitle = "Grey violin = null distribution (5000 permutations); bar = observed",
    x        = "Cluster",
    y        = "Valence contrast"
  )

print(p_perm_violin)

ggsave(
  filename    = "C:/Oded_data/single_cell_seq/campari_seq/p_perm_violin_subset.pdf",
  plot        = p_perm_violin,
  width       = 8, height = 6, units = "in",
  useDingbats = FALSE
)

# =========================================================
# 8) Optional null-distribution plot for one cluster
# =========================================================
plot_null_for_cluster <- function(cluster_id, perm_res) {
  obs_val <- perm_res$observed %>%
    filter(cluster == cluster_id) %>%
    pull(observed_contrast)
  
  null_df <- perm_res$permuted %>%
    filter(cluster == cluster_id)
  
  ggplot(null_df, aes(x = contrast)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    geom_vline(xintercept = obs_val, linetype = "dashed", linewidth = 0.8) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Permutation null for cluster", cluster_id),
      x = "Permuted valence contrast",
      y = "Count"
    )
}

# examples:
print(plot_null_for_cluster("18", perm_res))
print(plot_null_for_cluster("25", perm_res))
print(plot_null_for_cluster("24", perm_res))
