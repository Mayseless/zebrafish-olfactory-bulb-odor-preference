cluster_comparison_3 <- function(df1_in, df2_in, dataset1_name, dataset2_name,
                                 # new optional knobs
                                 do_bootstrap = FALSE, B_boot = 1000, boot_seed = 42,
                                 chisq_B = 10000) {
  # deps used originally
  library(cluster)
  library(ggplot2)
  library(reshape2)
  library(mclust)
  library(clValid)
  
  ## ---------- prep ----------
  df1 <- as.data.frame(df1_in); colnames(df1) <- "cluster"
  df2 <- as.data.frame(df2_in); colnames(df2) <- "cluster"
  df1$source <- "df1"; df2$source <- "df2"
  df <- rbind(df1, df2)
  
  tbl <- table(df$source, df$cluster)  # 2 x K
  
  ## ---------- global test (same as before) ----------
  res <- chisq.test(tbl, simulate.p.value = TRUE, B = chisq_B)
  stdres <- residuals(res, type = "standardized")
  
  # your original residual-based calls
  overrep_resid  <- which(stdres["df2", ] >  2) - 1
  underrep_resid <- which(stdres["df2", ] < -2) - 1
  
  ## ---------- per-cluster proportions & ratio ----------
  cluster_counts_df1 <- table(df1$cluster)
  df1_prop <- data.frame(cluster = names(cluster_counts_df1),
                         count.x = as.numeric(cluster_counts_df1),
                         prop.x  = as.numeric(cluster_counts_df1) / sum(cluster_counts_df1),
                         stringsAsFactors = FALSE)
  
  cluster_counts_df2 <- table(df2$cluster)
  df2_prop <- data.frame(cluster = names(cluster_counts_df2),
                         count.y = as.numeric(cluster_counts_df2),
                         prop.y  = as.numeric(cluster_counts_df2) / sum(cluster_counts_df2),
                         stringsAsFactors = FALSE)
  
  prop_ratio <- merge(df1_prop, df2_prop, by = "cluster", all = TRUE)
  prop_ratio[is.na(prop_ratio)] <- 0
  prop_ratio$cluster <- as.numeric(prop_ratio$cluster)
  prop_ratio <- prop_ratio[order(prop_ratio$cluster), ]
  
  # effect sizes
  eps <- .Machine$double.eps
  prop_ratio$ratio      <- prop_ratio$prop.y / pmax(prop_ratio$prop.x, eps)
  prop_ratio$log2_ratio <- log2(pmax(prop_ratio$ratio, eps))
  
  ## ---------- per-cluster Fisher exact (one-sided) + BH ----------
  N1 <- sum(prop_ratio$count.x)
  N2 <- sum(prop_ratio$count.y)
  fish_p_greater <- fish_p_less <- rep(NA_real_, nrow(prop_ratio))
  
  for (i in seq_len(nrow(prop_ratio))) {
    in1  <- prop_ratio$count.x[i]
    in2  <- prop_ratio$count.y[i]
    mtx  <- matrix(c(in2, N2 - in2,   # row 1 = df2 (CaMPARI)
                     in1, N1 - in1),  # row 2 = df1 (10x ref)
                   nrow = 2, byrow = TRUE,
                   dimnames = list(c("df2","df1"), c("in_k","not_in_k")))
    # over-representation in df2
    fish_p_greater[i] <- tryCatch(stats::fisher.test(mtx, alternative = "greater")$p.value, error = function(e) NA_real_)
    # under-representation in df2
    fish_p_less[i]    <- tryCatch(stats::fisher.test(mtx, alternative = "less")$p.value,    error = function(e) NA_real_)
  }
  prop_ratio$fisher_p_greater <- fish_p_greater
  prop_ratio$fisher_q_greater <- p.adjust(fish_p_greater, method = "BH")
  prop_ratio$fisher_p_less    <- fish_p_less
  prop_ratio$fisher_q_less    <- p.adjust(fish_p_less,    method = "BH")
  
  # call sets based on FDR 5%
  overrep_fdr  <- prop_ratio$cluster[prop_ratio$fisher_q_greater < 0.05 & is.finite(prop_ratio$log2_ratio) & prop_ratio$log2_ratio > 0]
  underrep_fdr <- prop_ratio$cluster[prop_ratio$fisher_q_less    < 0.05 & is.finite(prop_ratio$log2_ratio) & prop_ratio$log2_ratio < 0]
  
  ## ---------- optional: bootstrap CIs for observed ratio ----------
  if (isTRUE(do_bootstrap)) {
    set.seed(boot_seed)
    vec1 <- as.character(df1$cluster)
    vec2 <- as.character(df2$cluster)
    uniq <- sort(unique(c(vec1, vec2)))
    n1 <- length(vec1); n2 <- length(vec2)
    
    boot_mat <- vapply(seq_len(B_boot), function(b) {
      s1 <- sample(vec1, n1, replace = TRUE)
      s2 <- sample(vec2, n2, replace = TRUE)
      t1 <- table(factor(s1, levels = uniq))
      t2 <- table(factor(s2, levels = uniq))
      p1 <- as.numeric(t1) / sum(t1)
      p2 <- as.numeric(t2) / sum(t2)
      r  <- p2 / pmax(p1, eps)
      r[!is.finite(r)] <- 0
      r
    }, numeric(length(uniq)))
    # boot_mat: clusters x B
    obs_q05 <- apply(boot_mat, 1, stats::quantile, probs = 0.05, na.rm = TRUE)
    obs_q50 <- apply(boot_mat, 1, stats::quantile, probs = 0.50, na.rm = TRUE)
    obs_q95 <- apply(boot_mat, 1, stats::quantile, probs = 0.95, na.rm = TRUE)
    
    # align to numeric cluster in prop_ratio
    map_idx <- match(prop_ratio$cluster, as.numeric(uniq))
    prop_ratio$obs_q05 <- obs_q05[map_idx]
    prop_ratio$obs_q50 <- obs_q50[map_idx]
    prop_ratio$obs_q95 <- obs_q95[map_idx]
  }
  
  ## ---------- plot (kept as in your original) ----------
  p1 <- ggplot(prop_ratio, aes(x = factor(cluster))) +
    geom_bar(aes(y = prop.x, fill = "Dataset 1"), stat = "identity", alpha = .4, position = "dodge") +
    geom_bar(aes(y = prop.y, fill = "Dataset 2"), stat = "identity", alpha = .4, position = "dodge") +
    ylim(0, max(prop_ratio$prop.x, prop_ratio$prop.y) * 1.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Cluster", y = "Proportion", title = "Cluster Proportions by Dataset") +
    geom_text(aes(y = prop.y * 1.05, label = ifelse(ratio > 1, round(ratio, 2), "")), 
              position = position_dodge(width = 0.9), size = 3, alpha = 0.8) +
    scale_fill_manual(
      values = c("Dataset 1" = "red", "Dataset 2" = "blue"), 
      name = "Dataset", labels = c(dataset1_name, dataset2_name)
    )
  
  ## ---------- return ----------
  result <- list(
    plot        = p1,
    overrep     = overrep_resid,          # for backward compatibility (residual z > 2)
    underrep    = underrep_resid,         # for backward compatibility
    overrep_fdr = overrep_fdr,            # NEW: Fisher+BH over-represented clusters (df2 > df1)
    underrep_fdr= underrep_fdr,           # NEW: Fisher+BH under-represented clusters (df2 < df1)
    global_p    = res$p.value,
    residuals   = stdres,
    proportions = prop_ratio              # <- enriched table used by downstream plotting
  )
  return(result)
}
