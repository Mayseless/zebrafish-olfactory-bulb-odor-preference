
"""

Script to:
1) Load a precomputed cross-correlation matrix (CSV).
2) Perform hierarchical clustering with Seaborn's clustermap.
3) Assign cluster labels and map them to known odor preference scores.
4) Use permutation tests to assess significant clusters.
5) Perform additional analyses (odor consistency, chemical class testing,
   within-odor correlations, and more).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import fcluster
from scipy.stats import chi2_contingency
import joypy
from scipy.cluster.hierarchy import linkage
from collections import defaultdict


#%%
# ------------------------------
# 1. Define Preferences & Helper Functions
# ------------------------------

#  dictionary of median preference scores from behavior assay 
preference_scores = {
    "5bCholicAcid100uM": 0.487500,
    "AAmix100uM": 0.535000,
    "AMP100uM": 0.340269,
    "ATP100uM": 0.332500,
    "ATP10uM": 0.385000,
    "ATP500uM": 0.270833,
    "Adenosine100uM": 0.476667,
    "Ala100uM": 0.461250,
    "Ala10uM": 0.429167,
    "Ala1mM": 0.507500,
    "AmCl100uM": 0.604167,
    "Arg100uM": 0.570833,
    "Arg10uM": 0.428261,
    "Arg500uM": 0.616667,
    "BAmix100uM": 0.589583,
    "Cad100uM": 0.512949,
    "Cad1mM": 0.379167,
    "Cad300uM": 0.463314,
    "Cichlid": 0.700000,
    "CitricAcid100uM": 0.378500,
    "CitricAcid500uM": 0.188012,
    "DGlu100uM": 0.484167,
    "E3": 0.541250,
    "GCA300uM": 0.507500,
    "Glycine100uM": 0.595000,
    "IMP100uM": 0.514167,
    "ITP100uM": 0.476667,
    "LCys100uM": 0.570000,
    "LCys500uM": 0.551887,
    "LPro100uM": 0.567500,
    "LSer100uM": 0.497917,
    "Lmet100uM": 0.413333,
    "NaCl1mM": 0.554167,
    "NaCl200mM": 0.521997,
    "NaCl50mM": 0.467500,
    "phenylalanine100uM": 0.400833,
    "Spont": 0.390338,
    "TCA100uM": 0.566068,
    "TCA1mM": 0.568750,
    "TCA300uM": 0.610833,
    "TDCA300uM": 0.520000,
    "Thymine100uM": 0.560833,
    "Trp100uM": 0.496041,
    "adult": 0.704754,
    "Cytidine100uM": 0.695000,
    "food": 0.633333,
    "urea100uM": 0.501250
}

chemical_classes = {
    "Amino Acids": ["ala", "arg", "gly", "lcys", "lmet", "lpro", "lser", "phe", "trp", "phenylalanine","LCys100uM"],
    "Bile Acids": ["tca", "gca", "tdca"],
    "Nucleic Acids": ["atp", "amp", "adenosine", "itp", "imp", "cytidine", "thymine"],
    "Salts/Inorganics": ["amcl", "urea", "nacl"],
    "Amines": ["cadaverine", "cad"],
    "Complex Odors": ["aamix","bamix","cichlid", "adult", "food"],
    "Organic Acids": ["citric"],
    "Control": ["e3", "spont", "spontaneous"]
}



def get_preference_score(label, pref_dict):
    """
    Matches the provided label to keys in pref_dict using case-insensitive substring matching.
    Also accounts for a common minor spelling difference by replacing 'anic' with 'ic'.
    Returns the score if a match is found; otherwise, returns None.
    """
    label_clean = label.split('60sec')[0].strip().lower()
    for key, score in pref_dict.items():
        key_lower = key.lower()
        # Direct substring match in either direction
        if key_lower in label_clean or label_clean in key_lower:
            return score
        # Additional check for "anic" vs "ic" discrepancy
        if key_lower.replace("anic", "ic") in label_clean or label_clean.replace("anic", "ic") in key_lower:
            return score
    return None

def get_chemical_class(odor):
    odor_lc = odor.lower()
    for chem_class, substrings in chemical_classes.items():
        for sub in substrings:
            if sub.lower() in odor_lc:  # case-insensitive matching
                return chem_class
    return "Unmapped"




def permutation_test(cluster_scores, all_scores, n_permutations=10000, seed=42, significance_alpha=0.05):
    """
    For a given cluster (list of preference scores), randomly sample the same number of scores
    from the overall set repeatedly to build a null distribution.

    - If observed_mean >= overall_mean, the p-value is the fraction of permuted means >= observed_mean.
      The preliminary label is "attractive".
    - If observed_mean <  overall_mean, we do the symmetric comparison and label it "aversive".
    - If p_value >= significance_alpha, we override to "neutral".
    """
    np.random.seed(seed)
    observed_mean = np.mean(cluster_scores)
    cluster_size = len(cluster_scores)
    permuted_means = []

    for _ in range(n_permutations):
        sample = np.random.choice(all_scores, size=cluster_size, replace=False)
        permuted_means.append(np.mean(sample))

    permuted_means = np.array(permuted_means)
    overall_mean = np.mean(all_scores)

    # Determine direction first
    if observed_mean >= overall_mean:
        p_value = np.mean(permuted_means >= observed_mean)
        label = "attractive"
    else:
        p_value = np.mean(permuted_means <= observed_mean)
        label = "aversive"

    # Apply significance check
    if p_value >= significance_alpha:
        label = "neutral"

    return observed_mean, p_value, label, permuted_means

def compute_class_pattern_association(
    odor_df: pd.DataFrame,
    contig_od: pd.DataFrame | None = None,
    outdir: str | None = None,
    alpha: float = 0.05,
    min_count: int = 1,
):
    """
    Compute which chemical classes are over/under-represented in each odor *cluster*
    at the ODORANT level.

    Parameters
    ----------
    odor_df : DataFrame with columns ["odorant", "cluster", "chemical_class"]
    contig_od : optional contingency (rows=cluster, cols=chemical_class). If None, build from odor_df.
    outdir : optional folder to write CSVs
    alpha : FDR threshold
    min_count : cells with observed < min_count are flagged but still shown

    Returns
    -------
    assoc_long : long-form DataFrame with observed/expected, standardized residuals (z),
                 raw and FDR-adjusted p-values, direction (enriched/depleted), significance flag
    assoc_summary : per-cluster summary of significant enrichments/depletions
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2_contingency, norm
    from statsmodels.stats.multitest import multipletests
    from pathlib import Path

    # 1) Build contingency if needed
    if contig_od is None:
        if not {"cluster", "chemical_class"}.issubset(odor_df.columns):
            raise ValueError("odor_df must contain columns: 'cluster' and 'chemical_class'")
        contig_od = pd.crosstab(odor_df["cluster"], odor_df["chemical_class"])

    # 2) Global chi-square + expected counts
    chi2, p_global, dof, expected = chi2_contingency(contig_od, correction=False)
    obs = contig_od.values.astype(float)
    exp = expected.astype(float)

    # 3) Standardized residuals (Haberman)
    #    r_ij = (O - E) / sqrt(E * (1 - row_prop_i) * (1 - col_prop_j))
    N = obs.sum()
    row_sums = obs.sum(axis=1, keepdims=True)
    col_sums = obs.sum(axis=0, keepdims=True)
    row_prop = row_sums / N
    col_prop = col_sums / N
    denom = np.sqrt(exp * (1.0 - row_prop) * (1.0 - col_prop))
    with np.errstate(divide="ignore", invalid="ignore"):
        std_resid = (obs - exp) / denom
    std_resid = np.where(np.isfinite(std_resid), std_resid, np.nan)

    # 4) Per-cell p-values via normal approximation
    z = std_resid
    p_cell = 2.0 * (1.0 - norm.cdf(np.abs(z)))

    # 5) Assemble long-form table
    long = []
    row_names = contig_od.index.tolist()
    col_names = contig_od.columns.tolist()
    for i, rname in enumerate(row_names):
        for j, cname in enumerate(col_names):
            long.append({
                "cluster": rname,
                "chemical_class": cname,
                "observed": obs[i, j],
                "expected": exp[i, j],
                "std_resid": z[i, j],
                "p_cell_raw": p_cell[i, j],
            })
    assoc_long = pd.DataFrame(long)

    # 6) Multiple testing (FDR BH) over all cells
    mask = assoc_long["p_cell_raw"].notna()
    rej, qvals, *_ = multipletests(assoc_long.loc[mask, "p_cell_raw"], method="fdr_bh", alpha=alpha)
    assoc_long.loc[mask, "p_cell_fdr"] = qvals
    assoc_long["p_cell_fdr"] = assoc_long["p_cell_fdr"].astype(float)

    # 7) Flags and direction
    assoc_long["significant"] = assoc_long["p_cell_fdr"] < alpha
    assoc_long["direction"] = np.where(assoc_long["std_resid"] > 0, "enriched",
                                np.where(assoc_long["std_resid"] < 0, "depleted", "none"))
    assoc_long.loc[assoc_long["observed"] < min_count, ["significant"]] = False  # optional guard

    # 8) Per-cluster summary (top hits)
    def _summarize(df):
        sig = df[df["significant"]].copy()
        if sig.empty:
            return pd.Series({"n_sig": 0, "top_enriched": "", "top_depleted": ""})
        top_enr = (sig[sig["direction"] == "enriched"]
                   .sort_values("std_resid", ascending=False)
                   .head(3)["chemical_class"].tolist())
        top_dep = (sig[sig["direction"] == "depleted"]
                   .sort_values("std_resid", ascending=True)
                   .head(3)["chemical_class"].tolist())
        return pd.Series({
            "n_sig": int(sig.shape[0]),
            "top_enriched": ", ".join(top_enr) if top_enr else "",
            "top_depleted": ", ".join(top_dep) if top_dep else "",
        })

    assoc_summary = assoc_long.groupby("cluster").apply(_summarize).reset_index()

    # 9) Save (optional)
    if outdir is not None:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        assoc_long.sort_values(["cluster", "p_cell_fdr", "chemical_class"]).to_csv(
            out / "ChemClass_Association_long.csv", index=False
        )
        assoc_summary.to_csv(out / "ChemClass_Association_summary.csv", index=False)
        # Also spill the global stats
        with open(out / "ChemClass_Association_global.txt", "w", encoding="utf-8") as fh:
            fh.write(f"Chi2={chi2:.4f}, dof={dof}, p_global={p_global:.4g}, N={int(N)}\n")

    return assoc_long, assoc_summary


#%%
output_folder = "./campari-data-analysis-master/outputPlots_d10/"
#os.makedirs(output_folder, exist_ok=True)

STORE_FILENAME = "./campari-data-analysis-master/data/data.h5"
OUT_FILENAME = "./campari-data-analysis-master/out/3_output_dmax10.h5"
FILTER_SETTINGS = "dmin_0__dmax_10__cellvol_100_1000__topp_10"
BW = 5.0
method = "average" # ['complete', 'ward','single','weighted']
# ["Dch_matrix", "S_matrix", "XCM_matrix", "XC_matrix"]

matrix_str = "XCM_matrix"  # e.g., "XCM_matrix", "XC_matrix", etc.

#%%
# ---- Quick study-size summary: odors & larvae ----
import re

def _base_odor(name: str) -> str:
    """Collapse concentration suffixes to a base odor label."""
    s = name.replace("µ", "u")  # normalize micro symbol
    s = re.sub(r'[\s_@-]*\d+(\.\d+)?\s*(?:uM|um|mM|nM|pM|%)\s*$', '', s, flags=re.IGNORECASE)
    s = re.sub(r'[\s_@-]*(?:low|high|med|medium)\s*$', '', s, flags=re.IGNORECASE)
    return s.strip('_- ').strip()

row_info["larva"] = row_info["full_label"].str.split(",", n=1, expand=True)[0].str.strip()
row_info["base_odor"] = row_info["odor"].map(_base_odor)

n_odor_conditions = row_info["odor"].nunique()
n_base_odors      = row_info["base_odor"].nunique()
n_larvae          = row_info["larva"].nunique()

print("\n[SUMMARY]")
print(f"  Unique odor conditions (incl. concentrations): {n_odor_conditions}")
print(f"  Unique base odors (concentrations collapsed):  {n_base_odors}")
print(f"  Unique larvae:                                  {n_larvae}")

# ---- List all odors and concentrations grouped by base odor ----
print("\n[ODORS TESTED]")
for base, sub in row_info.groupby("base_odor"):
    conds = sorted(sub["odor"].unique())
    if len(conds) == 1:
        print(f"  {base}: {conds[0]}")
    else:
        print(f"  {base}: {', '.join(conds)}")
# ---- Animals per odor ----
animals_per_odor = (
    row_info.groupby("odor")["larva"]
    .nunique()
    .sort_values(ascending=False)
)

print("\n[ANIMALS PER ODOR]")
for odor, n_animals in animals_per_odor.items():
    print(f"  {odor}: {n_animals} larvae")

# also append to summary file
with open(os.path.join(output_folder, "SummaryCounts.txt"), "a", encoding="utf-8") as fh:
    fh.write("\n[ANIMALS PER ODOR]\n")
    for odor, n_animals in animals_per_odor.items():
        fh.write(f"{odor}: {n_animals} larvae\n")

# also save alongside the summary
with open(os.path.join(output_folder, "SummaryCounts.txt"), "w", encoding="utf-8") as fh:
    fh.write("[SUMMARY]\n")
    fh.write(f"Unique odor conditions (incl. concentrations): {n_odor_conditions}\n")
    fh.write(f"Unique base odors (concentrations collapsed):  {n_base_odors}\n")
    fh.write(f"Unique larvae:                                  {n_larvae}\n\n")
    fh.write("[ODORS TESTED]\n")
    for base, sub in row_info.groupby("base_odor"):
        conds = sorted(sub["odor"].unique())
        if len(conds) == 1:
            fh.write(f"{base}: {conds[0]}\n")
        else:
            fh.write(f"{base}: {', '.join(conds)}\n")


#%%

# --------------------
# Main Script 
# --------------------

if __name__ == "__main__":
    matrix_list = ["XCM_matrix"]

    # ------------- USER CONFIG -------------
    stab_cut = 0.50            # clusters below this bootstrap stability are flagged
    filter_unstable = False    # hide clusters with stability < stab_cut if True
    always_plot_heatmap = True # always export the co‑cluster heat‑map
    n_bootstraps = 1000         # bootstrap replicates
    palette_name = "tab20"      # colour palette for cluster bars
    dendro_ratio = (0.08, 0.08)  # (row, col) proportion for dendrograms (0–1)
    colors_ratio = (0.03, 0.03)  # (row, col) thickness of colour bars relative to heat‑map
    cbar_pos = (0.02, 0.82, 0.03, 0.15)  # x, y, width, height of colour‑bar

    fig_size_main = (13, 13)     # base clustermap size
    fig_size_high = (15, 15)     # coloured clustermap size
    # ---------------------------------------

    for matrix_str in matrix_list:
        print("\n=== Computing Clustermap with matrix =", matrix_str, "===")

        cluster_methods = ["average"]
        threshold_multipliers = [ 0.6]

        # --- 1) Load the Matrix ---
        store = pd.HDFStore(OUT_FILENAME)
        matrix = store.get(f"/{FILTER_SETTINGS}/{matrix_str}")
        store.close()
        print(f"Loaded matrix: {matrix_str} from {OUT_FILENAME}")
        print("Matrix shape:", matrix.shape)

        # Simplify row/col names
        def simplify_index(idx):
            parts = idx.split(',')
            return f"{parts[1].strip()}_{parts[0].strip()}"

        matrix_renamed = matrix.copy()
        matrix_renamed.index = [simplify_index(i) for i in matrix_renamed.index]
        matrix_renamed.columns = [simplify_index(c) for c in matrix_renamed.columns]

        matrix_dist = matrix_renamed   # 1‑corr if needed

        row_info = pd.DataFrame({
            "full_label": list(matrix.index),
            "odor": [x.split(',')[1] for x in matrix.index],
        })

        # ----------------------------
        # 2) iterate over linkage methods
        # ----------------------------
        for method in cluster_methods:
            print("\n=== Computing Clustermap with method =", method, "===")

            g = sns.clustermap(matrix_dist, metric="correlation", method=method,
                               cmap="inferno", figsize=(12, 12), xticklabels=True, yticklabels=True)
            plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize=4, rotation=45, ha="right")
            plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=4)
            plt.tight_layout()
            png_path = os.path.join(output_folder, f"{matrix_str}_{method}_clustermap.png")
            pdf_path = os.path.join(output_folder, f"{matrix_str}_{method}_clustermap.pdf")
            g.fig.savefig(png_path, dpi=600, bbox_inches="tight")
            g.fig.savefig(pdf_path, bbox_inches="tight")
            plt.show(g.fig)

            linkage_matrix = g.dendrogram_row.linkage

            # ----------------------------
            # thresholds
            # ----------------------------
            
            from sklearn.manifold import MDS
            from sklearn.metrics import silhouette_samples, silhouette_score
            
            for thresh_mult in threshold_multipliers:
                print(f"\n--- method={method}, threshold multiplier={thresh_mult} ---")
                threshold = thresh_mult * linkage_matrix[:, 2].max()
                cluster_labels = fcluster(linkage_matrix, t=threshold, criterion="distance")
            
                assignments = pd.DataFrame({
                    "odor": row_info["odor"],
                    "cluster": cluster_labels,
                })
                assignments["pref_score"] = assignments["odor"].apply(lambda o: get_preference_score(o, preference_scores))
                assignments_valid = assignments.dropna(subset=["pref_score"]).copy()
            
                # ----------------------------
                # Cluster-level permutation tests (valence)
                # ----------------------------
                all_pref = assignments_valid.pref_score.values
                overall_mean, overall_sd = all_pref.mean(), all_pref.std(ddof=1)
            
                results = {}
                for cid, grp in assignments_valid.groupby("cluster"):
                    sc = grp.pref_score.values
                    obs_mean, p_val, vlabel, perm_means = permutation_test(sc, all_pref)
                    results[cid] = {
                        "obs_mean": obs_mean,
                        "obs_z": (obs_mean - overall_mean) / overall_sd if overall_sd > 0 else np.nan,
                        "p": p_val,
                        "valence": vlabel,
                        "perm_z": (perm_means - overall_mean) / overall_sd if overall_sd > 0 else np.full_like(perm_means, np.nan),
                    }
            
                from statsmodels.stats.multitest import multipletests
                pvals = np.array([results[c]["p"] for c in sorted(results)])
                rej, qvals, *_ = multipletests(pvals, method="fdr_bh")
                q_by_cluster = {cid: q for cid, q in zip(sorted(results), qvals)}
            
                for cid in sorted(results):
                    res = results[cid]
                    print(
                        f"Cluster {cid}: mean={res['obs_mean']:.3f} "
                        f"Z={res['obs_z']:.2f} p={res['p']:.3g} → {res['valence']} , qval={q_by_cluster[cid]:.3g}"
                    )
            
                # ----------------------------
                # Define ODORANT-level labels early (dose-collapsed) so we can remove singletons
                # ----------------------------
                # Chemical class uses get_chemical_class (you already have it)
                assignments_valid["chemical_class"] = assignments_valid["odor"].apply(get_chemical_class)
            
                # Strip concentrations/units: ATP100uM -> ATP ; Cad1mM -> Cad ; tolerate spaces
                assignments_valid["odorant"] = (
                    assignments_valid["odor"]
                    .str.replace(r"(\d+(?:\.\d+)?)\s*(?:uM|µM|mM)", "", regex=True, case=False)
                    .str.replace(r"\s+", "", regex=True)
                )
            
                # Majority vote per odorant across trials/doses
                odorant2cluster = (
                    assignments_valid
                    .groupby(["odorant", "cluster"]).size().reset_index(name="n")
                    .sort_values(["odorant", "n"], ascending=[True, False])
                    .drop_duplicates(subset=["odorant"])
                    .set_index("odorant")["cluster"]
                )
            
                odorant_df = (
                    pd.DataFrame({
                        "cluster": odorant2cluster,
                        "chemical_class": odorant2cluster.index.map(lambda od: get_chemical_class(od))
                    })
                    .reset_index(names="odorant")
                )
            
                # ----------------------------
                # REMOVE SINGLETON CLUSTERS (clusters with <2 trials)
                # ----------------------------
                # Trial-level counts (rows/cols in the clustermap)
                trial_counts = pd.Series(cluster_labels).value_counts().sort_index()
                
                singleton_clusters = trial_counts[trial_counts < 2].index.tolist()
                keep_clusters = sorted([c for c in trial_counts.index if c not in singleton_clusters])
                
                print("[CLUSTER FILTER | TRIAL-LEVEL]")
                print(trial_counts.to_string())
                print(f"Dropping (n_trials < 2): {singleton_clusters}")
                print(f"Keeping: {keep_clusters}")
            
                # Filter results to kept clusters
                results_kept = {cid: results[cid] for cid in keep_clusters if cid in results}
                plot_clusters = sorted(results_kept.keys())
            
                # If you want to optionally drop unstable clusters too, do it here (AFTER singleton removal)
                # stability is computed later, so we keep this as a post-stability option (see below)
            
                # ----------------------------
                # Bootstrapping co-cluster stability (same as your logic, but we will *report* kept clusters)
                # ----------------------------
                print("Bootstrapping cluster stability…")
                n_items = matrix_dist.shape[0]
                take = np.zeros((n_items, n_items), dtype=int)
                same = np.zeros((n_items, n_items), dtype=int)
            
                for _ in range(n_bootstraps):
                    idx = np.random.choice(n_items, n_items, replace=True)
                    uniq = np.unique(idx)
                    sub_d = matrix_dist.iloc[uniq, uniq]
                    boot_link = linkage(sub_d, method=method, metric="correlation")
                    boot_lab = fcluster(boot_link, t=threshold, criterion="distance")
                    for a, ia in enumerate(uniq):
                        for b, ib in enumerate(uniq):
                            if b <= a:
                                continue
                            take[ia, ib] += 1; take[ib, ia] += 1
                            if boot_lab[a] == boot_lab[b]:
                                same[ia, ib] += 1; same[ib, ia] += 1
            
                cocluster_prob = np.divide(
                    same, take,
                    out=np.full_like(same, np.nan, dtype=float),
                    where=take != 0
                )
            
                # Stability per cluster (mean within-cluster co-cluster prob)
                stability = {}
                for cid in sorted(np.unique(cluster_labels)):
                    mem = np.where(cluster_labels == cid)[0]
                    if mem.size < 2:
                        stability[cid] = np.nan
                    else:
                        sub = cocluster_prob[np.ix_(mem, mem)]
                        stability[cid] = np.nanmean(sub[np.triu_indices_from(sub, 1)])
            
                # Write stability report (kept + dropped)
                out_stab = os.path.join(output_folder, f"ClusterStability_{method}_th{thresh_mult}.txt")
                print("\nBootstrap stability (mean co-cluster probability)")
                with open(out_stab, "w", encoding="utf-8") as fh:
                    fh.write("[KEPT CLUSTERS: n_odorants>=2]\n")
                    for cid in keep_clusters:
                        s = stability.get(cid, float("nan"))
                        line = f"  Cluster {cid:>2}:  {s:0.3f}"
                        print(line)
                        fh.write(line + "\n")
            
                    fh.write("\n[DROPPED SINGLETON CLUSTERS: n_odorants<2]\n")
                    for cid in singleton_clusters:
                        s = stability.get(cid, float("nan"))
                        line = f"  Cluster {cid:>2}:  {s:0.3f}"
                        print(line)
                        fh.write(line + "\n")
            
                # Optional: filter unstable clusters (AFTER stability exists)
                if filter_unstable:
                    plot_clusters = [cid for cid in plot_clusters if (stability.get(cid, np.nan) >= stab_cut)]
                    print(f"[FILTER] After stability cutoff ({stab_cut}): keeping clusters {plot_clusters}")
                    results_kept = {cid: results_kept[cid] for cid in plot_clusters}
            
                # ----------------------------
                # Colors for KEPT clusters only
                # ----------------------------
                if len(plot_clusters) == 0:
                    print("[WARN] No clusters left after singleton (and optional stability) filtering. Skipping plots/tests.")
                    continue
            
                pal = sns.color_palette(palette_name, len(plot_clusters))
                colmap = {cid: pal[i % len(pal)] for i, cid in enumerate(plot_clusters)}
            
                def edge_color(cid):
                    return {"attractive": "green", "aversive": "magenta", "neutral": "gray"}[results_kept[cid]["valence"]]
            
                # ----------------------------
                # Filter matrix to KEPT clusters (trial-level)
                # ----------------------------
                keep_mask = np.isin(cluster_labels, plot_clusters)
                matrix_filt = matrix_dist.loc[keep_mask, keep_mask]
            
                # ----------------------------
                # Filtered colorful clustermap + legend (KEPT clusters only)
                # ----------------------------
                if matrix_filt.shape[0] >= 2:
                    row_colors = [colmap[c] for c in cluster_labels[keep_mask]]
            
                    g_all = sns.clustermap(
                        matrix_filt,
                        metric="correlation",
                        method=method,
                        row_colors=row_colors,
                        col_colors=row_colors,
                        cmap="inferno",
                        figsize=fig_size_high,
                        xticklabels=True,
                        yticklabels=True,
                        dendrogram_ratio=dendro_ratio,
                        colors_ratio=colors_ratio,
                        cbar_pos=cbar_pos,
                    )
            
                    plt.setp(g_all.ax_heatmap.xaxis.get_majorticklabels(), fontsize=3, rotation=45, ha="right")
                    plt.setp(g_all.ax_heatmap.yaxis.get_majorticklabels(), fontsize=3)
            
                    g_all.ax_heatmap.set_title(f"Clusters (n_odorants>=2) | thr={thresh_mult}", fontsize=10)
                    handles = [plt.matplotlib.patches.Patch(color=colmap[c], label=f"Cluster {c}") for c in plot_clusters]
                    g_all.ax_heatmap.legend(
                        handles=handles,
                        title="Clusters",
                        bbox_to_anchor=(1.1, 0.5),
                        loc="center left",
                        frameon=False,
                        fontsize=6,
                        title_fontsize=6,
                    )
            
                    plt.tight_layout()
                    base = f"Clustermap_NoSingleton_{method}_th{thresh_mult}"
                    g_all.fig.savefig(os.path.join(output_folder, f"{base}.png"), dpi=600, bbox_inches="tight")
                    g_all.fig.savefig(os.path.join(output_folder, f"{base}.pdf"), bbox_inches="tight")
                    plt.show(g_all.fig)
                else:
                    print("[WARN] Not enough trials left for filtered clustermap.")
            
                # ----------------------------
                # Odors per cluster text output (fixed formatting) — KEPT clusters only
                # ----------------------------
                out_od = os.path.join(output_folder, f"Odors_perCluster_NoSingleton_{method}_th{thresh_mult}.txt")
                with open(out_od, "w", encoding="utf-8") as f:
                    for cid in plot_clusters:
                        odrs = assignments.loc[cluster_labels == cid, "odor"].unique()
                        f.write(f"Cluster {cid} (n={len(odrs)}):\n")
                        for od in odrs:
                            f.write(f"  - {od}\n")
                        f.write("\n")
            
                # ----------------------------
                # Co-cluster heatmap (filtered)
                # ----------------------------
                if always_plot_heatmap:
                    idx_kept = np.where(keep_mask)[0]
                    if idx_kept.size >= 2:
                        cocluster_filt = cocluster_prob[np.ix_(idx_kept, idx_kept)]
                        plt.figure(figsize=(6, 5))
                        sns.heatmap(cocluster_filt, cmap="viridis", vmin=0, vmax=1)
                        plt.title(f"Co-cluster probability (filtered) | thr={thresh_mult}")
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_folder, f"CoclusterProb_NoSingleton_{method}_th{thresh_mult}.png"), dpi=300)
                        plt.show()
                    else:
                        print("[WARN] Not enough items left for cocluster heatmap.")
            
                # ----------------------------
                # Cluster separation visualization 1: MDS embedding (2D), distance = 1 - r
                # ----------------------------
                idx_kept = np.where(keep_mask)[0]
                if idx_kept.size >= 3:
                    D_full = 1.0 - matrix_dist.values
                    D = D_full[np.ix_(idx_kept, idx_kept)]
                    labs = cluster_labels[idx_kept]
            
                    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0, normalized_stress="auto")
                    XY = mds.fit_transform(D)
            
                    plt.figure(figsize=(6.5, 5.5))
                    for cid in plot_clusters:
                        sel = (labs == cid)
                        if np.any(sel):
                            plt.scatter(XY[sel, 0], XY[sel, 1], s=18, alpha=0.85, label=f"Cluster {cid}", color=colmap[cid])
                    plt.title(f"MDS (distance=1-r), kept clusters | thr={thresh_mult}")
                    plt.xlabel("MDS1")
                    plt.ylabel("MDS2")
                    plt.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_folder, f"ClusterSeparation_MDS_NoSingleton_{method}_th{thresh_mult}.pdf"))
                    plt.savefig(os.path.join(output_folder, f"ClusterSeparation_MDS_NoSingleton_{method}_th{thresh_mult}.png"), dpi=300)
                    plt.show()
                else:
                    print("[SEPARATION] Not enough kept items for MDS (need >=3).")
            
                # ----------------------------
                # Cluster separation visualization 2: Silhouette (distance = 1 - r)
                # ----------------------------
                if idx_kept.size >= 3:
                    labs = cluster_labels[idx_kept]
                    if len(np.unique(labs)) >= 2:
                        sil = silhouette_samples(D, labs, metric="precomputed")
                        sil_score = silhouette_score(D, labs, metric="precomputed")
                        print(f"[SEPARATION] Silhouette score (kept clusters): {sil_score:.3f}")
            
                        sil_df = pd.DataFrame({"cluster": labs, "silhouette": sil})
                        plt.figure(figsize=(7.5, 4.2))
                        sns.boxplot(data=sil_df, x="cluster", y="silhouette", order=plot_clusters)
                        plt.axhline(0, ls="--", c="grey")
                        plt.title(f"Silhouette by cluster (distance=1-r) | score={sil_score:.3f}")
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_folder, f"ClusterSeparation_Silhouette_NoSingleton_{method}_th{thresh_mult}.pdf"))
                        plt.savefig(os.path.join(output_folder, f"ClusterSeparation_Silhouette_NoSingleton_{method}_th{thresh_mult}.png"), dpi=300)
                        plt.show()
            
                        sil_df.to_csv(
                            os.path.join(output_folder, f"ClusterSeparation_SilhouetteValues_NoSingleton_{method}_th{thresh_mult}.csv"),
                            index=False
                        )
                    else:
                        print("[SEPARATION] Silhouette requires >=2 clusters after filtering.")
                else:
                    print("[SEPARATION] Not enough items for silhouette (need >=3).")
            
                # ----------------------------
                # Beeswarm (filtered to KEPT clusters)
                # ----------------------------
                fig, ax = plt.subplots(figsize=(11, 5))
            
                for i, cid in enumerate(plot_clusters):
                    res = results_kept[cid]
                    xj = np.random.normal(i, 0.03, size=res["perm_z"].size)
                    ax.scatter(xj, res["perm_z"], color="grey", alpha=0.25, s=10)
            
                    lo, hi = np.percentile(res["perm_z"], [2.5, 97.5])
                    ax.hlines([lo, hi], i - 0.2, i + 0.2, color="black", lw=2)
            
                    ax.plot(i, res["obs_z"], marker="D", ms=12,
                            color=colmap[cid], markeredgecolor=edge_color(cid))
            
                ymin, ymax = ax.get_ylim()
                y_text = ymin - 1.0
            
                for i, cid in enumerate(plot_clusters):
                    res = results_kept[cid]
                    txtcol = {"attractive": "green", "aversive": "magenta", "neutral": "black"}[res["valence"]]
                    n_trials = int((cluster_labels == cid).sum())
                    n_odors = int(odorant_counts.get(cid, 0))
            
                    ax.text(i, y_text,
                            f"{res['valence']}\np={res['p']:.3g}\nN={n_trials} trials\n({n_odors} odorants)",
                            ha="center", color=txtcol, fontsize=12, fontweight='bold')
            
                ax.axhline(0, ls="--", c="grey")
                ax.set_xticks(range(len(plot_clusters)))
                ax.set_xticklabels(plot_clusters, fontsize=12, fontweight='bold')
                ax.set_ylabel("Preference Z-score", fontsize=14, fontweight='bold')
                ax.set_yticklabels(ax.get_yticks(), fontsize=12)
                ax.set_ylim(y_text - 0.5, ymax + 0.5)
            
                plt.tight_layout()
                plt.savefig(os.path.join(output_folder, f"BeeswarmZ_NoSingleton_{method}_th{thresh_mult}.pdf"))
                plt.show()
            
                # ----------------------------
                # Chemical class analyses (odorant-level) — KEPT clusters only
                # ----------------------------
                odorant_df_kept = odorant_df[odorant_df["cluster"].isin(plot_clusters)].copy()
            
                contig_od = pd.crosstab(odorant_df_kept["cluster"], odorant_df_kept["chemical_class"])
                chi2_od, p_od, dof_od, exp_od = chi2_contingency(contig_od)
                n_od = contig_od.values.sum()
                r, c = contig_od.shape
                cramer_v = np.sqrt(chi2_od / (n_od * (min(r, c) - 1))) if (n_od > 0 and min(r, c) > 1) else np.nan
            
                contig_od.to_csv(os.path.join(output_folder, f"ChemClass_byCluster_ODORANT_LEVEL_NoSingleton_{method}_th{thresh_mult}.csv"))
                print(f"[ODORANT-LEVEL | FILTERED] chi2={chi2_od:.2f}, dof={dof_od}, p={p_od:.3g}, Cramér's V={cramer_v:.3f}")
                print("Min expected:", exp_od.min() if exp_od.size else np.nan)
            
                # Fisher enrichment per (cluster, class), filtered
                from scipy.stats import fisher_exact
                fisher_rows = []
                for cl in contig_od.index:
                    for chem in contig_od.columns:
                        a = contig_od.loc[cl, chem]
                        b = contig_od.loc[cl].sum() - a
                        c_ = contig_od[chem].sum() - a
                        d = contig_od.values.sum() - (a + b + c_)
                        table = [[a, b], [c_, d]]
                        _, p_fisher = fisher_exact(table, alternative="greater")
                        fisher_rows.append({"cluster": cl, "chemical_class": chem, "p_fisher": p_fisher})
            
                fisher_df = pd.DataFrame(fisher_rows)
                rej_f, qvals_f, *_ = multipletests(fisher_df["p_fisher"], method="fdr_bh")
                fisher_df["q_fisher"] = qvals_f
                fisher_df.to_csv(os.path.join(output_folder, f"ChemClass_FisherEnrichment_NoSingleton_{method}_th{thresh_mult}.csv"), index=False)
            
                # Residual heatmap
                assoc_long, assoc_summary = compute_class_pattern_association(
                    odor_df=odorant_df_kept,
                    contig_od=contig_od,
                    outdir=output_folder,
                    alpha=0.05,
                    min_count=1
                )
            
                plt.figure(figsize=(7, 4.5))
                resid_mat = assoc_long.pivot(index="cluster", columns="chemical_class", values="std_resid")
                sns.heatmap(resid_mat, cmap="coolwarm", center=0, annot=True, fmt=".1f",
                            cbar_kws={"label": "Standardized residual (z)"})
                plt.title("Chi-square residuals: clusters × chemical classes (odorant-level, filtered)")
                plt.tight_layout()
                plt.savefig(os.path.join(output_folder, f"ChemClass_ResidualHeatmap_NoSingleton_{method}_th{thresh_mult}.pdf"))
                plt.show()
            
                # Permutation (odorant-level) shuffle classes across odorants
                rng = np.random.default_rng(42)
                n_perm = 10000
                obs_chi2 = chi2_od
                classes = odorant_df_kept["chemical_class"].to_numpy()
                clusters = odorant_df_kept["cluster"].to_numpy()
            
                perm_chi2 = np.empty(n_perm)
                for i in range(n_perm):
                    shuf_classes = rng.permutation(classes)
                    ct = pd.crosstab(clusters, shuf_classes)
                    perm_chi2[i] = chi2_contingency(ct)[0]
            
                p_perm = (np.sum(perm_chi2 >= obs_chi2) + 1) / (n_perm + 1)
                print(f"[ODORANT-LEVEL PERM | FILTERED] p_perm={p_perm:.4g}")
            
                # ----------------------------
                # Odor consistency (trial-level), unchanged but you may prefer filtering to kept clusters:
                # here we compute it on assignments_valid (all trials), and also a filtered view.
                # ----------------------------
                odor_groups = assignments_valid.groupby('odor')
                odor_consistency = odor_groups['cluster'].agg(
                    n_total='count',
                    n_majority=lambda x: x.value_counts().max()
                ).reset_index()
                odor_consistency['consistency_ratio'] = odor_consistency['n_majority'] / odor_consistency['n_total']
            
                # Filtered summary: only odors whose majority cluster is kept
                majority_cluster = (
                    assignments_valid.groupby("odor")["cluster"]
                    .agg(lambda x: x.value_counts().idxmax())
                )
                kept_odors = majority_cluster[majority_cluster.isin(plot_clusters)].index
                odor_consistency_kept = odor_consistency[odor_consistency["odor"].isin(kept_odors)].copy()
            
                cons = odor_consistency_kept["consistency_ratio"].to_numpy()
                if cons.size > 0:
                    cons_med = np.median(cons)
                    cons_iqr = np.percentile(cons, [25, 75])
                    cons_min, cons_max = float(np.min(cons)), float(np.max(cons))
                    print("\n[CONSISTENCY SUMMARY | FILTERED]")
                    print(f"  Majority-cluster ratio: median={cons_med:.2f} (IQR {cons_iqr[0]:.2f}–{cons_iqr[1]:.2f}), range={cons_min:.2f}–{cons_max:.2f}")
            
                    low = odor_consistency_kept.loc[odor_consistency_kept["consistency_ratio"] < 0.60, ["odor","n_total","consistency_ratio"]]
                    if not low.empty:
                        print("\n[LOW CONSISTENCY ODORS | FILTERED] ratio < 0.60")
                        for _, r0 in low.iterrows():
                            print(f"  {r0['odor']}: ratio={r0['consistency_ratio']:.2f} (N={int(r0['n_total'])})")
            
                    with open(os.path.join(output_folder, f"OdorConsistency_NoSingleton_{method}_th{thresh_mult}.txt"), "w", encoding="utf-8") as fh:
                        fh.write(f"Majority-cluster ratio (filtered): median={cons_med:.3f} "
                                 f"(IQR {cons_iqr[0]:.3f}–{cons_iqr[1]:.3f}), range={cons_min:.3f}–{cons_max:.3f}\n")
                        if not low.empty:
                            fh.write("Low-consistency (<0.60):\n")
                            for _, r0 in low.iterrows():
                                fh.write(f"  {r0['odor']}: ratio={r0['consistency_ratio']:.3f} (N={int(r0['n_total'])})\n")
            
                    plt.figure(figsize=(8, 4.5))
                    order = odor_consistency_kept.sort_values("consistency_ratio", ascending=True)
                    sns.barplot(data=order, x="consistency_ratio", y="odor", color="gray")
                    plt.xlabel("Majority-cluster ratio")
                    plt.ylabel("Odor (odorant–concentration)")
                    plt.xlim(0, 1.0)
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_folder, f"OdorConsistency_Bars_NoSingleton_{method}_th{thresh_mult}.pdf"))
                    plt.show()
            
                    plt.figure(figsize=(5, 4))
                    sns.histplot(odor_consistency_kept['consistency_ratio'], bins=10, kde=True)
                    plt.title(f"Consistency (filtered): method={method}, thresh={thresh_mult}")
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_folder, f"OdorConsistency_NoSingleton_{method}_th{thresh_mult}.png"), dpi=300)
                    plt.show()
                else:
                    print("[CONSISTENCY] No odors left in kept clusters for consistency plots.")
            
                # ----------------------------
                # Within-odor vs between-odor similarity (same as yours), unchanged
                # ----------------------------
                labels = pd.Series(row_info["odor"].values, index=matrix_renamed.index, name="odor")
                M = matrix_renamed.values
                same_mask = np.equal.outer(labels.values, labels.values)
                np.fill_diagonal(same_mask, False)
                diff_mask = ~same_mask
            
                within_vals = M[same_mask]
                between_vals = M[diff_mask]
            
                print("\n[WITHIN vs BETWEEN SIMILARITY]")
                print(f"  Within-odor mean={np.mean(within_vals):.3f} ± {np.std(within_vals, ddof=1):.3f} (SD)")
                print(f"  Between-odor mean={np.mean(between_vals):.3f} ± {np.std(between_vals, ddof=1):.3f} (SD)")
            
                rng = np.random.default_rng(0)
                n_perm = 10000
                obs_diff = np.mean(within_vals) - np.mean(between_vals)
                pool = np.concatenate([within_vals, between_vals])
                n_within = within_vals.size
            
                perm_diffs = np.empty(n_perm)
                for i in range(n_perm):
                    shuf = rng.permutation(pool)
                    perm_diffs[i] = shuf[:n_within].mean() - shuf[n_within:].mean()
            
                p_within_gt = (np.sum(perm_diffs >= obs_diff) + 1) / (n_perm + 1)
                print(f"  Permutation p (within > between): p={p_within_gt:.4g}")
            
                plt.figure(figsize=(4.2, 4))
                sns.violinplot(data=[within_vals, between_vals])
                plt.xticks([0, 1], ["Within-odor", "Between-odor"])
                plt.ylabel("Map similarity (Pearson r)")
                plt.tight_layout()
                plt.savefig(os.path.join(output_folder, f"WithinBetweenSimilarity_{method}_th{thresh_mult}.pdf"))
                plt.show()
            
                # ----------------------------
                # Print odorants by kept cluster (odorant-level)
                # ----------------------------
                odor_by_cluster = (
                    odorant_df_kept.groupby("cluster")["odorant"]
                    .apply(list)
                    .reset_index(name="odorants")
                )
                print("\n[ODORANTS BY CLUSTER | FILTERED]")
                print(odor_by_cluster.to_string(index=False))
                print("\n*** Analysis Complete ***")
            
            

#%%

# ----------------------------------------------------------
# EXPORT: Full cluster assignment table (publication-ready)
# ----------------------------------------------------------

cluster_table = pd.DataFrame({
    "matrix_label": matrix_renamed.index,
    "odor": row_info["odor"].values,
    "cluster": cluster_labels,
})

# Trial counts per cluster
trial_counts = pd.Series(cluster_labels).value_counts().to_dict()
cluster_table["cluster_size_n_trials"] = cluster_table["cluster"].map(trial_counts)

# Cluster stability (mean co-cluster probability)
cluster_table["cluster_stability"] = cluster_table["cluster"].map(stability)

# Mark singleton clusters (if any before filtering)
cluster_table["is_singleton_trial"] = cluster_table["cluster"].map(
    lambda c: trial_counts.get(c, 0) < 2
)

# If silhouette was computed:
if 'sil' in locals():
    cluster_table["silhouette_value"] = np.nan
    cluster_table.loc[keep_mask, "silhouette_value"] = sil

# Optional: valence label per cluster
if 'results_kept' in locals():
    valence_map = {cid: results_kept[cid]["valence"] for cid in results_kept}
    cluster_table["cluster_valence_label"] = cluster_table["cluster"].map(valence_map)

out_csv = os.path.join(output_folder, f"ClusterAssignments_FULL_{method}_th{thresh_mult}.csv")
cluster_table.to_csv(out_csv, index=False)

print(f"Full cluster assignment table exported → {out_csv}")
