import os
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import SimpleITK as sitk
from pathlib import Path

# ------------------ CONFIG ------------------
CSV_DIR = r"U:\Scientific Data\RG-AS04-Data01\Oded_Mayseless\imaging_data\proccessed imaging data\campari\202405_campari_analysis_CMTK\images\20241105_analysis\mapZbrain_camp_reg\202501_transformed_csv_files"
OB_MASK_PATH = r"C:\Oded_data\warpfield_reg\images\modified_olfactory_bulb_mask_dilated.tif"
OUTPUT_DIR = r"C:\Oded_data\results_KDE_from_preference_clusters_fixed"

KDE_SIGMA = 5.0

# Filtering params
VOL_TH = 100
VOL_MAX = 1000  # 
DISTANCE_THRESHOLD = 10
INTENSITY_COL = 'mean_intensity_channel_2'
QUANTILE = 0.9
DISTANCE_COLUMNS = ['mdG2', 'maG', 'mdG6', 'dG', 'vmG', 'lG', 'dlG', 'vpG']

STORE_FILENAME = "./campari-data-analysis-master/data/data.h5"
FSN = "dmin_0__dmax_10__cellvol_100_1000__topp_10"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------
# Build (odor, ds_idx) -> csv_filename map
# ------------------
def build_idx_to_csv_map(store_path: str, fsn: str) -> pd.DataFrame:
    store = pd.HDFStore(store_path, mode="r")
    try:
        odors = [node._v_name for node in store.get_node(f"/{fsn}") if node._v_name.startswith("st_")]
        rows = []
        for st_name in odors:
            odor = st_name[3:]
            fn_df = store.get(f"/{fsn}/{st_name}/filenames").rename(columns={"idx": "ds_idx"}).copy()
            fn_df["odor"] = odor
            rows.append(fn_df)
        return pd.concat(rows, ignore_index=True)  # ds_idx, csv_filename, odor
    finally:
        store.close()

idx_to_csv = build_idx_to_csv_map(STORE_FILENAME, FSN)

def parse_index_row(row_name: str):
    parts = row_name.split(',')
    idx = int(parts[0])
    odor = parts[1].strip()
    ds = parts[2].strip() if len(parts) > 2 else "ds_?"
    return idx, odor, ds

# ---- assumes matrix and cluster_labels are already in memory ----
rows_parsed = [parse_index_row(r) for r in matrix.index]
df_members = pd.DataFrame(rows_parsed, columns=["idx", "odor", "ds"])
df_members["cluster"] = cluster_labels
df_members["ds_idx"] = df_members["ds"].str.extract(r"(\d+)").astype(int)

df_members_with_csv = df_members.merge(idx_to_csv, on=["odor", "ds_idx"], how="left")

# check missing BEFORE casting to str
missing = df_members_with_csv[df_members_with_csv["csv_filename"].isna()]
if len(missing):
    print("⚠️ Rows that could not be matched:\n", missing[["odor", "ds"]].drop_duplicates())
else:
    print("✅ All df_members rows matched to original CSVs.")

# full paths (use CSV_DIR!)
data_path = Path(CSV_DIR)
df_members_with_csv["csv_path"] = df_members_with_csv["csv_filename"].map(
    lambda f: str(data_path / f) if pd.notna(f) else np.nan
)

# drop rows with missing or nonexistent files
df_members_with_csv = df_members_with_csv.dropna(subset=["csv_path"]).copy()
exists_mask = df_members_with_csv["csv_path"].map(os.path.exists)
if (~exists_mask).any():
    print("⚠️ Some csv_path do not exist; dropping them:")
    print(df_members_with_csv.loc[~exists_mask, ["odor", "ds", "csv_path"]].head(10).to_string(index=False))
df_members_with_csv = df_members_with_csv.loc[exists_mask].copy()

# ------------------
# Helpers
# ------------------
def load_sitk(path, dtype=np.float32):
    return sitk.GetArrayFromImage(sitk.ReadImage(path)).astype(dtype)

def save_volume_as_nrrd_like(vol, path, reference_img=None, dtype=np.float32):
    img = sitk.GetImageFromArray(vol.astype(dtype))
    if reference_img is not None:
        img.CopyInformation(reference_img)
    sitk.WriteImage(img, path)

def compute_kde(points_zyx, vol_shape, sigma):
    vol = np.zeros(vol_shape, dtype=np.float32)
    if len(points_zyx) > 0:
        z = np.clip(points_zyx[:, 0], 0, vol_shape[0]-1)
        y = np.clip(points_zyx[:, 1], 0, vol_shape[1]-1)
        x = np.clip(points_zyx[:, 2], 0, vol_shape[2]-1)
        vol[z, y, x] = 1.0
    kde = gaussian_filter(vol, sigma=sigma)
    m = kde.max()
    if m > 0:
        kde /= m
    return kde

def filter_and_extract_points(files, vol_shape):
    all_pts = []
    for f in files:
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"Skipping {f}: {e}")
            continue

        need = {'x','y','z','volume', INTENSITY_COL}
        if not need.issubset(df.columns):
            continue

        # NOTE: this is AND logic (volume + distance). Adjust if you want OR.
        df = df[(df['volume'] >= VOL_TH) & (df['volume'] < VOL_MAX)]
        df = df[df[DISTANCE_COLUMNS].min(axis=1) <= DISTANCE_THRESHOLD]

        if len(df) == 0:
            continue

        min_int = df[INTENSITY_COL].quantile(QUANTILE)
        max_int = df[INTENSITY_COL].quantile(0.99)
        df = df[(df[INTENSITY_COL] >= min_int) & (df[INTENSITY_COL] <= max_int)]

        df['x'] = df['x'].round().astype(int)
        df['y'] = df['y'].round().astype(int)
        df['z'] = df['z'].round().astype(int)

        valid = (
            (df['z'] >= 0) & (df['z'] < vol_shape[0]) &
            (df['y'] >= 0) & (df['y'] < vol_shape[1]) &
            (df['x'] >= 0) & (df['x'] < vol_shape[2])
        )

        coords = df.loc[valid, ['z', 'y', 'x']].to_numpy()
        if len(coords) > 0:
            all_pts.append(coords)

    return np.vstack(all_pts) if all_pts else np.empty((0,3), dtype=int)

# ------------------
# Generate one KDE NRRD per cluster
# ------------------
ob_mask = load_sitk(OB_MASK_PATH, dtype=bool)
vol_shape = ob_mask.shape
ref_img = sitk.ReadImage(OB_MASK_PATH)

unique_clusters = sorted(df_members_with_csv["cluster"].dropna().unique())
print(f"[INFO] Generating KDE NRRDs for {len(unique_clusters)} clusters: {unique_clusters}")

for cid in unique_clusters:
    cluster_df = df_members_with_csv[df_members_with_csv["cluster"] == cid]
    files = cluster_df["csv_path"].tolist()

    print(f"\n[INFO] Cluster {cid}: {len(files)} CSVs")
    pts = filter_and_extract_points(files, vol_shape)
    print(f"  points after filtering: {pts.shape[0]}")

    if pts.shape[0] == 0:
        print("  ⚠️ Skipping (no points)")
        continue

    kde = compute_kde(pts, vol_shape, KDE_SIGMA)
    kde_out = kde * ob_mask  # OB-restricted

    out_path = os.path.join(OUTPUT_DIR, f"cluster_{int(cid)}_KDE_sigma{KDE_SIGMA}_d{DISTANCE_THRESHOLD}.nrrd")
    save_volume_as_nrrd_like(kde_out, out_path, reference_img=ref_img)
    print(f"  ✅ Saved: {out_path}")

print(f"\n[DONE] Saved cluster KDE NRRDs to: {OUTPUT_DIR}")