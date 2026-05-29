
# %%
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from collections import OrderedDict
from tqdm import tqdm
import sys
sys.path.append('C:/Oded_data/campari_pipeline/gitlab/manualDownload/campari-data-analysis-master/')

import campari

# %%

STORE_FILENAME = "./campari-data-analysis-master/data/data.h5"
FSN = "dmin_0__dmax_10__cellvol_100_1000__topp_10"
BW = 5.0
OUT_FILENAME = "./campari-data-analysis-master/out/3_output_dmax10.h5"

if os.path.exists(OUT_FILENAME):
    ans = input(f"{OUT_FILENAME} exists. Delete it? [y/N] ").lower()
    if ans == "y":
        os.remove(OUT_FILENAME)


# %%
store = pd.HDFStore(STORE_FILENAME, 'r')

# %%
# print(store.keys())

# %%
store_output = pd.HDFStore(OUT_FILENAME)


# %%

# %% get stimulant names
stimulants = {}

flat_index = []  # (idx, stn, ds_n)
flat_index_str = []

idx = 0
for st in store.get_node(f"/{FSN}"):
    stn = st._v_name

    tmp = []
    for ds in store.get_node(f"/{FSN}/{stn}/data"):
        ds_n = ds._v_name
        tmp.append((ds_n, store.get(f"/{FSN}/{stn}/data/{ds_n}")))
        flat_index.append((idx, stn, ds_n))
        flat_index_str.append(f"{idx},{stn[3:]},{ds_n}")
        idx += 1
    stimulants[stn] = OrderedDict(tmp)

# %%
store.close()

# %% prepare a NumPY array
N = sum([len(x) for x in stimulants.values()])

# Chamfer distance
Dch_matrix = np.zeros((N, N))
# Overlap integral
S_matrix = np.ones((N, N))
# Cross correlation matrices
XC_matrix = np.zeros((N, N))
XCM_matrix = np.zeros((N, N))

# %%
for i, sti, dsi in tqdm(flat_index):
    print(f"{i+1}/{len(flat_index)}")
    ds_i = stimulants[sti][dsi]

    pc_i = np.array((ds_i.centroid_x, ds_i.centroid_y, ds_i.centroid_z)).T
    w_i = (ds_i.mean_intensity_channel_2 *
           ds_i.shape[0])/ds_i.mean_intensity_channel_2.sum()
    kde_i = KernelDensity(kernel='gaussian', bandwidth=BW).fit(
        pc_i, sample_weight=w_i)

    combined_i = pc_i, w_i, kde_i

    for j, stj, dsj in tqdm(flat_index[i:]):
        ds_j = stimulants[stj][dsj]

        pc_j = np.array((ds_j.centroid_x, ds_j.centroid_y, ds_j.centroid_z)).T
        w_j = (ds_j.mean_intensity_channel_2 *
               ds_j.shape[0])/ds_j.mean_intensity_channel_2.sum()
        kde_j = KernelDensity(kernel='gaussian', bandwidth=BW).fit(
            pc_j, sample_weight=w_j)

        combined_j = pc_j, w_j, kde_j

        # Chamfer distance
        Dch_matrix[i, j] = campari.chamfer_like_distance(
            pc_i, pc_j, averaging='mean')
        # Dch_matrix[i, j] = campari.chamfer_like_distance(pc_i, pc_j, averaging='median')
        Dch_matrix[j, i] = Dch_matrix[i, j]

        # Overlap integral
        R, R_masked, S_overlap = campari.calc_xcorr_overlap(
            combined_i, combined_j)

        XC_matrix[i, j] = R
        XC_matrix[j, i] = R
        XCM_matrix[i, j] = R_masked
        XCM_matrix[j, i] = R_masked
        S_matrix[i, j] = S_overlap
        S_matrix[j, i] = S_overlap

# %% save the distance matrix
df_Dch = pd.DataFrame(Dch_matrix, columns=flat_index_str, index=flat_index_str)
df_XC = pd.DataFrame(XC_matrix, columns=flat_index_str, index=flat_index_str)
df_XCM = pd.DataFrame(XCM_matrix, columns=flat_index_str, index=flat_index_str)
df_S = pd.DataFrame(S_matrix, columns=flat_index_str, index=flat_index_str)

# %%
store_output.put(f"/{FSN}/Dch_matrix", df_Dch)
store_output.put(f"/{FSN}/XC_matrix", df_XC)
store_output.put(f"/{FSN}/XCM_matrix", df_XCM)
store_output.put(f"/{FSN}/S_matrix", df_S)

# %%
store_output.close()

# %%
print("Done.")
