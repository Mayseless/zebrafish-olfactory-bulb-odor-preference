

# %%
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity

import os
from collections import OrderedDict
# from tqdm import tqdm
import sys
sys.path.append('C:/Oded_data/campari_pipeline/gitlab/manualDownload/campari-data-analysis-master/')

import campari

# %%
STORE_FILENAME = "./campari-data-analysis-master/data/data.h5"
FSN = "dmin_0__dmax_10__cellvol_100_1000__topp_10"  # Filter Settings Name

BW = 5.0  # bandwidth parameter

# OUT_FILENAME = f"out/output_{FSN}.h5"
OUT_FILENAME = f"./campari-data-analysis-master/out/2_output_dmax10.h5"

GRID_OFFSET, GRID_SPACING = 10.0, 5.0

# %%
if os.path.exists(OUT_FILENAME):
    ans = input(f"{OUT_FILENAME} exists. Delete it? [y/N] ").lower()
    if ans == "y":
        os.remove(OUT_FILENAME)

# %%
store = pd.HDFStore(STORE_FILENAME, 'r')

# %% get stimulant names
stimulants = {}

flat_index = []  # (idx, stn, ds_n)
flat_index_str = []

idx = 0
for st in store.get_node(f"/{FSN}"):
    stn = st._v_name  # stn := stimulant name

    print(stn)
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

# %% prepare summed up stimulants data
xyz_w_kdes = {}

for stn in stimulants.keys():
    dsns = stimulants[stn].keys()

    xyz_w_kdes_stn = {}

    for dsn in dsns:
        ds = stimulants[stn][dsn]
        
        xyz = np.array((ds.centroid_x, ds.centroid_y, ds.centroid_z)).T
        w = (ds.mean_intensity_channel_2 *
             ds.shape[0])/ds.mean_intensity_channel_2.sum()

        xyz_w_kdes_stn[dsn] = (xyz, w, KernelDensity(kernel='gaussian', bandwidth=BW).fit(
            xyz, sample_weight=w))

    xyz_w_kdes[stn] = xyz_w_kdes_stn

# %% Open an HDF5 store for output
store_output = pd.HDFStore(OUT_FILENAME)

# %% Calculate statistics / stimulant
for stn in stimulants.keys():
    # datasets names (dsns)
    dsns = tuple(stimulants[stn].keys())
    N_dsns = len(dsns)

    # alloc matrices
    xc = np.ones((N_dsns, N_dsns))
    xcm = np.ones((N_dsns, N_dsns))
    So = np.ones((N_dsns, N_dsns))
    rch = np.zeros((N_dsns, N_dsns))  # Chamfer distance

    # calculate matrices
    for i_dsn in range(len(dsns)):
        dsn_i = dsns[i_dsn]

        for j_dsn in range(i_dsn+1, len(dsns)):
            print(f"{stn}: ({i_dsn}, {j_dsn})")
            dsn_j = dsns[j_dsn]
            #
            R, R_masked, S_overlap = campari.calc_xcorr_overlap(
                xyz_w_kdes[stn][dsn_i], xyz_w_kdes[stn][dsn_j])

            r_chamfer = campari.chamfer_like_distance(
                xyz_w_kdes[stn][dsn_i][0], xyz_w_kdes[stn][dsn_j][0])

            xc[i_dsn, j_dsn] = R
            xcm[i_dsn, j_dsn] = R_masked
            So[i_dsn, j_dsn] = S_overlap
            rch[i_dsn, j_dsn] = r_chamfer

            xc[j_dsn, i_dsn] = xc[i_dsn, j_dsn]
            xcm[j_dsn, i_dsn] = xcm[i_dsn, j_dsn]
            So[j_dsn, i_dsn] = So[i_dsn, j_dsn]
            rch[j_dsn, i_dsn] = rch[i_dsn, j_dsn]

    df_xc = pd.DataFrame(xc, columns=dsns, index=dsns)
    df_xcm = pd.DataFrame(xcm, columns=dsns, index=dsns)
    df_So = pd.DataFrame(So, columns=dsns, index=dsns)
    df_rch = pd.DataFrame(rch, columns=dsns, index=dsns)

    # save to the store
    store_output.put(f"/{FSN}/{stn}/xc", df_xc)
    store_output.put(f"/{FSN}/{stn}/xcm", df_xcm)
    store_output.put(f"/{FSN}/{stn}/Soverlap", df_So)
    store_output.put(f"/{FSN}/{stn}/r_chamfer", df_rch)


# %%
store_output.close()
print("done")