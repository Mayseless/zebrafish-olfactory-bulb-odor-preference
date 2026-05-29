#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 13:24:27 2024

@author: mazur
"""

# %%
import os
import glob
import itertools
import json
import pandas as pd
from collections import OrderedDict

# %%


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


# %%
STORE_FILENAME = "./campari-data-analysis-master/data/data_AND.h5"

if os.path.exists(STORE_FILENAME):
    print(f"Removing file: {STORE_FILENAME}")
    os.remove(STORE_FILENAME)

# %%
params = {
    'DATA_PATH': "./202501_transformed_csv_files",
    'D_THRESHOLD': [10], #[10, 15],
    'LOWER_D_THRESHOLD': [0],
    'CELL_VOLUME_MIN': 100, # expected cell volume ~500 um^3: 4/3 * \pi * 5^3 = 523.6
    'CELL_VOLUME_MAX': 1000,
    'TOP_PERCENTAGE': [10] # , 10, 25]
}

# %%


def extract_odor_names_from_files(csv_directory):
    odor_names_from_files = set()
    for filename in os.listdir(csv_directory):
        if filename.endswith("_cell_distances_intensities.csv"):
            # Assume the odor name is before the first underscore
            odor_name = filename.split('60sec')[0]
            odor_names_from_files.add(odor_name)
    return odor_names_from_files
# %%


def load_clean(fname, lower_d_threshold, d_threshold, cell_volume_min, cell_volume_max, top_percentage):
    # glomeruli == common landmarks
    GLOMERULI = ['mdG2', 'maG', 'mdG6', 'dG', 'vmG', 'lG', 'dlG', 'vpG']

    if lower_d_threshold >= d_threshold:
        print("Lower threshold is greater than or equal to upper threshold. Exiting!")
        exit(-1)

    raw_df = pd.read_csv(fname)

    # Initialize mask as False for all rows
    # mask = pd.Series([False] * raw_df.shape[0])

    # initialize mask by cell volume
    mask = (raw_df['volume'] >= cell_volume_min) & (raw_df['volume'] < cell_volume_max)

    # Use bitwise OR to combine masks, so the mask is True if the value is within the thresholds
    for g in GLOMERULI:
        mask = mask | (raw_df[g] >= lower_d_threshold) & (raw_df[g] < d_threshold)

    df = raw_df[mask]

    # Calculate the intensity threshold for the top percentage
    ich2_threshold = df['mean_intensity_channel_2'].quantile(
        1 - top_percentage/100)
    df = df[df['mean_intensity_channel_2'] >= ich2_threshold]

    df = df[['cell_label_id', 'centroid_x', 'centroid_y',
             'centroid_z', 'volume', 'mean_intensity_channel_2']]

    return df



# %% STIMULANTS definition here
stimulant_names = list(extract_odor_names_from_files(params["DATA_PATH"]))

# %%
stimulants_combinations = list(zip(
    stimulant_names, stimulant_names)) + list(itertools.combinations(stimulant_names, 2))

# %%
stimulants_csvs = OrderedDict((stn, [os.path.basename(fn)
                                     for fn in sorted(glob.glob(f'{params["DATA_PATH"]}/{stn}*.csv'))])
                              for stn in stimulant_names)

# %%
store = pd.HDFStore(STORE_FILENAME)

# %%
params["STIMULANTS"] = stimulant_names
store.put("parameters", pd.DataFrame([("json", json.dumps(params))]))

# %%
for top_percentage in params["TOP_PERCENTAGE"]:
    cell_volume_min = params["CELL_VOLUME_MIN"]
    cell_volume_max = params["CELL_VOLUME_MAX"]
    
    for dmin in params["LOWER_D_THRESHOLD"]:
        for dmax in params["D_THRESHOLD"]:
            if dmin >= dmax:
                continue
            
            params_str = f"dmin_{dmin}__dmax_{dmax}__cellvol_{cell_volume_min}_{cell_volume_max}__topp_{top_percentage}"
            print(params_str)
            for st_name in stimulant_names:
                idx_csvname = []

                for idx, csv in enumerate(stimulants_csvs[st_name]):
                    idx_csvname.append([idx, csv])
                    fname = params['DATA_PATH'] + "/" + csv
                    df = load_clean(fname, dmin, dmax, cell_volume_min, cell_volume_max, top_percentage)
                    store.put(f"{params_str}/st_{st_name}/data/ds_{idx}", df)

                store.put(f"{params_str}/st_{st_name}/filenames",
                          pd.DataFrame(idx_csvname, columns=("idx", "csv_filename")))

# %%
store.close()

# %%
print("Done.")
