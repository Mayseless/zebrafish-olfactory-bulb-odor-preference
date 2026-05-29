# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import defaultdict
from scipy.stats import gaussian_kde
from matplotlib.animation import FuncAnimation, FFMpegWriter


def collect_points_by_odor(
    csv_folder,
    intensity_col="mean_intensity_channel_2",
    quantile=0.1,
    upper_quantile=0.99,
    distance_threshold=30,
    vol_th=100,
    min_intensity=0,
    flip_z=False,
    z_bump = 20,
    z_limits=(50, 250),
):
    """
    Loads all CSVs, groups them by odor, applies filtering, and returns
    one DataFrame per odor.
    """

    distance_columns = ["mdG2", "maG", "mdG6", "dG", "vmG", "lG", "dlG", "vpG"]

    files_by_odor = defaultdict(list)

    for csv_file in os.listdir(csv_folder):
        if csv_file.endswith("_cell_distances_intensities.csv"):
            odor_name = csv_file.split("60")[0]
            files_by_odor[odor_name].append(os.path.join(csv_folder, csv_file))

    points_by_odor = {}

    for odor_name, csv_files in files_by_odor.items():
        odor_points = []

        for csv_path in csv_files:
            cell_data_df = pd.read_csv(csv_path)

            lower_q = cell_data_df[intensity_col].quantile(quantile)
            lower_threshold = max(lower_q, min_intensity)
            upper_threshold = cell_data_df[intensity_col].quantile(upper_quantile)

            filtered = cell_data_df[
                (cell_data_df[intensity_col] >= lower_threshold)
                & (cell_data_df[intensity_col] < upper_threshold)
                & (cell_data_df[distance_columns].min(axis=1) <= distance_threshold)
                & (cell_data_df["volume"] >= vol_th)
            ].copy()

            # Flip Z coordinate if needed
            z_bump = z_bump

            if flip_z:
                filtered["centroid_z"] = z_limits[1] - (
                    filtered["centroid_z"] - z_limits[0]
                ) + z_bump

            filtered["source_file"] = os.path.basename(csv_path)
            odor_points.append(filtered)

        if len(odor_points) > 0:
            odor_df = pd.concat(odor_points, ignore_index=True)
        else:
            odor_df = pd.DataFrame(
                columns=["centroid_x", "centroid_y", "centroid_z", intensity_col]
            )

        points_by_odor[odor_name] = {
            "df": odor_df,
            "n_files": len(csv_files),
            "n_spots": len(odor_df),
        }

        print(f"{odor_name}: {len(odor_df)} spots from {len(csv_files)} files")

    return points_by_odor

def downsample_points_by_odor(points_by_odor, max_points_per_odor=1500, seed=1):
    """
    Makes a lightweight copy of points_by_odor with fewer points per odor.
    Useful for previewing movie settings quickly.
    """
    rng = np.random.default_rng(seed)
    small = {}

    for odor, info in points_by_odor.items():
        df = info["df"]

        if len(df) > max_points_per_odor:
            idx = rng.choice(df.index.to_numpy(), size=max_points_per_odor, replace=False)
            df_small = df.loc[idx].copy()
        else:
            df_small = df.copy()

        small[odor] = {
            "df": df_small,
            "n_files": info["n_files"],
            "n_spots": len(df_small),
            "n_spots_original": info["n_spots"],
        }

    return small

def make_odor_rotation_movie(
    points_by_odor,
    output_path,
    intensity_col="mean_intensity_channel_2",
    x_limits=(150, 350),
    y_limits=(150, 250),
    z_limits=(50, 250),
    elev=25,
    azim_start=-70,
    azim_end=290,
    frames_per_odor=48,
    fps=12,
    scatter_size=10,
    scatter_alpha=0.8,
    vmax=2500,
    make_density_projection=True,
    density_grid=30,
):
    """
    Makes a movie that flips through odors while rotating the view horizontally.

    elev:
        Fixed vertical viewing angle.

    azim_start / azim_end:
        Horizontal rotation range.

    frames_per_odor:
        Number of movie frames shown per odor.
        Higher = slower odor transitions and smoother rotation.

    fps:
        Frames per second in final movie.
    """

    odor_names = sorted(points_by_odor.keys())
    n_odors = len(odor_names)

    if n_odors == 0:
        raise ValueError("No odors found.")

    total_frames = n_odors * frames_per_odor

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Keep color scale fixed across odors
    vmax = 2500
    norm = None
    cmap = "plasma"    # Dummy scatter for a stable colorbar

    dummy = ax.scatter([], [], [], c=[], cmap=cmap, vmin = 0, vmax = 2500)
    colorbar = fig.colorbar(
        dummy,
        ax=ax,
        orientation="vertical",
        label="Activity A.U.",
        shrink=0.5,
    )
    colorbar.ax.set_position([0.85, 0.25, 0.05, 0.5])

    def draw_frame(frame_idx):
        ax.clear()

        odor_idx = frame_idx // frames_per_odor
        within_odor_frame = frame_idx % frames_per_odor

        odor_name = odor_names[odor_idx]
        odor_info = points_by_odor[odor_name]
        df = odor_info["df"]

        # Rotate smoothly within each odor
        frac = within_odor_frame / max(frames_per_odor - 1, 1)
        azim = azim_start + frac * (azim_end - azim_start)

        if len(df) > 0:
            x = df["centroid_x"].to_numpy()
            y = df["centroid_y"].to_numpy()
            z = df["centroid_z"].to_numpy()
            intensity = df[intensity_col].to_numpy()

            ax.scatter(
                x,
                y,
                z,
                alpha=scatter_alpha,
                c=intensity,
                cmap=cmap,
                vmax = 2500,
                s=scatter_size,
            )

            if make_density_projection and len(df) >= 5:
                try:
                    gx, gy, gz = np.mgrid[
                        x_limits[0]:x_limits[1]:complex(density_grid),
                        y_limits[0]:y_limits[1]:complex(density_grid),
                        z_limits[0]:z_limits[1]:complex(density_grid),
                    ]

                    positions = np.vstack([gx.ravel(), gy.ravel(), gz.ravel()])
                    xyz = np.vstack([x, y, z])

                    density = gaussian_kde(xyz)
                    values = density(positions)
                    density_values = values.reshape(gx.shape)

                    # Z-axis projection
                    ax.contour(
                        gx[:, :, 0],
                        gy[:, :, 0],
                        np.sum(density_values, axis=2),
                        zdir="z",
                        offset=z_limits[0],
                        cmap=cmap,
                        alpha=0.9,
                    )

                    # Y-axis projection
                    ax.contour(
                        gx[:, :, 0],
                        np.sum(density_values, axis=1),
                        gz[0, :, :],
                        zdir="y",
                        offset=y_limits[1],
                        cmap=cmap,
                        alpha=0.9,
                    )

                    # X-axis projection
                    ax.contour(
                        np.sum(density_values, axis=0),
                        gy[0, :, :],
                        gz[0, :, :],
                        zdir="x",
                        offset=x_limits[0],
                        cmap=cmap,
                        alpha=0.9,
                    )

                except Exception as e:
                    print(f"Skipping density for {odor_name}: {e}")

        ax.set_xlim(x_limits)
        ax.set_ylim(y_limits)
        ax.set_zlim(z_limits)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        ax.view_init(elev=elev, azim=azim)

        ax.set_title(
            f"{odor_name}\n"
            f"n = {odor_info['n_files']}",
            fontsize=24,
        )

        # Optional: reduce visual clutter
        ax.grid(False)

        return []

    animation = FuncAnimation(
        fig,
        draw_frame,
        frames=total_frames,
        interval=1000 / fps,
        blit=False,
    )

    writer = FFMpegWriter(fps=fps, bitrate=3000)
    animation.save(output_path, writer=writer)

    plt.close(fig)

    print(f"Saved movie to:\n{output_path}")


# %%
# Example usage

#csv_folder = r"C:\Oded_data\campari_pipeline\gitlab\manualDownload\202501_transformed_csv_files"
output_folder = os.path.join(csv_folder, "output_movies")
os.makedirs(output_folder, exist_ok=True)

points_by_odor = collect_points_by_odor(
    csv_folder=csv_folder,
    intensity_col="mean_intensity_channel_2",
    quantile=0.1,
    upper_quantile=0.99,
    distance_threshold=30,
    vol_th=100,
    min_intensity=0,
    flip_z=True,
    z_bump= 20,
    z_limits=(50, 250),
)

movie_path = os.path.join(output_folder, "odor_scatter_rotation_movie.mp4")

make_odor_rotation_movie(
    points_by_odor=points_by_odor,
    output_path=movie_path,
    intensity_col="mean_intensity_channel_2",
    x_limits=(150, 350),
    y_limits=(150, 250),
    z_limits=(50, 250),

    # point of view
    elev=25,
    azim_start=-70,
    azim_end=290,

    # timing
    frames_per_odor=36,
    fps=12,

    # plot appearance
    scatter_size=10,
    scatter_alpha=0.8,
    vmax=2500,

    # set False if the movie is too slow to render
    make_density_projection=True,
    density_grid=30,
)