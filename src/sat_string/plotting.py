"""Plots and H.264 animation regenerated exclusively from saved histories."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from typing import Mapping

import imageio_ffmpeg

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sat_string_matplotlib")
)
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
from numpy.typing import NDArray

from .atmosphere import background_density
from .config import SimulationConfig
from .disturbance import relative_density_increment, rotating_coordinate_m


REQUIRED_PLOT_FILENAMES = (
    "gap_error_spacetime.png",
    "max_gap_error.png",
    "mode_spacetime.png",
    "gap_and_mode_spacetime.png",
    "relative_altitude.png",
    "altitude_cost_vs_baseline.png",
    "minimum_gap.png",
    "density_input_spacetime.png",
)


class AnimationGenerationError(RuntimeError):
    """Raised only when MP4 encoding fails after numerical artifacts exist."""


def _save_figure(figure: plt.Figure, path: Path, dpi: int) -> None:
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _diverging_norm(values: NDArray[np.float64]) -> TwoSlopeNorm:
    limit = float(np.max(np.abs(values)))
    if limit == 0.0:
        limit = 1.0
    return TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)


def _spacetime(
    axis: plt.Axes,
    values: NDArray,
    time_h: NDArray[np.float64],
    *,
    cmap,
    norm=None,
    ylabel: str,
):
    rows = values.shape[1]
    image = axis.imshow(
        values.T,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=(time_h[0], time_h[-1], 0.5, rows + 0.5),
        cmap=cmap,
        norm=norm,
    )
    axis.set_ylabel(ylabel)
    axis.set_yticks(np.arange(1, rows + 1))
    return image


def generate_standard_plots(
    controlled: Mapping[str, NDArray],
    undisturbed: Mapping[str, NDArray],
    low_drag_baseline: Mapping[str, NDArray],
    config: SimulationConfig,
    output_directory: str | Path,
    *,
    dpi: int = 300,
) -> list[Path]:
    """Generate all eight required static v1 figures."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    time = np.asarray(controlled["time_s"], dtype=float)
    time_h = time / 3600.0
    error = np.asarray(controlled["gap_error_m"], dtype=float)
    modes = np.asarray(controlled["mode"], dtype=np.int8)
    paths: list[Path] = []

    figure, axis = plt.subplots(figsize=(10, 5))
    image = _spacetime(
        axis,
        error,
        time_h,
        cmap="RdBu_r",
        norm=_diverging_norm(error),
        ylabel="Gap number",
    )
    axis.set_xlabel("Time (h)")
    figure.colorbar(image, ax=axis, label="Gap error (m)")
    path = output / "gap_error_spacetime.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    raw_max = np.max(np.abs(error), axis=0)
    effect = error - np.asarray(undisturbed["gap_error_m"], dtype=float)
    effect_max = np.max(np.abs(effect), axis=0)
    gaps = np.arange(1, error.shape[1] + 1)
    figure, axis = plt.subplots(figsize=(9, 5))
    axis.plot(gaps, raw_max, "o-", label="Controlled, disturbed")
    axis.plot(gaps, effect_max, "s--", label="Disturbed - undisturbed")
    axis.set(xlabel="Gap number", ylabel="Maximum absolute gap error (m)")
    axis.set_xticks(gaps)
    axis.grid(True, alpha=0.3)
    axis.legend()
    path = output / "max_gap_error.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    mode_cmap = ListedColormap(["white", "#f4a261", "#1f2937", "#2a9d8f"])
    mode_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], mode_cmap.N)
    figure, axis = plt.subplots(figsize=(10, 5))
    image = _spacetime(
        axis,
        modes,
        time_h,
        cmap=mode_cmap,
        norm=mode_norm,
        ylabel="Satellite number",
    )
    axis.set_xlabel("Time (h)")
    colorbar = figure.colorbar(image, ax=axis, ticks=[0, 1, 2, 3])
    colorbar.ax.set_yticklabels(["L", "L_TO_H", "H", "H_TO_L"])
    path = output / "mode_spacetime.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    top = _spacetime(
        axes[0],
        error,
        time_h,
        cmap="RdBu_r",
        norm=_diverging_norm(error),
        ylabel="Gap number",
    )
    bottom = _spacetime(
        axes[1],
        modes,
        time_h,
        cmap=mode_cmap,
        norm=mode_norm,
        ylabel="Satellite number",
    )
    axes[1].set_xlabel("Time (h)")
    figure.colorbar(top, ax=axes[0], label="Gap error (m)")
    mode_bar = figure.colorbar(bottom, ax=axes[1], ticks=[0, 1, 2, 3])
    mode_bar.ax.set_yticklabels(["L", "L_TO_H", "H", "H_TO_L"])
    path = output / "gap_and_mode_spacetime.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    radii = np.asarray(controlled["a_m"], dtype=float)
    relative = radii - np.mean(radii, axis=1, keepdims=True)
    spread = np.ptp(relative, axis=1)
    figure, axis = plt.subplots(figsize=(10, 5))
    for index in range(relative.shape[1]):
        axis.plot(time_h, relative[:, index], linewidth=0.8, alpha=0.75)
    axis.plot(time_h, spread, color="black", linewidth=2.0, label="Altitude spread")
    axis.set(xlabel="Time (h)", ylabel="Relative altitude / spread (m)")
    axis.grid(True, alpha=0.3)
    axis.legend()
    path = output / "relative_altitude.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    altitude_loss = -(
        radii - np.asarray(low_drag_baseline["a_m"], dtype=float)
    )
    figure, axis = plt.subplots(figsize=(10, 5))
    for index in range(altitude_loss.shape[1]):
        axis.plot(time_h, altitude_loss[:, index], linewidth=0.9, label=f"Sat {index + 1}")
    axis.set(xlabel="Time (h)", ylabel="Altitude loss vs all-low baseline (m)")
    axis.grid(True, alpha=0.3)
    if altitude_loss.shape[1] <= 12:
        axis.legend(ncol=2, fontsize="small")
    path = output / "altitude_cost_vs_baseline.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    minimum_gap = np.min(np.asarray(controlled["gap_distance_m"], dtype=float), axis=1)
    figure, axis = plt.subplots(figsize=(10, 5))
    axis.plot(time_h, minimum_gap / 1000.0)
    axis.set(xlabel="Time (h)", ylabel="Minimum adjacent spacing (km)")
    axis.grid(True, alpha=0.3)
    path = output / "minimum_gap.png"
    _save_figure(figure, path, dpi)
    paths.append(path)

    background = background_density(
        np.asarray(controlled["altitude_m"], dtype=float),
        config.formation.reference_altitude_m,
        config.atmosphere.reference_density_kg_m3,
        config.atmosphere.effective_scale_height_m,
    )
    relative_density = np.asarray(controlled["density_kg_m3"], dtype=float) / background - 1.0
    figure, axis = plt.subplots(figsize=(10, 5))
    image = _spacetime(
        axis,
        relative_density,
        time_h,
        cmap="magma",
        ylabel="Satellite number",
    )
    axis.set_xlabel("Time (h)")
    figure.colorbar(image, ax=axis, label="Relative density increment (-)")
    path = output / "density_input_spacetime.png"
    _save_figure(figure, path, dpi)
    paths.append(path)
    return paths


def nearest_unique_frame_indices(
    saved_time_s: NDArray[np.float64], frame_step_s: float
) -> NDArray[np.int64]:
    """Map requested frame times to nearest saved samples without duplicates."""
    time = np.asarray(saved_time_s, dtype=float)
    targets = np.arange(time[0], time[-1] + 0.5 * frame_step_s, frame_step_s)
    if targets.size == 0 or targets[-1] < time[-1] - 1e-9:
        targets = np.append(targets, time[-1])
    right = np.searchsorted(time, targets, side="left")
    right = np.clip(right, 0, len(time) - 1)
    left = np.clip(right - 1, 0, len(time) - 1)
    choose_right = np.abs(time[right] - targets) < np.abs(time[left] - targets)
    indices = np.where(choose_right, right, left)
    return np.unique(indices).astype(np.int64)


def generate_formation_animation(
    arrays: Mapping[str, NDArray],
    config: SimulationConfig,
    output_path: str | Path,
) -> Path:
    """Encode the required formation and moving-density MP4 from saved data."""
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    time = np.asarray(arrays["time_s"], dtype=float)
    longitude = np.asarray(arrays["lambda_rad"], dtype=float)
    modes = np.asarray(arrays["mode"], dtype=np.int8)
    indices = nearest_unique_frame_indices(
        time, config.visualization.animation_frame_step_s
    )
    positions = np.array(
        [
            rotating_coordinate_m(
                longitude[index],
                time[index],
                config.derived.reference_radius_m,
                config.derived.reference_mean_motion_rad_s,
            )
            for index in indices
        ]
    )
    centers = config.disturbance.center_s_m + config.disturbance.propagation_speed_m_s * (
        time[indices] - config.disturbance.start_time_s
    )
    margin = max(3.0 * config.derived.disturbance_sigma_m, 0.05 * np.ptp(positions))
    x_min = float(min(np.min(positions), np.min(centers)) - margin)
    x_max = float(max(np.max(positions), np.max(centers)) + margin)
    grid_m = np.linspace(x_min, x_max, 600)

    figure, (formation_axis, density_axis) = plt.subplots(
        2, 1, figsize=(10, 5.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )
    formation_axis.set_ylim(-1.0, 1.0)
    formation_axis.set_yticks([])
    formation_axis.set_ylabel("Satellite train")
    density_axis.set(xlabel="Rotating-frame along-track position (km)", ylabel="Density\nincrement (-)")
    density_axis.set_xlim(x_min / 1000.0, x_max / 1000.0)
    density_axis.set_ylim(
        min(0.0, 1.1 * config.disturbance.amplitude_fraction),
        max(0.05, 1.1 * config.disturbance.amplitude_fraction),
    )
    scatter = formation_axis.scatter([], [], s=55, edgecolors="black", linewidths=1.0)
    density_line, = density_axis.plot([], [], color="#9d174d", linewidth=2.0)
    center_line = density_axis.axvline(0.0, color="gray", linestyle="--", linewidth=1.0)
    time_text = formation_axis.text(0.02, 0.88, "", transform=formation_axis.transAxes)
    legend = [
        Line2D([0], [0], marker="o", color="none", markeredgecolor="black", markerfacecolor=color, label=label)
        for label, color in (
            ("L", "white"),
            ("L_TO_H", "#f4a261"),
            ("H", "#1f2937"),
            ("H_TO_L", "#2a9d8f"),
        )
    ]
    formation_axis.legend(handles=legend, loc="lower center", ncol=4, fontsize="small")
    mode_colors = np.array(["white", "#f4a261", "#1f2937", "#2a9d8f"])

    def update(frame_number: int):
        sample_index = int(indices[frame_number])
        sample_time = float(time[sample_index])
        sample_positions = positions[frame_number] / 1000.0
        scatter.set_offsets(np.column_stack((sample_positions, np.zeros_like(sample_positions))))
        scatter.set_facecolors(mode_colors[modes[sample_index]])
        increment = relative_density_increment(
            grid_m,
            sample_time,
            enabled=config.disturbance.enabled,
            amplitude_fraction=config.disturbance.amplitude_fraction,
            center_at_start_m=config.disturbance.center_s_m,
            start_time_s=config.disturbance.start_time_s,
            duration_s=config.disturbance.duration_s,
            sigma_m=config.derived.disturbance_sigma_m,
            propagation_speed_m_s=config.disturbance.propagation_speed_m_s,
        )
        density_line.set_data(grid_m / 1000.0, increment)
        center = config.disturbance.center_s_m + config.disturbance.propagation_speed_m_s * (
            sample_time - config.disturbance.start_time_s
        )
        center_line.set_xdata([center / 1000.0, center / 1000.0])
        time_text.set_text(f"Elapsed time: {sample_time / 3600.0:.3f} h")
        return scatter, density_line, center_line, time_text

    # Keep axes, labels, and legend in a cached background. Only the four changing
    # artists are redrawn for each frame before their RGBA buffer is piped to FFmpeg.
    # This is materially faster than asking Matplotlib to redraw a 4,321-frame
    # standard animation from scratch.
    for artist in (scatter, density_line, center_line, time_text):
        artist.set_animated(True)
    figure.set_dpi(80)
    figure.canvas.draw()
    background = figure.canvas.copy_from_bbox(figure.bbox)
    width, height = figure.canvas.get_width_height()
    writer = None
    try:
        imageio_ffmpeg.get_ffmpeg_exe()
        writer = imageio_ffmpeg.write_frames(
            str(destination),
            (width, height),
            pix_fmt_in="rgba",
            pix_fmt_out="yuv420p",
            fps=config.visualization.animation_fps,
            codec="libx264",
            quality=6,
            output_params=["-preset", "veryfast", "-movflags", "+faststart"],
            ffmpeg_log_level="warning",
        )
        writer.send(None)
        for frame_number in range(len(indices)):
            update(frame_number)
            figure.canvas.restore_region(background)
            formation_axis.draw_artist(scatter)
            formation_axis.draw_artist(time_text)
            density_axis.draw_artist(density_line)
            density_axis.draw_artist(center_line)
            figure.canvas.blit(figure.bbox)
            writer.send(np.asarray(figure.canvas.buffer_rgba(), dtype=np.uint8))
        writer.close()
        writer = None
    except Exception as exc:
        raise AnimationGenerationError(
            f"Failed to generate H.264 animation {destination}: {exc}"
        ) from exc
    finally:
        if writer is not None:
            writer.close()
        plt.close(figure)
    if not destination.is_file() or destination.stat().st_size == 0:
        raise AnimationGenerationError(f"FFmpeg did not create a usable file: {destination}")
    return destination
