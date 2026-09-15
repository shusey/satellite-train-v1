"""Run-directory and simulation artifact input/output."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
from typing import Any, Mapping

import numpy as np

from .config import SimulationConfig, dump_normalized_config
from .simulator import SimulationResult


def prepare_output_directory(path: str | Path, *, overwrite: bool = False) -> Path:
    """Create a clean output directory without silently replacing existing data."""
    destination = Path(path).resolve()
    if destination.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output directory already exists: {destination}. Use --overwrite to replace it."
            )
        protected = {Path.cwd().resolve(), Path.home().resolve()}
        if destination.parent == destination or destination in protected:
            raise ValueError(f"Refusing to overwrite protected directory: {destination}")
        if destination.is_dir() and (destination / ".git").exists():
            raise ValueError(f"Refusing to overwrite Git repository root: {destination}")
        if destination.is_dir():
            shutil.rmtree(destination)
        else:
            destination.unlink()
    destination.mkdir(parents=True)
    return destination


def write_json(path: str | Path, data: Mapping[str, Any]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def save_case(
    result: SimulationResult,
    config: SimulationConfig,
    output_directory: str | Path,
) -> Path:
    """Save required arrays, normalized metadata, and switch events for one case."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination / "results.npz", **result.arrays())
    metadata = {
        "format_version": 1,
        "case": {
            "name": result.case.name,
            "control_enabled": result.case.control_enabled,
            "disturbance_enabled": result.case.disturbance_enabled,
        },
        "configuration": config.to_dict(),
        "recovery_time_resolution_s": config.simulation.output_step_s,
        "array_shapes": {
            key: list(value.shape) for key, value in result.arrays().items()
        },
    }
    write_json(destination / "metadata.json", metadata)
    dump_normalized_config(config, destination / "normalized_config.yaml")
    fieldnames = [
        "satellite_id",
        "command_time",
        "start_time",
        "complete_time",
        "from_mode",
        "to_mode",
    ]
    with (destination / "switch_events.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(event.to_dict() for event in result.switch_events)
    return destination


def load_case_arrays(directory: str | Path) -> dict[str, np.ndarray]:
    """Load all stored arrays into memory without running a simulation."""
    with np.load(Path(directory) / "results.npz", allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def load_metadata(directory: str | Path) -> dict[str, Any]:
    with (Path(directory) / "metadata.json").open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_switch_events(directory: str | Path) -> list[dict[str, Any]]:
    """Load switch event records with numeric time fields restored."""
    path = Path(directory) / "switch_events.csv"
    events: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            events.append(
                {
                    "satellite_id": int(row["satellite_id"]),
                    "command_time": float(row["command_time"]),
                    "start_time": float(row["start_time"]),
                    "complete_time": float(row["complete_time"]),
                    "from_mode": row["from_mode"],
                    "to_mode": row["to_mode"],
                }
            )
    return events
