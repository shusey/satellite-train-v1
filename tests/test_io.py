from pathlib import Path

import numpy as np
import pytest

from sat_string.io import (
    load_case_arrays,
    load_metadata,
    prepare_output_directory,
    save_case,
)
from sat_string.simulator import CaseDefinition, simulate_case


REQUIRED_ARRAYS = {
    "time_s",
    "a_m",
    "lambda_rad",
    "altitude_m",
    "area_m2",
    "mode",
    "density_kg_m3",
    "q_m_s2",
    "gap_error_m",
    "gap_rate_m_s",
    "gap_distance_m",
}


def test_case_round_trip_preserves_required_arrays(tmp_path, short_config):
    result = simulate_case(short_config, CaseDefinition("roundtrip", False, True))
    output = prepare_output_directory(tmp_path / "case")
    save_case(result, short_config, output)
    loaded = load_case_arrays(output)
    assert set(loaded) == REQUIRED_ARRAYS
    for key, expected in result.arrays().items():
        assert np.array_equal(loaded[key], expected)
    metadata = load_metadata(output)
    assert metadata["case"]["name"] == "roundtrip"
    assert metadata["recovery_time_resolution_s"] == 5.0
    assert (output / "switch_events.csv").is_file()
    assert (output / "normalized_config.yaml").is_file()


def test_existing_output_is_not_silently_overwritten(tmp_path):
    output = prepare_output_directory(tmp_path / "existing")
    marker = output / "marker.txt"
    marker.write_text("keep", encoding="utf-8")
    with pytest.raises(FileExistsError, match="--overwrite"):
        prepare_output_directory(output)
    assert marker.read_text(encoding="utf-8") == "keep"


def test_explicit_overwrite_replaces_output(tmp_path):
    output = prepare_output_directory(tmp_path / "replace")
    (output / "old.txt").write_text("old", encoding="utf-8")
    replaced = prepare_output_directory(output, overwrite=True)
    assert replaced.is_dir()
    assert not (replaced / "old.txt").exists()


def test_overwrite_refuses_workspace_root():
    with pytest.raises(ValueError, match="protected directory"):
        prepare_output_directory(Path.cwd(), overwrite=True)
