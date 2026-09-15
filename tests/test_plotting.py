from dataclasses import replace

import numpy as np

from sat_string.io import prepare_output_directory, save_case
from sat_string.plotting import (
    REQUIRED_PLOT_FILENAMES,
    generate_formation_animation,
    nearest_unique_frame_indices,
)
from sat_string.report import build_standard_report
from sat_string.simulator import CaseDefinition, simulate_case, simulate_standard_suite


def test_nearest_frames_do_not_duplicate_saved_samples():
    indices = nearest_unique_frame_indices(np.array([0.0, 10.0, 20.0]), 3.0)
    assert np.array_equal(indices, np.unique(indices))
    assert indices[0] == 0
    assert indices[-1] == 2


def test_report_regenerates_metrics_plots_and_h264_from_saved_data(tmp_path, short_config):
    run = prepare_output_directory(tmp_path / "suite")
    suite = simulate_standard_suite(short_config)
    for name, result in suite.items():
        save_case(result, short_config, run / name)
    aggregate = build_standard_report(run, short_config, dpi=72)
    assert "comparison" in aggregate
    assert (run / "metrics.json").is_file()
    assert (run / "comparison_results.npz").is_file()
    for filename in REQUIRED_PLOT_FILENAMES:
        assert (run / filename).stat().st_size > 0
    assert (run / "formation_animation.mp4").stat().st_size > 0


def test_moving_disturbance_animation_encodes_without_resimulation(tmp_path, short_config):
    moving_disturbance = replace(
        short_config.disturbance, propagation_speed_m_s=-1500.0
    )
    moving_config = replace(short_config, disturbance=moving_disturbance)
    result = simulate_case(
        moving_config,
        CaseDefinition("moving", False, True),
    )
    output = generate_formation_animation(
        result.arrays(), moving_config, tmp_path / "moving.mp4"
    )
    assert output.stat().st_size > 0
    density_peak_satellite = np.argmax(result.density_kg_m3, axis=1)
    assert density_peak_satellite.shape == result.time_s.shape
