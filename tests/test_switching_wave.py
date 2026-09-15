from dataclasses import replace
import json

import numpy as np

from sat_string.attitude import SwitchEvent
from sat_string.switching_wave import (
    analyze_switching_wave,
    fit_distance_time,
    first_l_to_h_start_times,
    generate_aggregate_plots,
    save_sweep_csv,
    sweep_csv_row,
    threshold_scan_values,
)


def _event(satellite_id, start_time, from_mode="L", to_mode="H"):
    return SwitchEvent(
        satellite_id=satellite_id,
        command_time=start_time - 1.0,
        start_time=start_time,
        complete_time=start_time + 2.0,
        from_mode=from_mode,
        to_mode=to_mode,
    )


def test_threshold_grid_is_exactly_the_requested_union():
    values_in_nanounits = np.rint(threshold_scan_values() / 1.0e-9).astype(int)
    expected = sorted(set(range(100, 201, 10)) | set(range(130, 151, 2)))
    assert values_in_nanounits.tolist() == expected
    assert len(values_in_nanounits) == 19


def test_first_l_to_h_time_uses_start_and_ignores_unstarted_or_reverse_events():
    events = [
        _event(1, 20.0),
        _event(1, 10.0),
        _event(2, 12.0, from_mode="H", to_mode="L"),
        _event(3, 101.0),
    ]
    times = first_l_to_h_start_times(events, 3, simulation_end_time_s=100.0)
    assert times[0] == 10.0
    assert np.isnan(times[1])
    assert np.isnan(times[2])


def test_fit_returns_nan_below_three_points_and_full_diagnostics_otherwise():
    insufficient = fit_distance_time([1.0, 2.0], [10.0, 20.0])
    assert insufficient["point_count"] == 2
    assert np.isnan(insufficient["speed_m_s"])
    assert np.isnan(insufficient["r_squared"])
    assert np.isnan(insufficient["rmse_m"])

    fitted = fit_distance_time([1.0, 2.0, 3.0], [15.0, 25.0, 35.0])
    assert np.isclose(fitted["speed_m_s"], 10.0)
    assert np.isclose(fitted["intercept_m"], 5.0)
    assert np.isclose(fitted["r_squared"], 1.0)
    assert np.isclose(fitted["rmse_m"], 0.0, atol=1e-12)


def test_wave_analysis_excludes_three_nearest_satellites_and_computes_deltas(short_config):
    config = replace(
        short_config,
        formation=replace(short_config.formation, satellite_count=15),
        simulation=replace(short_config.simulation, duration_s=200.0),
    )
    center_id = 8
    events = [_event(center_id, 105.0)]
    for offset in range(1, 8):
        events.append(_event(center_id - offset, 100.0 + 10.0 * offset))
        events.append(_event(center_id + offset, 100.0 + 12.0 * offset))

    analysis = analyze_switching_wave(events, config, center_satellite_id=center_id)
    assert analysis["delay_s"] == 85.0
    assert analysis["left"]["fit_satellite_ids"] == [4, 3, 2, 1]
    assert analysis["right"]["fit_satellite_ids"] == [12, 13, 14, 15]
    assert analysis["left"]["fit"]["point_count"] == 4
    assert np.isclose(analysis["left"]["fit"]["speed_m_s"], 10_000.0)
    assert np.isclose(analysis["right"]["fit"]["speed_m_s"], 100_000.0 / 12.0)
    assert analysis["left"]["delta_t_s"] == [10.0] * 6
    assert analysis["right"]["delta_t_s"] == [12.0] * 6
    assert analysis["left"]["std_delta_t_s"] == 0.0
    assert analysis["left"]["cv_delta_t"] == 0.0


def test_partial_propagation_does_not_bridge_missing_satellite(short_config):
    config = replace(
        short_config,
        formation=replace(short_config.formation, satellite_count=11),
        simulation=replace(short_config.simulation, duration_s=100.0),
    )
    events = [_event(5, 30.0), _event(3, 70.0), _event(2, 80.0), _event(1, 90.0)]
    analysis = analyze_switching_wave(events, config, center_satellite_id=6)
    left = analysis["left"]
    assert left["fit"]["point_count"] == 2
    assert np.isnan(left["fit"]["speed_m_s"])
    assert left["delta_t_s"] == [10.0, 10.0]
    assert left["delta_t_records"][0]["inner_satellite_id"] == 3
    assert left["delta_t_records"][0]["outer_satellite_id"] == 2


def test_csv_contains_individual_deltas_and_aggregate_plots(tmp_path, short_config):
    events = [_event(2, 30.0), _event(1, 40.0), _event(4, 35.0), _event(5, 50.0)]
    analysis = analyze_switching_wave(events, short_config, center_satellite_id=3)
    row = sweep_csv_row(1.0e-7, 0.9e-7, analysis)
    csv_path = save_sweep_csv([row], tmp_path / "switching_wave_sweep.csv")
    assert csv_path.is_file()
    assert json.loads(row["delta_t_left"]) == [10.0]
    assert json.loads(row["delta_t_right"]) == [15.0]
    paths = generate_aggregate_plots([row], tmp_path, dpi=50)
    assert len(paths) == 5
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths)
