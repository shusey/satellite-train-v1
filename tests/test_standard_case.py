from dataclasses import replace

import numpy as np

from sat_string.simulator import CaseDefinition, simulate_case


def test_short_case_step_size_convergence(short_config):
    one_second = replace(
        short_config,
        simulation=replace(short_config.simulation, integration_step_s=1.0),
    )
    half_second = replace(
        short_config,
        simulation=replace(short_config.simulation, integration_step_s=0.5),
    )
    case = CaseDefinition("convergence", True, True)
    coarse = simulate_case(one_second, case)
    fine = simulate_case(half_second, case)
    coarse_max_error = float(np.max(np.abs(coarse.gap_error_m)))
    fine_max_error = float(np.max(np.abs(fine.gap_error_m)))
    coarse_final_altitude = float(np.mean(coarse.altitude_m[-1]))
    fine_final_altitude = float(np.mean(fine.altitude_m[-1]))
    assert np.isclose(coarse_max_error, fine_max_error, rtol=0.01, atol=1.0)
    assert np.isclose(coarse_final_altitude, fine_final_altitude, rtol=0.01, atol=1.0)


def test_required_array_shapes_and_finiteness(short_config):
    result = simulate_case(
        short_config, CaseDefinition("shape_regression", True, True)
    )
    sample_count = result.time_s.size
    satellite_count = short_config.formation.satellite_count
    for name in (
        "a_m",
        "lambda_rad",
        "altitude_m",
        "area_m2",
        "mode",
        "density_kg_m3",
        "q_m_s2",
    ):
        value = getattr(result, name)
        assert value.shape == (sample_count, satellite_count)
        assert np.all(np.isfinite(value))
    for name in ("gap_error_m", "gap_rate_m_s", "gap_distance_m"):
        value = getattr(result, name)
        assert value.shape == (sample_count, satellite_count - 1)
        assert np.all(np.isfinite(value))

