from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from sat_string.config import config_from_dict, load_config


BASELINE = Path(__file__).parents[1] / "configs" / "baseline.yaml"


@pytest.fixture()
def baseline_dict():
    return yaml.safe_load(BASELINE.read_text(encoding="utf-8"))


def test_baseline_loads_and_computes_derived_values():
    config = load_config(BASELINE)
    assert config.formation.satellite_count == 21
    assert config.derived.reference_radius_m == pytest.approx(6.778137e6)
    assert config.derived.target_spacing_rad == pytest.approx(1.0e5 / 6.778137e6)
    assert config.derived.controller_kp_s2_inv == pytest.approx(1.0 / 21600.0**2)
    assert config.derived.controller_kd_s_inv == pytest.approx(2.0 / 21600.0)
    assert config.disturbance.center_s_m == 0.0


@pytest.mark.parametrize(
    ("section", "key", "value", "message"),
    [
        ("formation", "satellite_count", 1, "at least 2"),
        ("spacecraft", "area_high_m2", 0.005, "area_high_m2"),
        ("controller", "threshold_off_m_s2", 5e-7, "thresholds"),
        ("controller", "command_delay_s", -1.0, "nonnegative"),
        ("disturbance", "amplitude_fraction", -1.0, "greater than -1"),
        ("disturbance", "fwhm_m", 0.0, "positive"),
        ("simulation", "output_step_s", 0.5, "at least integration_step"),
        ("visualization", "animation_frame_step_s", 5.0, "at least the output step"),
    ],
)
def test_invalid_config_is_rejected(baseline_dict, section, key, value, message):
    data = deepcopy(baseline_dict)
    data[section][key] = value
    with pytest.raises(ValueError, match=message):
        config_from_dict(data)


def test_manual_disturbance_center_is_preserved(baseline_dict):
    data = deepcopy(baseline_dict)
    data["disturbance"]["center_mode"] = "manual"
    data["disturbance"]["center_s_m"] = 1234.5
    assert config_from_dict(data).disturbance.center_s_m == pytest.approx(1234.5)


def test_nonfinite_input_is_rejected(baseline_dict):
    data = deepcopy(baseline_dict)
    data["atmosphere"]["reference_density_kg_m3"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        config_from_dict(data)
