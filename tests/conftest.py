from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from sat_string.config import config_from_dict


@pytest.fixture()
def short_config():
    baseline_path = Path(__file__).parents[1] / "configs" / "baseline.yaml"
    data = deepcopy(yaml.safe_load(baseline_path.read_text(encoding="utf-8")))
    data["formation"]["satellite_count"] = 5
    data["controller"].update(
        {
            "design_time_s": 200.0,
            "control_period_s": 7.0,
            "command_delay_s": 1.3,
            "slew_time_s": 6.5,
            "dwell_time_s": 9.5,
        }
    )
    data["disturbance"].update(
        {
            "start_time_s": 20.0,
            "duration_s": 20.0,
            "fwhm_m": 100_000.0,
            "propagation_speed_m_s": 0.0,
        }
    )
    data["simulation"].update(
        {
            "duration_s": 60.0,
            "integration_step_s": 2.7,
            "output_step_s": 5.0,
        }
    )
    data["recovery"]["hold_time_s"] = 10.0
    data["visualization"].update(
        {
            "animation_enabled": True,
            "animation_frame_step_s": 10.0,
            "animation_fps": 5,
        }
    )
    return config_from_dict(data)

