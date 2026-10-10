import numpy as np
import pytest
from sat_string.necessity import exceedance_metrics, improvement, validate_comparison_config


def test_equality_is_allowed_and_terminal_has_no_duration():
    m = exceedance_metrics([0, 10, 20], [1, 1, 2], 1)
    assert m['ever_exceeded'] and m['first_exceedance_time_s'] == 20
    assert m['exceedance_duration_s'] == m['longest_continuous_exceedance_s'] == 0


def test_nonuniform_duration_and_separate_episodes():
    m = exceedance_metrics([0, 1, 5, 8, 10], [2, 0, 2, 2, 0], 1)
    assert m['exceedance_duration_s'] == 6
    assert m['exceedance_fraction'] == .6
    assert m['longest_continuous_exceedance_s'] == 5
    assert m['first_exceedance_time_s'] == 0


@pytest.mark.parametrize('values,expected', [([0, 0, 0], 0), ([2, 2, 2], 20)])
def test_all_inside_or_all_outside(values, expected):
    m = exceedance_metrics([0, 10, 20], values, 1)
    assert m['exceedance_duration_s'] == expected
    assert m['longest_continuous_exceedance_s'] == expected


@pytest.mark.parametrize('times,errors,tolerance', [
    ([0], [0], 1), ([0, 0], [1, 2], 1), ([0, 1], [0, np.nan], 1),
    ([0, 1], [-1, 2], 1), ([0, 1], [1, 2], -1), ([0, 1], [1], 1)])
def test_invalid_exceedance_input(times, errors, tolerance):
    with pytest.raises(ValueError): exceedance_metrics(times, errors, tolerance)


def test_improvement_preserves_absolute_values():
    assert improvement(100, 500) == dict(all_low_m=500, controlled_m=100,
        reduction_m=400, reduction_percent=80)
    assert improvement(2, 1)['reduction_m'] == -1
    assert improvement(0, 0)['reduction_percent'] is None


def test_config_rejects_physics_change():
    ref = dict(controller={}, derived={'controller_kp_s2_inv': 2}, physical={'mass': 4})
    cfg = dict(controller={'kp_s2_inv': 6}, derived={'controller_kp_s2_inv': 6}, physical={'mass': 4})
    validate_comparison_config(cfg, ref, 3)
    cfg['physical']['mass'] = 5
    with pytest.raises(ValueError): validate_comparison_config(cfg, ref, 3)
