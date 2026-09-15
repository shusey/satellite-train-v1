import numpy as np
import pytest

from sat_string.integrator import integrate_interval, next_step_end, rk4_step


def test_rk4_fourth_order_accuracy():
    result = integrate_interval(
        lambda _t, y: y,
        0.0,
        np.array([1.0]),
        1.0,
        0.1,
    )
    assert result[0] == pytest.approx(np.e, rel=1e-6)


def test_rk4_reevaluates_stage_times_and_states():
    calls = []

    def rhs(time, state):
        calls.append((time, state.copy()))
        return np.array([time + state[0]])

    rk4_step(rhs, 2.0, np.array([3.0]), 0.4)
    assert [call[0] for call in calls] == pytest.approx([2.0, 2.2, 2.2, 2.4])
    assert not np.array_equal(calls[0][1], calls[1][1])
    assert not np.array_equal(calls[1][1], calls[2][1])


def test_noncommensurate_event_clips_step_exactly():
    assert next_step_end(0.0, 1.0, 10.0, [0.37, 2.1]) == pytest.approx(0.37)
    assert next_step_end(0.37, 1.0, 10.0, [0.37, 1.13]) == pytest.approx(1.13)


def test_integrate_interval_hits_partial_final_step():
    calls = []

    def rhs(time, state):
        calls.append(time)
        return np.ones_like(state)

    result = integrate_interval(rhs, 0.0, np.array([0.0]), 2.3, 1.0)
    assert result[0] == pytest.approx(2.3)
    assert max(calls) == pytest.approx(2.3)


def test_nonfinite_result_raises():
    with pytest.raises(FloatingPointError, match="NaN/Inf"):
        rk4_step(lambda _t, _y: np.array([np.inf]), 0.0, np.array([1.0]), 1.0)
