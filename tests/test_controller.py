import numpy as np
import pytest

from sat_string.controller import (
    HIGH_MODE,
    LOW_MODE,
    NO_COMMAND,
    select_commands,
    switching_function,
)


def test_open_path_switching_function_signs():
    error = np.array([1.0, 3.0, -2.0])
    rate = np.array([0.5, 0.0, -0.5])
    kp, kd = 2.0, 4.0
    edge = kp * error + kd * rate
    q = switching_function(error, rate, kp, kd)
    assert q == pytest.approx([edge[0], edge[1] - edge[0], edge[2] - edge[1], -edge[2]])


def test_hysteresis_includes_exact_thresholds():
    commands = select_commands(
        np.array([4e-7, 2e-7, 3e-7, 1e-7]),
        np.array([LOW_MODE, HIGH_MODE, LOW_MODE, HIGH_MODE]),
        np.ones(4, dtype=bool),
        4e-7,
        2e-7,
    )
    assert np.array_equal(
        commands,
        np.array([HIGH_MODE, LOW_MODE, NO_COMMAND, LOW_MODE], dtype=np.int8),
    )


def test_disallowed_satellites_receive_no_command():
    commands = select_commands(
        np.array([1.0, -1.0]),
        np.array([LOW_MODE, HIGH_MODE]),
        np.array([False, False]),
        4e-7,
        2e-7,
    )
    assert np.all(commands == NO_COMMAND)


def test_command_selection_is_scan_order_independent():
    q = np.array([1e-6, -1e-6, 5e-7, 0.0])
    modes = np.array([LOW_MODE, HIGH_MODE, LOW_MODE, HIGH_MODE])
    allowed = np.array([True, True, False, True])
    forward = select_commands(q, modes, allowed, 4e-7, 2e-7)
    reverse = select_commands(q[::-1], modes[::-1], allowed[::-1], 4e-7, 2e-7)[::-1]
    assert np.array_equal(forward, reverse)
