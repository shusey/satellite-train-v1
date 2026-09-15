import numpy as np
import pytest

from sat_string.attitude import H, H_TO_L, L, L_TO_H, NO_TARGET, AttitudeStateMachine


def make_machine(delay=0.0):
    return AttitudeStateMachine(
        satellite_count=2,
        area_low_m2=0.01,
        area_high_m2=0.03,
        slew_time_s=120.0,
        dwell_time_s=300.0,
        command_delay_s=delay,
    )


def test_area_interpolation_endpoints_midpoint_and_zero_endpoint_slopes():
    machine = make_machine()
    machine.apply_commands(10.0, np.array([H, NO_TARGET]))
    assert machine.mode[0] == L_TO_H
    assert machine.areas(10.0)[0] == pytest.approx(0.01)
    assert machine.areas(70.0)[0] == pytest.approx(0.02)
    assert machine.areas(130.0)[0] == pytest.approx(0.03)
    epsilon = 1e-3
    start_slope = (machine.areas(10.0 + epsilon)[0] - machine.areas(10.0)[0]) / epsilon
    end_slope = (machine.areas(130.0)[0] - machine.areas(130.0 - epsilon)[0]) / epsilon
    assert abs(start_slope) < 1e-8
    assert abs(end_slope) < 1e-8


def test_no_recommand_during_transition_or_dwell():
    machine = make_machine()
    machine.apply_commands(0.0, np.array([H, NO_TARGET]))
    assert not machine.command_allowed(60.0)[0]
    machine.process_completions(120.0)
    assert machine.mode[0] == H
    assert not machine.command_allowed(419.999)[0]
    assert machine.command_allowed(420.0)[0]
    machine.apply_commands(420.0, np.array([L, NO_TARGET]))
    assert machine.mode[0] == H_TO_L


def test_delayed_command_is_pending_and_display_mode_stays_stable():
    machine = make_machine(delay=15.0)
    machine.apply_commands(5.0, np.array([H, NO_TARGET]))
    assert machine.mode[0] == L
    assert not machine.command_allowed(10.0)[0]
    assert machine.process_pending_starts(19.999) == []
    assert machine.process_pending_starts(20.0) == [0]
    assert machine.mode[0] == L_TO_H
    assert machine.events[0].command_time == 5.0
    assert machine.events[0].start_time == 20.0
    assert machine.events[0].complete_time == 140.0


def test_next_event_reports_pending_start_then_completion():
    machine = make_machine(delay=5.5)
    machine.apply_commands(1.0, np.array([H, NO_TARGET]))
    assert machine.next_event_time(1.0) == pytest.approx(6.5)
    machine.process_pending_starts(6.5)
    assert machine.next_event_time(6.5) == pytest.approx(126.5)


def test_transition_event_is_counted_once():
    machine = make_machine()
    machine.apply_commands(0.0, np.array([H, NO_TARGET]))
    machine.process_completions(120.0)
    assert len(machine.events) == 1
    assert machine.events[0].from_mode == "L"
    assert machine.events[0].to_mode == "H"

