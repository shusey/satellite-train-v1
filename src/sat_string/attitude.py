"""Finite-slew binary attitude state machine with delay, hysteresis, and dwell."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


L = np.int8(0)
L_TO_H = np.int8(1)
H = np.int8(2)
H_TO_L = np.int8(3)
NO_TARGET = np.int8(-1)
MODE_NAMES = {0: "L", 1: "L_TO_H", 2: "H", 3: "H_TO_L"}


def smoothstep(progress: ArrayLike) -> NDArray[np.float64]:
    """C1-continuous cubic interpolation on [0, 1]."""
    clipped = np.clip(np.asarray(progress, dtype=float), 0.0, 1.0)
    return 3.0 * clipped**2 - 2.0 * clipped**3


@dataclass(frozen=True)
class SwitchEvent:
    satellite_id: int
    command_time: float
    start_time: float
    complete_time: float
    from_mode: str
    to_mode: str

    def to_dict(self) -> dict[str, int | float | str]:
        return asdict(self)


class AttitudeStateMachine:
    """Vectorized state of all satellite attitude actuators."""

    def __init__(
        self,
        satellite_count: int,
        area_low_m2: float,
        area_high_m2: float,
        slew_time_s: float,
        dwell_time_s: float,
        command_delay_s: float = 0.0,
    ) -> None:
        self.area_low_m2 = float(area_low_m2)
        self.area_high_m2 = float(area_high_m2)
        self.slew_time_s = float(slew_time_s)
        self.dwell_time_s = float(dwell_time_s)
        self.command_delay_s = float(command_delay_s)
        self.mode = np.full(satellite_count, L, dtype=np.int8)
        self.transition_start_time = np.full(satellite_count, np.nan)
        self.transition_complete_time = np.full(satellite_count, np.nan)
        self.next_allowed_time = np.zeros(satellite_count)
        self.pending_target = np.full(satellite_count, NO_TARGET, dtype=np.int8)
        self.pending_start_time = np.full(satellite_count, np.nan)
        self.events: list[SwitchEvent] = []

    @property
    def satellite_count(self) -> int:
        return int(self.mode.size)

    def command_allowed(self, time_s: float, atol: float = 1e-9) -> NDArray[np.bool_]:
        stable = (self.mode == L) | (self.mode == H)
        no_pending = self.pending_target == NO_TARGET
        return stable & no_pending & (time_s >= self.next_allowed_time - atol)

    def areas(self, time_s: float) -> NDArray[np.float64]:
        """Return actual projected areas, including smooth finite slews."""
        area = np.where(self.mode == H, self.area_high_m2, self.area_low_m2).astype(float)
        transitioning = (self.mode == L_TO_H) | (self.mode == H_TO_L)
        if np.any(transitioning):
            progress = (
                time_s - self.transition_start_time[transitioning]
            ) / self.slew_time_s
            blend = smoothstep(progress)
            rising = self.mode[transitioning] == L_TO_H
            transition_area = np.where(
                rising,
                self.area_low_m2
                + (self.area_high_m2 - self.area_low_m2) * blend,
                self.area_high_m2
                - (self.area_high_m2 - self.area_low_m2) * blend,
            )
            area[transitioning] = transition_area
        return area

    def _start_transition(self, index: int, target: int, start_time_s: float) -> None:
        if target == H and self.mode[index] == L:
            self.mode[index] = L_TO_H
        elif target == L and self.mode[index] == H:
            self.mode[index] = H_TO_L
        else:
            raise RuntimeError(
                f"Invalid transition for satellite {index + 1}: mode={self.mode[index]}, target={target}"
            )
        self.transition_start_time[index] = start_time_s
        self.transition_complete_time[index] = start_time_s + self.slew_time_s

    def apply_commands(self, time_s: float, targets: ArrayLike) -> None:
        """Register a simultaneously selected command vector."""
        target_array = np.asarray(targets, dtype=np.int8)
        if target_array.shape != self.mode.shape:
            raise ValueError("target command vector has the wrong shape")
        allowed = self.command_allowed(time_s)
        requested = target_array != NO_TARGET
        if np.any(requested & ~allowed):
            indices = np.flatnonzero(requested & ~allowed) + 1
            raise RuntimeError(f"Commands are not allowed for satellites {indices.tolist()}")

        for index in np.flatnonzero(requested):
            target = int(target_array[index])
            stable_mode = int(self.mode[index])
            if target not in (int(L), int(H)) or target == stable_mode:
                raise ValueError(f"Invalid target mode {target} for satellite {index + 1}")
            start_time = time_s + self.command_delay_s
            self.events.append(
                SwitchEvent(
                    satellite_id=index + 1,
                    command_time=float(time_s),
                    start_time=float(start_time),
                    complete_time=float(start_time + self.slew_time_s),
                    from_mode=MODE_NAMES[stable_mode],
                    to_mode=MODE_NAMES[target],
                )
            )
            if self.command_delay_s == 0.0:
                self._start_transition(index, target, start_time)
            else:
                self.pending_target[index] = target
                self.pending_start_time[index] = start_time

    def process_pending_starts(self, time_s: float, atol: float = 1e-9) -> list[int]:
        """Apply all delayed starts due at this event time."""
        due = (self.pending_target != NO_TARGET) & (
            self.pending_start_time <= time_s + atol
        )
        started: list[int] = []
        for index in np.flatnonzero(due):
            target = int(self.pending_target[index])
            start_time = float(self.pending_start_time[index])
            self.pending_target[index] = NO_TARGET
            self.pending_start_time[index] = np.nan
            self._start_transition(index, target, start_time)
            started.append(index)
        return started

    def process_completions(self, time_s: float, atol: float = 1e-9) -> list[int]:
        """Finalize all slews due at this event time and begin their dwell clocks."""
        transitioning = (self.mode == L_TO_H) | (self.mode == H_TO_L)
        due = transitioning & (self.transition_complete_time <= time_s + atol)
        completed: list[int] = []
        for index in np.flatnonzero(due):
            completion = float(self.transition_complete_time[index])
            self.mode[index] = H if self.mode[index] == L_TO_H else L
            self.transition_start_time[index] = np.nan
            self.transition_complete_time[index] = np.nan
            self.next_allowed_time[index] = completion + self.dwell_time_s
            completed.append(index)
        return completed

    def next_event_time(self, after_time_s: float, atol: float = 1e-9) -> float:
        """Return the next pending-start or completion time, or infinity."""
        candidates = np.concatenate(
            (
                self.pending_start_time[np.isfinite(self.pending_start_time)],
                self.transition_complete_time[np.isfinite(self.transition_complete_time)],
            )
        )
        future = candidates[candidates > after_time_s + atol]
        return float(np.min(future)) if future.size else float("inf")

