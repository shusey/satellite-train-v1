"""Configuration loading, normalization, validation, and derived quantities."""

from __future__ import annotations

from dataclasses import MISSING, asdict, dataclass, fields
from math import isfinite, log, pi, sqrt
from pathlib import Path
from typing import Any, Mapping, get_args, get_origin, get_type_hints
import warnings

import yaml


@dataclass(frozen=True)
class ConstantsConfig:
    earth_radius_m: float
    mu_m3_s2: float


@dataclass(frozen=True)
class SpacecraftConfig:
    mass_kg: float
    drag_coefficient: float
    area_low_m2: float
    area_high_m2: float


@dataclass(frozen=True)
class FormationConfig:
    satellite_count: int
    reference_altitude_m: float
    target_spacing_m: float
    topology: str


@dataclass(frozen=True)
class ControllerConfig:
    design_time_s: float
    threshold_on_m_s2: float
    threshold_off_m_s2: float
    control_period_s: float
    command_delay_s: float
    slew_time_s: float
    dwell_time_s: float
    reference_on_error_m: float | None = None
    reference_on_rate_m_s: float | None = None
    reference_off_error_m: float | None = None
    reference_off_rate_m_s: float | None = None


@dataclass(frozen=True)
class AtmosphereConfig:
    reference_density_kg_m3: float
    effective_scale_height_m: float


@dataclass(frozen=True)
class DisturbanceConfig:
    enabled: bool
    type: str
    start_time_s: float
    duration_s: float
    amplitude_fraction: float
    center_mode: str
    center_s_m: float
    fwhm_m: float
    propagation_speed_m_s: float


@dataclass(frozen=True)
class IntegrationConfig:
    duration_s: float
    integration_step_s: float
    output_step_s: float
    integrator: str


@dataclass(frozen=True)
class RecoveryConfig:
    position_tolerance_m: float
    velocity_tolerance_m_s: float
    hold_time_s: float


@dataclass(frozen=True)
class VisualizationConfig:
    animation_enabled: bool
    animation_frame_step_s: float
    animation_fps: int


@dataclass(frozen=True)
class DerivedConfig:
    reference_radius_m: float
    reference_mean_motion_rad_s: float
    reference_period_s: float
    target_spacing_rad: float
    controller_kp_s2_inv: float
    controller_kd_s_inv: float
    disturbance_sigma_m: float


@dataclass(frozen=True)
class SimulationConfig:
    constants: ConstantsConfig
    spacecraft: SpacecraftConfig
    formation: FormationConfig
    controller: ControllerConfig
    atmosphere: AtmosphereConfig
    disturbance: DisturbanceConfig
    simulation: IntegrationConfig
    recovery: RecoveryConfig
    visualization: VisualizationConfig
    derived: DerivedConfig

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON/YAML-serializable normalized configuration."""
        return asdict(self)


_SECTIONS: dict[str, type[Any]] = {
    "constants": ConstantsConfig,
    "spacecraft": SpacecraftConfig,
    "formation": FormationConfig,
    "controller": ControllerConfig,
    "atmosphere": AtmosphereConfig,
    "disturbance": DisturbanceConfig,
    "simulation": IntegrationConfig,
    "recovery": RecoveryConfig,
    "visualization": VisualizationConfig,
}


def _section(data: Mapping[str, Any], name: str) -> Any:
    if name not in data or not isinstance(data[name], Mapping):
        raise ValueError(f"Missing or invalid configuration section: {name}")
    cls = _SECTIONS[name]
    hints = get_type_hints(cls)
    converted: dict[str, Any] = {}
    try:
        known = {field.name for field in fields(cls)}
        unknown = set(data[name]) - known
        if unknown:
            raise TypeError(f"unexpected fields: {', '.join(sorted(unknown))}")
        for field in fields(cls):
            if field.name not in data[name]:
                if field.default is MISSING and field.default_factory is MISSING:
                    raise TypeError(f"missing required field: {field.name}")
                continue
            value = data[name][field.name]
            expected = hints[field.name]
            args = get_args(expected) if get_origin(expected) is not None else ()
            nullable = type(None) in args
            base = next((arg for arg in args if arg is not type(None)), expected)
            if value is None and nullable:
                converted[field.name] = None
            elif base is float:
                if isinstance(value, bool):
                    raise TypeError(f"{field.name} must be numeric")
                converted[field.name] = float(value)
            elif base is int:
                if isinstance(value, bool) or float(value) != int(float(value)):
                    raise TypeError(f"{field.name} must be an integer")
                converted[field.name] = int(float(value))
            elif base is bool:
                if not isinstance(value, bool):
                    raise TypeError(f"{field.name} must be boolean")
                converted[field.name] = value
            elif base is str:
                if not isinstance(value, str):
                    raise TypeError(f"{field.name} must be a string")
                converted[field.name] = value
            else:
                converted[field.name] = value
        return cls(**converted)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid fields in configuration section {name}: {exc}") from exc


def _validate_finite(config: SimulationConfig) -> None:
    def walk(value: Any, path: str) -> None:
        if isinstance(value, bool) or value is None or isinstance(value, str):
            return
        if isinstance(value, Mapping):
            for key, item in value.items():
                walk(item, f"{path}.{key}" if path else key)
            return
        if isinstance(value, (int, float)) and not isfinite(float(value)):
            raise ValueError(f"Configuration value must be finite: {path}={value!r}")

    walk(config.to_dict(), "")


def _validate(config: SimulationConfig, raw: Mapping[str, Any]) -> None:
    c = config
    if c.formation.satellite_count < 2:
        raise ValueError("formation.satellite_count must be at least 2")
    if c.formation.topology != "open_path":
        raise ValueError("v1 supports only formation.topology=open_path")
    if not (c.spacecraft.area_high_m2 > c.spacecraft.area_low_m2 > 0.0):
        raise ValueError("spacecraft areas must satisfy area_high_m2 > area_low_m2 > 0")
    if c.spacecraft.mass_kg <= 0.0 or c.spacecraft.drag_coefficient <= 0.0:
        raise ValueError("mass and drag coefficient must be positive")
    if not (
        c.controller.threshold_on_m_s2
        > c.controller.threshold_off_m_s2
        >= 0.0
    ):
        raise ValueError("controller thresholds must satisfy on > off >= 0")
    positive_times = {
        "controller.design_time_s": c.controller.design_time_s,
        "controller.control_period_s": c.controller.control_period_s,
        "controller.slew_time_s": c.controller.slew_time_s,
        "controller.dwell_time_s": c.controller.dwell_time_s,
        "disturbance.duration_s": c.disturbance.duration_s,
        "simulation.duration_s": c.simulation.duration_s,
        "simulation.integration_step_s": c.simulation.integration_step_s,
        "simulation.output_step_s": c.simulation.output_step_s,
        "recovery.hold_time_s": c.recovery.hold_time_s,
        "visualization.animation_frame_step_s": c.visualization.animation_frame_step_s,
    }
    for name, value in positive_times.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    if c.controller.command_delay_s < 0.0:
        raise ValueError("controller.command_delay_s must be nonnegative")
    if c.disturbance.start_time_s < 0.0:
        raise ValueError("disturbance.start_time_s must be nonnegative")
    if c.recovery.position_tolerance_m < 0.0 or c.recovery.velocity_tolerance_m_s < 0.0:
        raise ValueError("recovery tolerances must be nonnegative")
    if c.constants.earth_radius_m <= 0.0 or c.constants.mu_m3_s2 <= 0.0:
        raise ValueError("physical constants must be positive")
    if c.formation.reference_altitude_m <= 0.0 or c.formation.target_spacing_m <= 0.0:
        raise ValueError("reference altitude and target spacing must be positive")
    if c.atmosphere.reference_density_kg_m3 <= 0.0 or c.atmosphere.effective_scale_height_m <= 0.0:
        raise ValueError("atmosphere density and scale height must be positive")
    if c.disturbance.amplitude_fraction <= -1.0:
        raise ValueError("disturbance.amplitude_fraction must be greater than -1")
    if c.disturbance.fwhm_m <= 0.0:
        raise ValueError("disturbance.fwhm_m must be positive")
    if c.disturbance.type != "single_pulse":
        raise ValueError("v1 supports only disturbance.type=single_pulse")
    if c.disturbance.center_mode not in {"center_of_train", "manual"}:
        raise ValueError("disturbance.center_mode must be center_of_train or manual")
    if c.simulation.integrator.lower() != "rk4":
        raise ValueError("v1 supports only simulation.integrator=rk4")
    if c.simulation.output_step_s < c.simulation.integration_step_s:
        raise ValueError("simulation.output_step_s must be at least integration_step_s")
    if c.visualization.animation_frame_step_s < c.simulation.output_step_s:
        raise ValueError("animation frame step must be at least the output step")
    if c.simulation.duration_s <= c.disturbance.start_time_s + c.disturbance.duration_s:
        raise ValueError("simulation must end after the disturbance ends")
    if c.visualization.animation_fps <= 0:
        raise ValueError("visualization.animation_fps must be positive")

    optional = raw.get("controller", {})
    pairs = (
        ("reference_on_error_m", "reference_on_rate_m_s", c.controller.threshold_on_m_s2),
        ("reference_off_error_m", "reference_off_rate_m_s", c.controller.threshold_off_m_s2),
    )
    for e_key, v_key, expected in pairs:
        if optional.get(e_key) is not None and optional.get(v_key) is not None:
            calculated = c.derived.controller_kp_s2_inv * float(optional[e_key]) + c.derived.controller_kd_s_inv * float(optional[v_key])
            tolerance = max(1e-12, 0.05 * abs(expected))
            if abs(calculated - expected) > tolerance:
                warnings.warn(
                    f"Controller reference values imply {calculated:.6g}, not threshold {expected:.6g}",
                    UserWarning,
                    stacklevel=2,
                )
    _validate_finite(c)


def config_from_dict(data: Mapping[str, Any]) -> SimulationConfig:
    """Build and validate a normalized configuration from a mapping."""
    sections = {name: _section(data, name) for name in _SECTIONS}
    for section_name, section in sections.items():
        for key, value in asdict(section).items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                if not isfinite(float(value)):
                    raise ValueError(
                        f"Configuration value must be finite: {section_name}.{key}={value!r}"
                    )
    if sections["disturbance"].center_mode == "center_of_train":
        disturbance_dict = asdict(sections["disturbance"])
        disturbance_dict["center_s_m"] = 0.0
        sections["disturbance"] = DisturbanceConfig(**disturbance_dict)

    reference_radius = sections["constants"].earth_radius_m + sections["formation"].reference_altitude_m
    mean_motion = sqrt(sections["constants"].mu_m3_s2 / reference_radius**3)
    design_time = sections["controller"].design_time_s
    derived = DerivedConfig(
        reference_radius_m=reference_radius,
        reference_mean_motion_rad_s=mean_motion,
        reference_period_s=2.0 * pi / mean_motion,
        target_spacing_rad=sections["formation"].target_spacing_m / reference_radius,
        controller_kp_s2_inv=1.0 / design_time**2,
        controller_kd_s_inv=2.0 / design_time,
        disturbance_sigma_m=sections["disturbance"].fwhm_m / (2.0 * sqrt(2.0 * log(2.0))),
    )
    config = SimulationConfig(**sections, derived=derived)
    _validate(config, data)
    return config


def load_config(path: str | Path) -> SimulationConfig:
    """Load a YAML configuration file and return its normalized representation."""
    source = Path(path)
    with source.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, Mapping):
        raise ValueError("The YAML document must contain a mapping at its root")
    return config_from_dict(data)


def dump_normalized_config(config: SimulationConfig, path: str | Path) -> None:
    """Write the exact normalized configuration used by a run."""
    destination = Path(path)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        yaml.safe_dump(config.to_dict(), stream, sort_keys=False, allow_unicode=True)
