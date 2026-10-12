"""Reproduce literature-review conversions, without running orbital simulations.

This checks definitions and arithmetic. It does not validate a population sigma
or convert orbit-determination uncertainty into physical formation residuals.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import yaml


def calculate() -> dict:
    here = Path(__file__).resolve().parent
    repo = here.parents[1]
    config = yaml.safe_load((repo / "configs/baseline.yaml").read_text(encoding="utf-8"))
    evidence = json.loads((here / "evidence.json").read_text(encoding="utf-8"))
    conditions = json.loads((here / "recommended_conditions.json").read_text(encoding="utf-8"))
    ids = [row["id"] for row in evidence["references"]]
    assert len(ids) == len(set(ids))
    for row in evidence["references"]:
        assert row["url"].startswith("https://")
        assert row["location_checked"] and row["quantity_type"] and row["limitations"]
    for field in ("references_a", "references_b", "references_requirements"):
        assert set(conditions[field]).issubset(ids)

    sc = config["spacecraft"]
    constants = config["constants"]
    # PyYAML treats some scientific-notation scalars as strings; parse the
    # numeric settings explicitly, as the simulator's configuration loader does.
    r = float(constants["earth_radius_m"]) + float(config["formation"]["reference_altitude_m"])
    mu = float(constants["mu_m3_s2"])
    n = math.sqrt(mu / r**3)
    speed = math.sqrt(mu / r)
    b = sc["drag_coefficient"] * sc["area_low_m2"] / sc["mass_kg"]
    assert math.isclose(b, 0.0055, abs_tol=1e-15)
    beta = 1 / b
    # Exact inverse-BC conversion: beta_i smaller by 1% gives B_i bigger
    # by 1/0.99-1, rather than exactly 1%.
    beta_i, beta_j = beta * 0.99, beta
    relative_b = beta_j / beta_i - 1
    assert math.isclose(relative_b, ((1/beta_i) - (1/beta_j)) / (1/beta_j))
    assert relative_b > 0.01
    mass_examples = []
    for fraction in (0.001, 0.003, 0.01):
        mass_examples.append({"relative_difference": fraction,
                              "first_order_mass_difference_g": sc["mass_kg"] * fraction * 1000,
                              "exact_B_change_for_heavier_satellite": 1/(1+fraction)-1})
    angle_examples = []
    side_to_low_ratio = sc["area_high_m2"] / sc["area_low_m2"]
    for fraction in (0.003, 0.01):
        angle = fraction / side_to_low_ratio
        projected_ratio = math.cos(angle) + side_to_low_ratio * abs(math.sin(angle))
        assert abs((projected_ratio-1) - fraction) < 1e-5
        angle_examples.append({"area_fraction": fraction,
                               "angle_deg_first_order": math.degrees(angle),
                               "exact_area_fraction_at_angle": projected_ratio-1})
    speed_705 = math.sqrt(mu / (float(constants["earth_radius_m"]) + 705000))
    # Centering preserves differences. Use a deterministic vector to check it.
    profile = [-0.7, 0.1, 1.3, -0.2]
    centered = [v - sum(profile)/len(profile) for v in profile]
    assert math.isclose(sum(centered), 0, abs_tol=1e-14)
    assert all(math.isclose(profile[i+1]-profile[i], centered[i+1]-centered[i])
               for i in range(len(profile)-1))
    report = (here / "README.md").read_text(encoding="utf-8")
    summary = report.split("# 詳細技術付録", 1)[0]
    assert 1500 <= len(summary) <= 3200, len(summary)
    return {
        "type": "arithmetic_and_evidence_ledger_audit_not_new_physics_run",
        "reference_count": len(ids), "summary_characters_including_markdown": len(summary),
        "B_definition": "Cd*A/m [m^2/kg]", "B0_m2_kg": b,
        "inverse_B0_kg_m2": beta, "inverse_B_conversion_example": relative_b,
        "mass_difference_examples": mass_examples,
        "low_face_rectangular_geometry_examples": angle_examples,
        "baseline_mean_motion_rad_s": n, "baseline_speed_m_s": speed,
        "drift_per_day_for_pair_delta_a_0p1m_m": 1.5*n*0.1*86400,
        "typical_pair_drift_per_day_for_generated_sigma_a_0p1m_m": 1.5*n*math.sqrt(2)*0.1*86400,
        "typical_pair_sigma_for_generated_B_sigma_0p003": math.sqrt(2)*0.003,
        "generated_B_sigma_for_pair_sigma_0p003": 0.003/math.sqrt(2),
        "atrain_orbital_distance_conversion_m": {"target": speed_705*12.5,
                                                 "tolerance": speed_705*2.5},
        "requirement_sensitivity": [{"error_m": e, "percent_of_100km": e/100000*100,
                                     "orbital_time_difference_s": e/speed}
                                    for e in (10,100,500,1000,5000)],
        "TDX_Tsx_B_difference_from_2024_masses": 1293/1256-1,
        "ANSER_max_flap_area_difference_relative_FLW2": 1800/2100-1,
        "validated": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write a new audit JSON; refuse to replace an existing file")
    args = parser.parse_args()
    result = calculate()
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(text)
    print(text)
