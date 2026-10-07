"""The offline CLI must preserve inputs and reject an existing output directory."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from sat_string.io import save_case
from sat_string.simulator import simulate_standard_suite


def test_offline_cli_sensitivity_and_output_protection(tmp_path, short_config):
    suite = simulate_standard_suite(short_config)
    for name, result in suite.items():
        save_case(result, short_config, tmp_path / name)
    source = tmp_path / "controlled_disturbed"
    low = tmp_path / "all_low_drag_baseline"
    original_hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()}
    script = Path(__file__).parents[1] / "scripts/evaluate_saved_case.py"
    command = [sys.executable, str(script), "--run", str(source), "--all-low", str(low),
               "--output", str(tmp_path / "evaluation"), "--epsilon-e", "50", "1", "--hold-window", "5"]
    first = subprocess.run(command, capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    sensitivity = json.loads((tmp_path / "evaluation/settling_sensitivity.json").read_text())
    assert [item["settling"]["epsilon_e_m"] for item in sensitivity["evaluations"]] == [50, 1]
    before = (tmp_path / "evaluation/metrics.json").read_bytes()
    second = subprocess.run(command, capture_output=True, text=True)
    assert second.returncode != 0
    assert "Output directory already exists" in second.stderr
    assert (tmp_path / "evaluation/metrics.json").read_bytes() == before
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == expected for p, expected in original_hashes.items())


def test_offline_cli_rejects_mismatched_disturbance_reference(tmp_path, short_config):
    suite = simulate_standard_suite(short_config)
    for name, result in suite.items():
        save_case(result, short_config, tmp_path / name)
    script = Path(__file__).parents[1] / "scripts/evaluate_saved_case.py"
    result = subprocess.run([sys.executable, str(script), "--run", str(tmp_path / "controlled_undisturbed"),
                             "--all-low", str(tmp_path / "all_low_drag_baseline"), "--output", str(tmp_path / "evaluation")],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "matching disturbance enable flags" in result.stderr
    assert not (tmp_path / "evaluation").exists()
