"""Golden integrity and checker sensitivity; full rerun uses the separate CLI."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from sat_string.io import load_case_arrays
from sat_string.regression import compare_baseline
from sat_string.simulator import STANDARD_CASES


REFERENCE = Path(__file__).parent / "fixtures" / "regression_24h"


def test_frozen_reference_integrity():
    manifest = json.loads((REFERENCE / "reference.json").read_text(encoding="utf-8"))
    for name, expected in manifest["sha256"].items():
        assert hashlib.sha256((REFERENCE / name).read_bytes()).hexdigest() == expected


@pytest.fixture
def reference_copy(tmp_path):
    for case in STANDARD_CASES:
        shutil.copytree(REFERENCE / case.name, tmp_path / case.name)
        shutil.copyfile(REFERENCE / "config.yaml", tmp_path / case.name / "normalized_config.yaml")
    return tmp_path


def test_regression_reference_passes(reference_copy):
    assert compare_baseline(reference_copy, REFERENCE)["passed"]


def test_regression_detects_history_change(reference_copy):
    path = reference_copy / "controlled_disturbed"
    arrays = load_case_arrays(path)
    arrays["gap_error_m"][100, 0] += .01
    np.savez_compressed(path / "results.npz", **arrays)
    report = compare_baseline(reference_copy, REFERENCE)
    assert not report["passed"]
    assert not report["cases"]["controlled_disturbed"]["histories"]["gap_error_m"]["passed"]


def test_regression_detects_event_change(reference_copy):
    path = reference_copy / "controlled_disturbed" / "switch_events.csv"
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("20540.0", "20540.01"), encoding="utf-8")
    assert not compare_baseline(reference_copy, REFERENCE)["cases"]["controlled_disturbed"]["events_match"]
