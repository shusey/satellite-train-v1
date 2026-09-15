from pathlib import Path
import subprocess
import sys

import yaml


PROJECT_ROOT = Path(__file__).parents[1]


def test_run_case_cli_and_output_protection(tmp_path, short_config):
    config_path = tmp_path / "short.yaml"
    config_path.write_text(
        yaml.safe_dump(short_config.to_dict(), sort_keys=False), encoding="utf-8"
    )
    output = tmp_path / "cli_case"
    command = [
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "run_case.py"),
        "--config",
        str(config_path),
        "--output",
        str(output),
    ]
    first = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert first.returncode == 0, first.stderr
    assert (output / "results.npz").is_file()
    assert (output / "metrics.json").is_file()
    assert "estimated_remaining" in (output / "run.log").read_text(encoding="utf-8")

    second = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert second.returncode != 0
    assert "--overwrite" in second.stderr


def test_standard_suite_and_report_cli_use_saved_histories(tmp_path, short_config):
    config_path = tmp_path / "suite_short.yaml"
    config_path.write_text(
        yaml.safe_dump(short_config.to_dict(), sort_keys=False), encoding="utf-8"
    )
    output = tmp_path / "cli_suite"
    suite_command = [
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "run_standard_suite.py"),
        "--config",
        str(config_path),
        "--output",
        str(output),
        "--skip-animation",
    ]
    suite_run = subprocess.run(
        suite_command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert suite_run.returncode == 0, suite_run.stderr
    result_files = [
        output / name / "results.npz"
        for name in (
            "controlled_disturbed",
            "controlled_undisturbed",
            "all_low_drag_baseline",
        )
    ]
    assert all(path.is_file() for path in result_files)
    mtimes = [path.stat().st_mtime_ns for path in result_files]
    controlled_log = (
        output / "controlled_disturbed" / "run.log"
    ).read_text(encoding="utf-8")
    assert "progress=" in controlled_log
    assert "status=completed" in controlled_log

    report_run = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "make_report.py"),
            "--run",
            str(output),
            "--skip-animation",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert report_run.returncode == 0, report_run.stderr
    assert [path.stat().st_mtime_ns for path in result_files] == mtimes
    assert (output / "metrics.json").is_file()
    assert (output / "gap_error_spacetime.png").is_file()
