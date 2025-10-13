from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import textwrap

import pytest


def _create_fake_pytest_module(base_dir: Path) -> Path:
    """Return a directory containing a minimal stand-in for :mod:`pytest`."""

    package_dir = base_dir / "pytest"
    package_dir.mkdir(parents=True, exist_ok=True)

    init_path = package_dir / "__init__.py"
    init_path.write_text(
        "\n".join(
            (
                "'''A lightweight pytest stub used in tests.'''",
                '__all__ = ["main"]',
                '__version__ = "fake-0.0"',
                "",
                "def main() -> int:",
                "    from . import __main__ as _main",
                "",
                "    return _main.main()",
                "",
            )
        ),
        encoding="utf-8",
    )

    main_path = package_dir / "__main__.py"
    main_path.write_text(
        "\n".join(
            (
                "from __future__ import annotations",
                "",
                "import json",
                "import sys",
                "from pathlib import Path",
                "",
                "_DEFAULT_TEST = {",
                '    "nodeid": "tests/unit/test_sample.py::test_placeholder",',
                '    "outcome": "passed",',
                '    "setup": {"duration": 0.001},',
                '    "call": {"duration": 0.002},',
                '    "teardown": {"duration": 0.001},',
                "}",
                "",
                "def _parse_args(argv: list[str]) -> tuple[Path | None, Path | None, Path | None]:",
                "    json_report: Path | None = None",
                "    coverage_xml: Path | None = None",
                "    log_file: Path | None = None",
                "    i = 0",
                "    while i < len(argv):",
                "        arg = argv[i]",
                "        if arg == \"--json-report-file\" and i + 1 < len(argv):",
                "            json_report = Path(argv[i + 1])",
                "            i += 1",
                "        elif arg.startswith(\"--cov-report=\"):",
                "            value = arg.split(\"=\", 1)[1]",
                "            if value.startswith(\"xml:\"):",
                "                coverage_xml = Path(value.split(\":\", 1)[1])",
                "            elif value.startswith(\"html:\"):",
                "                html_dir = Path(value.split(\":\", 1)[1])",
                "                html_dir.mkdir(parents=True, exist_ok=True)",
                "                (html_dir / \"index.html\").write_text(\"<html>fake</html>\", encoding=\"utf-8\")",
                "        elif arg == \"--log-file\" and i + 1 < len(argv):",
                "            log_file = Path(argv[i + 1])",
                "            i += 1",
                "        i += 1",
                "    return json_report, coverage_xml, log_file",
                "",
                "def main() -> int:",
                "    json_report, coverage_xml, log_file = _parse_args(sys.argv[1:])",
                "",
                "    if log_file is not None:",
                "        log_file.parent.mkdir(parents=True, exist_ok=True)",
                "        log_file.write_text(\"fake pytest log\\n\", encoding=\"utf-8\")",
                "",
                "    if json_report is not None:",
                "        json_report.parent.mkdir(parents=True, exist_ok=True)",
                "        payload = {",
                '            "created": "2020-01-01T00:00:00",',
                '            "duration": 0.01,',
                '            "collectors": ["json-report"],',
                '            "tests": [_DEFAULT_TEST],',
                '            "summary": {"total": 1, "passed": 1, "failed": 0, "skipped": 0},',
                "        }",
                "        json_report.write_text(json.dumps(payload), encoding=\"utf-8\")",
                "",
                "    if coverage_xml is not None:",
                "        coverage_xml.parent.mkdir(parents=True, exist_ok=True)",
                "        coverage_xml.write_text('<coverage line-rate=\"1.0\"></coverage>', encoding=\"utf-8\")",
                "",
                "    return 0",
                "",
                'if __name__ == "__main__":',
                "    raise SystemExit(main())",
                "",
            )
        ),
        encoding="utf-8",
    )

    return base_dir


@pytest.mark.e2e
def test_run_tests_script__generates_reports(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    repo_copy = tmp_path / "repo"

    shutil.copytree(
        repo_root,
        repo_copy,
        ignore=shutil.ignore_patterns(
            "__pycache__",
            "*.pyc",
            "*.pyo",
            ".pytest_cache",
            ".mypy_cache",
            ".git",
            "output",
            "data",
            "reports",
        ),
    )

    fake_pytest_dir = tmp_path / "fake_pytest"
    _create_fake_pytest_module(fake_pytest_dir)

    outputs_dir = tmp_path / "artifacts"
    outputs_dir.mkdir()
    report_path = outputs_dir / "report.json"
    summary_path = outputs_dir / "summary.md"

    env = os.environ.copy()
    python_path_entries = [str(fake_pytest_dir), str(repo_copy)]
    existing_pythonpath = env.get("PYTHONPATH")
    if existing_pythonpath:
        python_path_entries.append(existing_pythonpath)
    env["PYTHONPATH"] = os.pathsep.join(python_path_entries)

    command = [
        sys.executable,
        "scripts/run_tests.py",
        "--pytest-timeout",
        "30",
        "--json",
        str(report_path),
        "--markdown",
        str(summary_path),
    ]

    completed = subprocess.run(
        command,
        cwd=repo_copy,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout

    assert report_path.is_file(), f"Structured JSON report missing at {report_path!s}"
    assert summary_path.is_file(), f"Markdown summary missing at {summary_path!s}"

    report_data = json.loads(report_path.read_text(encoding="utf-8"))

    assert report_data["summary"]["total"] == 1
    assert report_data["summary"]["passed"] == 1
    assert report_data["summary"]["success_rate"] == pytest.approx(1.0)
    assert report_data["meta"]["exit_code"] == 0
    assert any(test["status"] == "passed" for test in report_data["tests"])

    summary_text = summary_path.read_text(encoding="utf-8")
    assert "# Test Summary" in summary_text
    assert "Success rate: 100.00%" in summary_text
    assert "| total | passed | failed |" in summary_text

