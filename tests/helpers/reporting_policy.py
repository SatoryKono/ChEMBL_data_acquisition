"""Utilities enforcing the reporting policy for the pytest suite."""

from __future__ import annotations

import platform
import subprocess
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

DEFAULT_REPO_SLUG = "SatoryKono/ChEMBL_data_acquisition"
SUCCESS_RATE_THRESHOLD = 0.95

PIPELINE_SCENARIOS: Mapping[str, str] = {
    "csv_loading": "Загрузка входных CSV и валидация схемы/типов/обязательных колонок",
    "normalization": "Нормализация и предобработка (включая кодировки, разделители)",
    "enrichment": "Обогащение данными из словарей/справочников (в т.ч. отсутствие соответствий)",
    "transformation_rules": "Правила трансформации и расчета флагов/категорий",
    "missing_data": "Обработка пропусков, дублей, конфликтных значений",
    "logging": "Детализация логирования и уровни WARN/ERROR при аномалиях",
    "assembly": "Итоговая сборка результатов: структура, сортировка, инварианты",
    "export": "Постобработка и экспорт: корректность форматов/имён/путей",
    "degradation": "Деградационные кейсы: частичные данные, пустые файлы, неверный заголовок",
    "idempotence": "Идемпотентность: повторный запуск на тех же входах даёт тот же результат",
}


SUMMARY_FIELDS = ("passed", "failed", "skipped", "xfailed", "xpassed", "error")
SUCCESS_STATUSES = frozenset({"passed", "xpassed"})
STATUS_PRIORITY = {
    "error": 50,
    "failed": 40,
    "xpassed": 30,
    "passed": 20,
    "xfailed": 10,
    "skipped": 0,
}


@dataclass
class ScenarioRecord:
    """Track which tests cover a given pipeline scenario."""

    name: str
    description: str
    tests: set[str] = field(default_factory=set)

    def register(self, nodeid: str) -> None:
        self.tests.add(nodeid)

    def covered(self, outcomes: Mapping[str, str]) -> bool:
        return any(outcomes.get(nodeid) in SUCCESS_STATUSES for nodeid in self.tests)

    def to_json(self, covered: bool) -> dict[str, Any]:
        return {
            "description": self.description,
            "tests": sorted(self.tests),
            "covered": covered,
        }


class ReportingPolicyState:
    """Collect per-test results and enforce policy constraints."""

    def __init__(
        self,
        *,
        scenarios: Mapping[str, str] | None = None,
        success_threshold: float = SUCCESS_RATE_THRESHOLD,
    ) -> None:
        scenario_map = scenarios or PIPELINE_SCENARIOS
        self.scenarios: dict[str, ScenarioRecord] = {
            name: ScenarioRecord(name=name, description=description)
            for name, description in scenario_map.items()
        }
        self.success_threshold = success_threshold
        self.node_status: dict[str, str] = {}
        self.summary: dict[str, int | float] = {field: 0 for field in SUMMARY_FIELDS}
        self.summary["total"] = 0
        self.success_rate: float = 1.0
        self.duration_sec: float = 0.0
        self.failure_reasons: list[str] = []
        self.missing_scenarios: list[str] = []

    # registration -------------------------------------------------
    def register_item(self, nodeid: str, scenario_names: Iterable[str]) -> None:
        for raw_name in scenario_names:
            name = str(raw_name)
            if name not in self.scenarios:
                known = ", ".join(sorted(self.scenarios))
                raise ValueError(
                    f"Unknown pipeline scenario '{name}'. Known scenarios: {known}"
                )
            self.scenarios[name].register(nodeid)

    def record_status(self, nodeid: str, status: str) -> None:
        priority = STATUS_PRIORITY.get(status, -1)
        previous = self.node_status.get(nodeid)
        if previous is None:
            self.node_status[nodeid] = status
            return
        if priority >= STATUS_PRIORITY.get(previous, -1):
            self.node_status[nodeid] = status

    # synthesis ----------------------------------------------------
    def finalise(self, duration_sec: float) -> None:
        self.duration_sec = max(float(duration_sec), 0.0)
        for summary_field in SUMMARY_FIELDS:
            self.summary[summary_field] = 0
        self.summary["total"] = len(self.node_status)

        for status in self.node_status.values():
            if status in SUMMARY_FIELDS:
                self.summary[status] += 1

        executed = self.summary["total"] - self.summary["skipped"]
        successes = self.summary["passed"] + self.summary["xfailed"]
        if executed > 0:
            self.success_rate = max(0.0, min(successes / executed, 1.0))
        else:
            self.success_rate = 1.0

        missing: list[str] = []
        for name, record in self.scenarios.items():
            if not record.tests or not record.covered(self.node_status):
                missing.append(name)
        self.missing_scenarios = sorted(missing)

        reasons: list[str] = []
        if self.success_rate + 1e-9 < self.success_threshold:
            pct = self.success_rate * 100.0
            threshold_pct = self.success_threshold * 100.0
            reasons.append(
                f"Success rate {pct:.2f}% is below the {threshold_pct:.2f}% threshold"
            )
        if self.missing_scenarios:
            reasons.append(
                "Missing pipeline scenarios: "
                + ", ".join(self.missing_scenarios)
            )
        self.failure_reasons = reasons

    # reporting ----------------------------------------------------
    def build_summary(self) -> dict[str, Any]:
        summary = {
            "total": int(self.summary["total"]),
            "success_rate": round(self.success_rate, 4),
        }
        summary.update({field: int(self.summary[field]) for field in SUMMARY_FIELDS})
        return summary

    def build_meta(
        self,
        *,
        repo: str | None = None,
        commit: str | None = None,
        branch: str | None = None,
        python_version: str | None = None,
        pytest_version: str | None = None,
    ) -> dict[str, Any]:
        repo_value = repo or DEFAULT_REPO_SLUG
        meta = {
            "repo": repo_value,
            "commit": commit or "",
            "branch": branch or "",
            "ts_utc": datetime.now(UTC).isoformat(),
            "duration_sec": round(self.duration_sec, 3),
            "python": python_version or platform.python_version(),
            "pytest": pytest_version or pytest.__version__,
            "success_threshold": self.success_threshold,
            "success_rate": round(self.success_rate, 4),
        }
        meta["pipeline_scenarios"] = {
            name: record.to_json(name not in self.missing_scenarios)
            for name, record in self.scenarios.items()
        }
        return meta

    def update_json_report(
        self,
        report: dict[str, Any],
        *,
        repo: str | None = None,
        commit: str | None = None,
        branch: str | None = None,
        python_version: str | None = None,
        pytest_version: str | None = None,
    ) -> None:
        report["summary"] = self.build_summary()
        report["meta"] = self.build_meta(
            repo=repo,
            commit=commit,
            branch=branch,
            python_version=python_version,
            pytest_version=pytest_version,
        )
        report.setdefault("tests", [])

    # helpers ------------------------------------------------------
    def scenario_description(self, name: str) -> str:
        return self.scenarios[name].description

    def gate_failed(self) -> bool:
        return bool(self.failure_reasons)


def git_output(*args: str, cwd: Path | None = None) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):  # pragma: no cover
        return ""
    return completed.stdout.strip()


class ReportingPolicyPlugin:
    """Pytest plugin wiring :class:`ReportingPolicyState` into the session."""

    def __init__(self, config: pytest.Config) -> None:
        self.config = config
        self.state = ReportingPolicyState()
        self._start_time: float = 0.0

    # pytest hook implementations ---------------------------------
    def pytest_sessionstart(self, session: pytest.Session) -> None:  # noqa: D401
        self._start_time = time.perf_counter()

    def pytest_collection_modifyitems(
        self,
        session: pytest.Session,
        config: pytest.Config,
        items: list[pytest.Item],
    ) -> None:  # noqa: D401
        for item in items:
            marker = item.get_closest_marker("pipeline_scenario")
            if marker is None:
                continue
            if not marker.args:
                raise pytest.UsageError(
                    f"pipeline_scenario marker on {item.nodeid} requires at least one scenario name"
                )
            self.state.register_item(item.nodeid, marker.args)

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:  # noqa: D401
        status = _determine_status(report)
        if status is not None:
            self.state.record_status(report.nodeid, status)

    def pytest_sessionfinish(
        self, session: pytest.Session, exitstatus: int | pytest.ExitCode
    ) -> None:  # noqa: D401
        duration = time.perf_counter() - self._start_time
        self.state.finalise(duration)
        if self.state.gate_failed():
            if exitstatus == 0 or exitstatus == pytest.ExitCode.OK:
                session.exitstatus = pytest.ExitCode.TESTS_FAILED

    def pytest_terminal_summary(  # noqa: D401
        self, terminalreporter: pytest.TerminalReporter
    ) -> None:
        terminalreporter.section("Reporting policy", sep="=")
        threshold_pct = self.state.success_threshold * 100.0
        terminalreporter.write_line(
            f"Success rate: {self.state.success_rate * 100.0:.2f}% (threshold {threshold_pct:.2f}%)"
        )
        if self.state.missing_scenarios:
            terminalreporter.write_line("Missing pipeline scenarios:")
            for name in self.state.missing_scenarios:
                description = self.state.scenario_description(name)
                terminalreporter.write_line(f"- {name}: {description}")
        else:
            terminalreporter.write_line("All pipeline scenarios covered.")
        if self.state.failure_reasons:
            terminalreporter.write_line("Policy gate failed:")
            for reason in self.state.failure_reasons:
                terminalreporter.write_line(f"- {reason}")
        else:
            terminalreporter.write_line("Policy gate satisfied.")

    @pytest.hookimpl(optionalhook=True)
    def pytest_json_modifyreport(self, json_report: dict[str, Any]) -> None:  # noqa: D401
        repo = git_output("config", "--get", "remote.origin.url") or DEFAULT_REPO_SLUG
        commit = git_output("rev-parse", "HEAD")
        branch = git_output("rev-parse", "--abbrev-ref", "HEAD")
        self.state.update_json_report(
            json_report,
            repo=repo,
            commit=commit,
            branch=branch,
            python_version=platform.python_version(),
            pytest_version=pytest.__version__,
        )


def _determine_status(report: pytest.TestReport) -> str | None:
    if report.skipped:
        return "xfailed" if getattr(report, "wasxfail", False) else "skipped"
    if report.failed:
        return "failed" if report.when == "call" else "error"
    if report.when == "call":
        return "xpassed" if getattr(report, "wasxfail", False) else "passed"
    return None

