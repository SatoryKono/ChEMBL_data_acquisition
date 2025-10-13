"""Ensure the reporting policy scenarios and test markers stay in sync."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.helpers.reporting_policy import PIPELINE_SCENARIOS

TESTS_ROOT = Path(__file__).resolve().parents[1]


def _collect_pipeline_scenarios_from_tests() -> set[str]:
    scenario_names: set[str] = set()
    for path in TESTS_ROOT.rglob("test_*.py"):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            for decorator in node.decorator_list:
                names = _extract_marker_scenarios(decorator, path)
                scenario_names.update(names)
    return scenario_names


def _extract_marker_scenarios(node: ast.AST, path: Path) -> set[str]:
    if not isinstance(node, ast.Call):
        return set()
    func = node.func
    if not (
        isinstance(func, ast.Attribute)
        and func.attr == "pipeline_scenario"
        and isinstance(func.value, ast.Attribute)
        and func.value.attr == "mark"
        and isinstance(func.value.value, ast.Name)
        and func.value.value.id == "pytest"
    ):
        return set()

    scenarios: set[str] = set()
    for arg in node.args:
        if not isinstance(arg, ast.Constant) or not isinstance(arg.value, str):
            raise AssertionError(
                "pipeline_scenario marker must use literal string arguments; "
                f"found unsupported argument in {path}:{getattr(arg, 'lineno', '?')}"
            )
        scenarios.add(arg.value)
    return scenarios


@pytest.mark.unit
def test_reporting_policy_pipeline_scenarios__in_sync() -> None:
    declared = set(PIPELINE_SCENARIOS)
    referenced = _collect_pipeline_scenarios_from_tests()

    missing_in_policy = referenced - declared
    missing_in_tests = declared - referenced

    assert not missing_in_policy, (
        "Found pipeline_scenario markers without policy declaration: "
        + ", ".join(sorted(missing_in_policy))
    )
    assert not missing_in_tests, (
        "Reporting policy scenarios lack covering tests: "
        + ", ".join(sorted(missing_in_tests))
    )
