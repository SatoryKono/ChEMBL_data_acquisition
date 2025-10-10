"""Unit tests for the success rate extraction helper."""

from __future__ import annotations

import pytest

from scripts.check_test_success_rate import _extract_success_rate


def test_extract_success_rate__fallback_counts_skipped_and_xfailed() -> None:
    summary = {"total": 10, "passed": 7, "skipped": 2, "xfailed": 1}

    rate = _extract_success_rate(summary)

    assert rate == pytest.approx((7 + 1) / 8)


def test_extract_success_rate__fallback_denominator_floor() -> None:
    summary = {"total": 4, "passed": 4, "skipped": 4, "xfailed": 0}

    rate = _extract_success_rate(summary)

    assert rate == 1.0


def test_extract_success_rate__prefers_existing_success_rate() -> None:
    summary = {
        "success_rate": 0.42,
        "total": 10,
        "passed": 10,
        "skipped": 10,
        "xfailed": 0,
    }

    rate = _extract_success_rate(summary)

    assert rate == pytest.approx(0.42)
