"""Unit tests for :mod:`library.postprocessing.helpers`."""

from pathlib import Path

import pytest

from library.postprocessing.helpers import normalise_export_basename


@pytest.mark.unit
def test_normalise_export_basename__preserves_regular_names() -> None:
    source = Path("output.targets_20250101.csv")
    assert normalise_export_basename(source) == "output.targets_20250101.csv"


@pytest.mark.unit
def test_normalise_export_basename__drops_atomic_prefix_and_suffix() -> None:
    source = Path(".output.targets_20250101.csv.tmp")
    assert normalise_export_basename(source) == "output.targets_20250101.csv"


@pytest.mark.unit
def test_normalise_export_basename__restores_csv_suffix_order() -> None:
    source = Path("output.targets_20250101.csv_normalized.tmp")
    assert normalise_export_basename(source) == "output.targets_20250101_normalized.csv"
