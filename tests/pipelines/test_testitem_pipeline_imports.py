"""Tests covering import helpers for :mod:`library.testitem_pipeline`."""

from __future__ import annotations

import importlib

import pytest


def test_legacy_namespace_no_longer_exposes_loader() -> None:
    """The historic namespace proxies to the modern pipeline without loaders."""

    module = importlib.import_module("library.testitem_pipeline")

    with pytest.raises(AttributeError):
        getattr(module, "_load_local_module")
