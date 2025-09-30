"""Packaging regression tests for bundled resources."""

from __future__ import annotations

import importlib.resources as importlib_resources
from pathlib import Path

from library.config import Config, ResourcesCfg


def _assert_path_exists(path: Path) -> None:
    assert path.exists(), f"Expected resource path to exist: {path}"


def test_dictionary_package_is_importable() -> None:
    """The ``dictionary`` package should be accessible via importlib.resources."""

    package = importlib_resources.files("dictionary")
    with importlib_resources.as_file(package) as path:
        assert path.is_dir()
        _assert_path_exists(path)


def test_resources_cfg_default_paths_exist() -> None:
    """Default :class:`ResourcesCfg` paths should resolve to real files."""

    cfg = ResourcesCfg()
    for field in (
        "dictionary_dir",
        "iuphar_target_csv",
        "iuphar_family_csv",
        "uniprot_data_dir",
        "targets_type_csv",
    ):
        value = getattr(cfg, field)
        _assert_path_exists(value)

    cfg_with_relative = ResourcesCfg(dictionary_dir=Path("dictionary"))
    _assert_path_exists(cfg_with_relative.dictionary_dir)

    config = Config()
    enrichment_sources = config.testitem_molecule_enrichment.sources
    _assert_path_exists(enrichment_sources.molecule_catalog_path)
    _assert_path_exists(enrichment_sources.molecule_hierarchy_path)
    target_cfg = config.sources.chembl.pipelines.target
    _assert_path_exists(target_cfg.uniprot.data_dir)
    _assert_path_exists(target_cfg.iuphar.target_csv)
    _assert_path_exists(target_cfg.iuphar.family_csv)
    _assert_path_exists(target_cfg.all.data_dir)
