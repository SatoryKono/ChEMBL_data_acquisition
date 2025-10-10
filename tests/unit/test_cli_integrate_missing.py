from pathlib import Path

import pandas as pd
import pytest


def _cli():
    from library.pipelines.testitem import cli

    return cli


@pytest.mark.unit
def test_integrate_missing_identifiers__appends_placeholders_in_order() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL2"],
            "value": [10],
        }
    )
    missing = ["CHEMBL3"]
    requested = ["CHEMBL2", "CHEMBL1", "CHEMBL3"]

    result = _cli().integrate_missing_identifiers(
        df,
        missing_ids=missing,
        requested_ids=requested,
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL2", "CHEMBL3"]
    assert pd.isna(result.loc[1, "value"])


@pytest.mark.unit
def test_integrate_missing_identifiers__restores_requested_order() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL3", "CHEMBL1"],
        }
    )
    missing: list[str] = []
    requested = ["CHEMBL1", "CHEMBL2", "CHEMBL3"]

    result = _cli().integrate_missing_identifiers(
        df,
        missing_ids=missing,
        requested_ids=requested,
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL3"]


@pytest.mark.unit
def test_integrate_missing_identifiers__ignores_blank_requested_ids() -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    result = _cli().integrate_missing_identifiers(
        df,
        missing_ids=["CHEMBL2"],
        requested_ids=["", None, "CHEMBL1", "CHEMBL2"],
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]


@pytest.mark.unit
def test_integrate_missing_identifiers__missing_ids_appended_at_end() -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL3"]})
    result = _cli().integrate_missing_identifiers(
        df,
        missing_ids=["CHEMBL1"],
        requested_ids=["CHEMBL3"],
    )

    assert list(result["molecule_chembl_id"]) == ["CHEMBL3", "CHEMBL1"]


@pytest.mark.unit
def test_finalize_output__optional_column_present_in_later_chunk(
    tmp_path: Path,
    sample_input_csv: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    second_chunk = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL2"],
            "pref_name": ["Example"],
        }
    )
    output_path = tmp_path / "optional.csv"

    warnings: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        warnings.append((event, fields))

    cli = _cli()

    monkeypatch.setattr(cli.logger, "warning", capture_warning)

    from library.config import Config
    from library.pipelines.testitem.catalog import ParentLookupStats
    from library.resources import dictionaries as dictionary_resources

    def _fake_resource_path(name: str, base_dir: Path | None = None) -> Path:
        root = tmp_path / "resources"
        path = root / f"{name}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.write_text("{}", encoding="utf-8")
        return path

    def _fake_get_resource(name: str, base_dir: Path | None = None):
        path = _fake_resource_path(name, base_dir)
        return dictionary_resources.DictionaryResource(
            name=name,
            relative_path=Path(f"{name}.json"),
            path=path,
            version="test",
            sha256="dummy",
            generator=Path("generator"),
        )

    monkeypatch.setattr(dictionary_resources, "get_resource_path", _fake_resource_path)
    monkeypatch.setattr(dictionary_resources, "get_resource", _fake_get_resource)

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    cfg.sources.pubchem.user_agent = "pubchem-contact@example.org"
    cfg.system.doc_quality.enable = False

    stats = ParentLookupStats(
        source="test",
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
        failed_ids=(),
        hierarchy_attached=0,
        fallback_attached=0,
        no_parent=0,
    )

    exit_code = cli.finalize_output(
        [first_chunk, second_chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=lambda: stats,
        input_csv=sample_input_csv,
    )

    assert exit_code == 0

    optional_warnings = [
        fields for event, fields in warnings if event == "optional_columns_missing"
    ]
    for payload in optional_warnings:
        columns = payload.get("columns", [])
        assert "pref_name" not in columns

    final = pd.read_csv(output_path)
    assert "pref_name" in final.columns
    assert final.loc[1, "pref_name"] == "Example"
