from __future__ import annotations

from scripts import get_target_data as gtd


def test_organism_csv_option_removed(tmp_path, monkeypatch) -> None:
    """The ``all`` command no longer exposes an ``--organism-csv`` option."""
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "sources:\n"
        "  chembl:\n"
        "    pipelines:\n"
        "      target:\n"
        "        all:\n"
        "          organism_csv: dictionary/organism.csv\n"
        "    api:\n"
        "      timeout_read: 30\n"
        "local:\n"
        "  io:\n"
        "    csv_sep: ','\n"
        "    csv_encoding: utf8\n"
        "  resources:\n"
        "    dictionary_dir: dictionary\n"
        "    iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "    iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "    targets_type_csv: dictionary/targets_type.csv\n"
        "system:\n"
        "  log:\n"
        "    level: INFO\n"
    )
    captured: dict[str, object] = {}

    def fake_run_all(cfg, args) -> int:  # type: ignore[unused-argument]
        captured["has_organism_csv"] = hasattr(args, "organism_csv")
        return 0

    monkeypatch.setattr(gtd, "run_all", fake_run_all)
    rc = gtd.main(["all", "--config", str(cfg_path)])
    assert rc == 0
    assert not captured["has_organism_csv"]
