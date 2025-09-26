import io
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from pydantic import BaseModel, Field, ValidationError

from library.cli import LoggerConfig, configure_logger
from library.config import (
    Config,
    ConfigError,
    build_alias_map,
    ensure_dirs,
    load_config,
)
from scripts import get_target_data as target_cli


@pytest.fixture(autouse=True)
def _user_agent_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Provide a default user agent for configuration loading."""

    monkeypatch.setenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT",
        "test-agent/1.0 (mailto:test@example.org)",
    )


def test_config_to_dict(tmp_path: Path) -> None:
    """``Config.to_dict`` should mirror :meth:`model_dump`."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path)
    assert cfg.to_dict() == cfg.model_dump()


def test_load_minimal_config(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path)
    assert cfg.api.rps == 5
    assert cfg.openalex.rps == 4


def test_missing_config_raises(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="configuration file not found"):
        load_config(tmp_path / "missing.yaml")


def test_schema_file_rejected(tmp_path: Path) -> None:
    """Passing the JSON schema instead of a config file should fail."""

    path = tmp_path / "schema.json"
    path.write_text('{"$defs": {}}')
    with pytest.raises(ConfigError, match="configuration schema"):
        load_config(path)


def test_env_overrides_yaml(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n"
        "  chembl:\n"
        "    api:\n"
        "      rps: 1\n"
        "  openalex:\n"
        "    rps: 1\n"
    )
    monkeypatch.setenv("CHEMBL_DA__SOURCES__CHEMBL__API__RPS", "3")
    monkeypatch.setenv("CHEMBL_DA__SOURCES__OPENALEX__RPS", "6")
    cfg = load_config(path)
    assert cfg.api.rps == 3
    assert cfg.openalex.rps == 6


def test_alias_env_overrides(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure shorthand environment variables map to the correct fields."""

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n"
        "  chembl:\n"
        "    api:\n"
        "      chembl_base: https://www.ebi.ac.uk/chembl/api/data\n"
        "      timeout_connect: 1\n"
        "      timeout_read: 1\n"
    )

    monkeypatch.setenv("CHEMBL_DA_BASE", "https://example.org")
    monkeypatch.setenv("CHEMBL_DA_TIMEOUT_CONNECT", "10")
    monkeypatch.setenv("CHEMBL_DA_TIMEOUT_READ", "20")
    cfg = load_config(path)

    assert cfg.api.chembl_base == "https://example.org"
    assert cfg.api.timeout_connect == 10
    assert cfg.api.timeout_read == 20
    assert cfg.api.rps == 5


def test_retry_and_log_aliases(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """New environment variable aliases should override retry and log defaults."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")

    monkeypatch.setenv("CHEMBL_DA_RETRY_MAX_ATTEMPTS", "10")
    monkeypatch.setenv("CHEMBL_DA_RETRY_BACKOFF_FACTOR", "2.0")
    monkeypatch.setenv("CHEMBL_DA_LOG_FORMAT", "%(levelname)s")
    monkeypatch.setenv("CHEMBL_DA_LOG_DATEFMT", "%d/%m/%Y")

    cfg = load_config(path)

    assert cfg.retry.max_attempts == 10
    assert cfg.retry.backoff_factor == 2.0
    assert cfg.log.format == "%(levelname)s"
    assert cfg.log.datefmt == "%d/%m/%Y"


def test_cache_dir_alias(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The ``CHEMBL_DA_CACHE_DIR`` alias should override the cache path."""

    cache = tmp_path / "cache"
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_CACHE_DIR", str(cache))

    cfg = load_config(path)

    assert cfg.io.cache_dir == cache


def test_molecule_catalog_cache_alias(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cache = tmp_path / "catalog.json"
    monkeypatch.setenv("CHEMBL_DA_MOLECULE_CATALOG_CACHE", str(cache))

    cfg = load_config(path)

    assert cfg.molecule_catalog.cache_path == cache


def test_chembl_cache_maxsize_from_yaml(tmp_path: Path) -> None:
    """Custom ``chembl.cache_maxsize`` should override the default."""

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n  chembl:\n    cache:\n      cache_maxsize: 42\n"
    )

    cfg = load_config(path)

    assert cfg.chembl.cache_maxsize == 42


def test_mailto_aliases_override_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """New mailto aliases should override the default email addresses."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")

    monkeypatch.setenv("CHEMBL_DA_OPENALEX_MAILTO", "openalex@example.org")
    monkeypatch.setenv("CHEMBL_DA_CROSSREF_MAILTO", "crossref@example.org")

    cfg = load_config(path)

    assert cfg.openalex.mailto == "openalex@example.org"
    assert cfg.crossref.mailto == "crossref@example.org"


def test_base_aliases_override_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """New base aliases should override the default base URLs."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")

    monkeypatch.setenv("CHEMBL_DA_OPENALEX_BASE", "https://example.org/openalex")
    monkeypatch.setenv("CHEMBL_DA_CROSSREF_BASE", "https://example.org/crossref")
    monkeypatch.setenv("CHEMBL_DA_UNIPROT_BASE", "https://example.org/uniprot")
    monkeypatch.setenv("CHEMBL_DA_IUPHAR_BASE", "https://example.org/iuphar")
    monkeypatch.setenv("CHEMBL_DA_PUBCHEM_BASE", "https://example.org/pubchem")

    cfg = load_config(path)

    assert cfg.openalex.base == "https://example.org/openalex"
    assert cfg.crossref.base == "https://example.org/crossref"
    assert cfg.uniprot.base == "https://example.org/uniprot"
    assert cfg.iuphar.base == "https://example.org/iuphar"
    assert cfg.pubchem.base == "https://example.org/pubchem"


@pytest.mark.parametrize("section", ["pubmed", "semantic_scholar"])
def test_pubmed_semantic_bad_base(section: str, tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(f"sources:\n  {section}:\n    base: https://\n")
    with pytest.raises(ValueError, match=f"{section}.base"):
        load_config(path)


def test_cli_overrides_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("sources:\n  chembl:\n    api:\n      rps: 1\n")

    monkeypatch.setenv("CHEMBL_DA__SOURCES__CHEMBL__API__RPS", "2")
    cfg = load_config(path, cli_overrides={"sources.chembl.api.rps": 4})
    assert cfg.api.rps == 4


def test_cli_path_override(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    out = tmp_path / "out"
    cfg = load_config(path, cli_overrides={"local.io.output_dir": str(out)})
    assert cfg.io.output_dir == out


def test_doc_type_cli_override(tmp_path: Path) -> None:
    """CLI overrides should update document type settings."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(
        path,
        cli_overrides={
            "system.doc_type.weights.pubmed": 8,
            "system.doc_type.thresholds.review": 2,
        },
    )
    assert cfg.doc_type.weights["pubmed"] == 8
    assert cfg.doc_type.thresholds["review"] == 2


def test_type_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("sources:\n  chembl:\n    api:\n      rps: fast\n")
    with pytest.raises(ValidationError, match="api.rps"):
        load_config(path)


def test_schema_negative_value(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("sources:\n  openalex:\n    rps: -1\n")
    with pytest.raises(ValidationError):
        load_config(path)


def test_schema_list_item_type(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "system:\n  retry:\n    status_forcelist: [429, '500']\n"
    )
    with pytest.raises(ValidationError):
        load_config(path)


def test_missing_dirs_raise(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(tmp_path / "out"))
    monkeypatch.setenv("CHEMBL_DA__LOCAL__IO__CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("CHEMBL_DA__LOCAL__IO__EXIST_OK", "false")
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    with pytest.raises(FileNotFoundError):
        load_config(path)


def test_invalid_bool_env_var(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure unknown boolean strings raise :class:`ValueError`."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA__LOCAL__IO__EXIST_OK", "maybe")
    with pytest.raises(ValidationError, match="Invalid boolean value"):
        load_config(path)


def test_ensure_dirs_creates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "out"
    cache = tmp_path / "cache"
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(out))
    monkeypatch.setenv("CHEMBL_DA__IO__CACHE_DIR", str(cache))
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path)
    assert not out.exists() and not cache.exists()
    ensure_dirs(cfg)
    assert out.is_dir() and cache.is_dir()



def test_unknown_key_warning_non_strict(tmp_path: Path) -> None:

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "unknown: 1\n"
        "sources:\n"
        "  chembl:\n"
        "    api:\n"
        "      rps: 1\n"
    )
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    load_config(path, strict=False)
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    msg = record.get("msg", "") or record.get("event", "")
    assert "Unknown configuration key" in msg


def test_unknown_key_error(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("unknown: 1\n")
    with pytest.raises(ValueError, match="Unknown configuration key"):
        load_config(path)


def test_config_type_coercion() -> None:
    """The default configuration should load without type errors in strict mode."""

    path = Path(__file__).resolve().parents[1] / "config.yaml"
    cfg = load_config(path, strict=True)
    assert isinstance(cfg.pubchem.delay, float)


def test_default_resource_paths_exist() -> None:
    """Default resource paths should exist on disk."""

    cfg_path = Path(__file__).resolve().parents[1] / "config.yaml"
    cfg = load_config(cfg_path)
    resources = cfg.resources
    project_root = cfg_path.parent
    for field in (
        "dictionary_dir",
        "iuphar_target_csv",
        "iuphar_family_csv",
        "uniprot_data_dir",
    ):
        resource_path = getattr(resources, field)
        full_path = (
            resource_path
            if resource_path.is_absolute()
            else project_root / resource_path
        )
        assert full_path.exists(), f"Missing default resource: {full_path}"
def test_yaml_error_includes_path(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("sources:\n  chembl:\n    api: [\n")  # malformed YAML
    with pytest.raises(ConfigError) as excinfo:
        load_config(path)
    msg = str(excinfo.value)
    assert str(path) in msg
    assert "while parsing" in msg


def test_user_agent_must_include_contact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n"
        "  chembl:\n"
        "    api:\n"
        "      user_agent: chembl-da/0.1\n"
        "  openalex:\n"
        "    mailto: info@example.org\n"
        "  crossref:\n"
        "    mailto: info@example.org\n"
    )
    # Ensure the invalid YAML value is used rather than the environment
    # override provided by the autouse fixture.
    monkeypatch.delenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT", raising=False
    )
    with pytest.raises(ValidationError, match="user_agent"):
        load_config(path)


def test_user_agent_default_and_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Default user agent applies and can be overridden."""

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n"
        "  openalex:\n"
        "    mailto: info@example.org\n"
        "  crossref:\n"
        "    mailto: info@example.org\n"
    )
    monkeypatch.delenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT", raising=False
    )
    cfg = load_config(path)
    assert cfg.api.user_agent == "chembl-da/0.1 (mailto:contact@example.org)"

    monkeypatch.setenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT",
        "cli-agent/1.0 (mailto:test@example.org)",
    )
    cfg = load_config(path)
    assert cfg.api.user_agent == "cli-agent/1.0 (mailto:test@example.org)"

    monkeypatch.delenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT", raising=False
    )
    cfg = load_config(
        path,
        cli_overrides={
            "sources.chembl.api.user_agent": "override/1 (mailto:me@example.org)"
        },
    )
    assert cfg.api.user_agent == "override/1 (mailto:me@example.org)"


def test_openalex_mailto_required(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n"
        "  openalex:\n"
        "    mailto: ''\n"
        "  crossref:\n"
        "    mailto: info@example.org\n"
    )
    with pytest.raises(ValidationError, match="openalex.mailto"):
        load_config(path)


def test_crossref_mailto_format(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "sources:\n  crossref:\n    mailto: not-an-email\n"
    )
    with pytest.raises(ValidationError, match="crossref.mailto"):
        load_config(path)


@pytest.mark.parametrize(
    ("snippet", "match"),
    [
        (
            "sources:\n  chembl:\n    api:\n      chembl_base: https://\n",
            "sources.chembl.api.chembl_base",
        ),
        ("sources:\n  openalex:\n    base: https://\n", "sources.openalex.base"),
        ("sources:\n  crossref:\n    base: https://\n", "sources.crossref.base"),
        (
            "sources:\n  uniprot:\n    api:\n      base: https://\n",
            "sources.uniprot.api.base",
        ),
        ("sources:\n  iuphar:\n    base: https://\n", "sources.iuphar.base"),
        ("sources:\n  pubchem:\n    base: https://\n", "sources.pubchem.base"),
    ],
)
def test_invalid_urls_raise(tmp_path: Path, snippet: str, match: str) -> None:
    """Invalid base URLs should trigger :class:`ValueError`."""

    path = tmp_path / "cfg.yaml"
    path.write_text(snippet)
    with pytest.raises(ValidationError, match=match):
        load_config(path)


def test_unknown_env_var_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA__FOO__BAR", "1")
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    load_config(path)
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    msg = record.get("msg", "") or record.get("event", "")
    assert "Environment variable CHEMBL_DA__FOO__BAR ignored" in msg


def test_new_field_auto_alias() -> None:
    """Adding a field to the config should expose an alias automatically."""

    class Extra(BaseModel):
        foo: int = 1

    class ExtendedConfig(Config):
        extra: Extra = Field(default_factory=Extra)

    aliases = build_alias_map(ExtendedConfig)
    assert aliases["CHEMBL_DA_EXTRA_FOO"] == ["extra", "foo"]


def test_log_level_valid(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("system:\n  log:\n    level: warn\n")
    cfg = load_config(path)
    assert cfg.log.level == "warn"


def test_log_level_invalid(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("system:\n  log:\n    level: verbose\n")
    with pytest.raises(ValidationError) as exc:
        load_config(path)
    valid = ", ".join(sorted(logging.getLevelNamesMapping()))
    assert f"log.level must be one of {valid}, got 'verbose'" in str(exc.value)


def test_log_level_valid_no_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fallback mapping should validate known log levels."""

    path = tmp_path / "cfg.yaml"
    path.write_text("system:\n  log:\n    level: warn\n")
    monkeypatch.delattr(logging, "getLevelNamesMapping", raising=False)
    cfg = load_config(path)
    assert cfg.log.level == "warn"


def test_target_chembl_defaults_match_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Default target configuration should work with the CLI helpers."""

    cfg = load_config(Path("config.yaml"))
    assert cfg.target.chembl.column == "target_chembl_id"

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n")
    output_csv = tmp_path / "out.csv"

    class DummyClient:
        def __init__(self, *_: object, **__: object) -> None:
            pass

        def __enter__(self) -> "DummyClient":
            return self

        def __exit__(
            self,
            exc_type: object,
            exc: object,
            traceback: object,
        ) -> bool:
            return False

    def fake_get_targets(
        ids: object,
        *,
        cfg: object,
        client: object,
        mapping_cfg: object,

        chunk_size: object | None = None,

        timeout: object,
    ) -> pd.DataFrame:
        return pd.DataFrame({"target_chembl_id": list(ids)})

    monkeypatch.setattr(target_cli, "ChemblClient", DummyClient)
    monkeypatch.setattr(target_cli.cl, "get_targets", fake_get_targets)

    args = SimpleNamespace(input_csv=input_csv, output_csv=output_csv, limit=None)

    exit_code = target_cli.run_chembl(cfg, args)

    assert exit_code == 0
    assert output_csv.exists()
    header = output_csv.read_text(encoding=cfg.io.csv_encoding).splitlines()[0]
    expected_header = cfg.io.csv_sep.join(
        ["pipeline_version", "target_chembl_id", "timestamp_utc"]
    )
    assert header == expected_header


def test_log_level_invalid_no_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fallback mapping should reject unknown levels."""

    path = tmp_path / "cfg.yaml"
    path.write_text("system:\n  log:\n    level: verbose\n")
    monkeypatch.delattr(logging, "getLevelNamesMapping", raising=False)
    with pytest.raises(ValidationError) as exc:
        load_config(path)
    level_names = {name.upper(): level for name, level in logging._nameToLevel.items()}
    valid = ", ".join(sorted(level_names))
    assert f"log.level must be one of {valid}, got 'verbose'" in str(exc.value)
