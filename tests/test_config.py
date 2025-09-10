from pathlib import Path
import logging

import pytest

from jsonschema import ValidationError

from library.config import ConfigError, ensure_dirs, load_config


def test_load_minimal_config(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path)
    assert cfg.api.rps == 5
    assert cfg.openalex.rps == 4


def test_missing_config_raises(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="configuration file not found"):
        load_config(tmp_path / "missing.yaml")


def test_env_overrides_yaml(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 1\nopenalex:\n  rps: 1\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "3")
    monkeypatch.setenv("CHEMBL_DA__OPENALEX__RPS", "6")
    cfg = load_config(path)
    assert cfg.api.rps == 3
    assert cfg.openalex.rps == 6


def test_alias_env_overrides(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure shorthand environment variables map to the correct fields."""

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "api:\n"
        "  chembl_base: https://www.ebi.ac.uk/chembl/api/data\n"
        "  timeout_connect: 1\n"
        "  timeout_read: 1\n"
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

    cfg = load_config(path)

    assert cfg.retry.max_attempts == 10
    assert cfg.retry.backoff_factor == 2.0
    assert cfg.log.format == "%(levelname)s"


def test_cli_overrides_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 1\n")

    monkeypatch.setenv("CHEMBL_DA__API__RPS", "2")
    cfg = load_config(path, cli_overrides={"api.rps": 4})
    assert cfg.api.rps == 4


def test_cli_path_override(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    out = tmp_path / "out"
    cfg = load_config(path, cli_overrides={"io.output_dir": str(out)})
    assert cfg.io.output_dir == out


def test_type_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("api:\n  rps: fast\n")
    with pytest.raises(TypeError, match="api.rps must be int"):
        load_config(path)


def test_schema_negative_value(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("openalex:\n  rps: -1\n")
    with pytest.raises(ValidationError):
        load_config(path)


def test_schema_list_item_type(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("retry:\n  status_forcelist: [429, '500']\n")
    with pytest.raises(ValidationError):
        load_config(path)


def test_missing_dirs_raise(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(tmp_path / "out"))
    monkeypatch.setenv("CHEMBL_DA__IO__CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("CHEMBL_DA__IO__EXIST_OK", "false")
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    with pytest.raises(FileNotFoundError):
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


def test_unknown_key_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("unknown: 1\napi:\n  rps: 1\n")
    with caplog.at_level(logging.WARNING, logger="library.config"):
        load_config(path)
    assert "Unknown configuration key" in caplog.text


def test_unknown_key_error(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("unknown: 1\n")
    with pytest.raises(ValueError, match="Unknown configuration key"):
        load_config(path, strict=True)


def test_yaml_error_includes_path(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api: [\n")  # malformed YAML
    with pytest.raises(ConfigError) as excinfo:
        load_config(path)
    msg = str(excinfo.value)
    assert str(path) in msg
    assert "while parsing" in msg


def test_user_agent_must_include_contact(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "api:\n  user_agent: chembl-da/0.1\n"
        "openalex:\n  mailto: info@example.org\n"
        "crossref:\n  mailto: info@example.org\n"
    )
    with pytest.raises(ValueError, match="user_agent"):
        load_config(path)


def test_openalex_mailto_required(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text(
        "openalex:\n  mailto: ''\n" "crossref:\n  mailto: info@example.org\n"
    )
    with pytest.raises(ValueError, match="openalex.mailto"):
        load_config(path)


def test_crossref_mailto_format(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("crossref:\n  mailto: not-an-email\n")
    with pytest.raises(ValueError, match="crossref.mailto"):
        load_config(path)


@pytest.mark.parametrize(
    ("snippet", "match"),
    [
        ("api:\n  chembl_base: https://\n", "api.chembl_base"),
        ("openalex:\n  base: https://\n", "openalex.base"),
        ("crossref:\n  base: https://\n", "crossref.base"),
        ("uniprot:\n  base: https://\n", "uniprot.base"),
        ("iuphar:\n  base: https://\n", "iuphar.base"),
        ("pubchem:\n  base: https://\n", "pubchem.base"),
    ],
)
def test_invalid_urls_raise(tmp_path: Path, snippet: str, match: str) -> None:
    """Invalid base URLs should trigger :class:`ValueError`."""

    path = tmp_path / "cfg.yaml"
    path.write_text(snippet)
    with pytest.raises(ValueError, match=match):
        load_config(path)

 
def test_unknown_env_var_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA__FOO__BAR", "1")
    with caplog.at_level(logging.WARNING, logger="library.config"):
        load_config(path)
    assert "Environment variable CHEMBL_DA__FOO__BAR ignored" in caplog.text
 
def test_log_level_valid(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("log:\n  level: warn\n")
    cfg = load_config(path)
    assert cfg.log.level == "warn"


def test_log_level_invalid(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("log:\n  level: verbose\n")
    with pytest.raises(ValueError, match="log.level"):
        load_config(path)
 
