# ruff: noqa: E402, I001

from __future__ import annotations

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import import_module
from pathlib import Path, PureWindowsPath

# NOTE:
#   ``python scripts/get_data.py`` executed from Windows adds ``scripts``
#   directory to ``sys.path`` but not the project root.  The import below then
#   fails because Python looks for ``scripts`` package inside the ``scripts``
#   directory itself.  Import the helper lazily and, on failure, prepend the
#   project root so the package resolves when the script is executed directly.
try:  # pragma: no cover - exercised via CLI invocation
    from scripts._bootstrap import ensure_project_root
except ModuleNotFoundError as exc:  # pragma: no cover - import guard
    if exc.name not in {"scripts", "scripts._bootstrap"}:
        raise
    current_dir = Path(__file__).resolve().parent
    project_root = current_dir.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    sys.modules.pop("scripts", None)
    from scripts._bootstrap import ensure_project_root

ensure_project_root(__file__)

from library.config import DEFAULT_CONFIG_PATH, load_config  # noqa: E402
from library.io.path_utils import OUTPUT_DIR as LIB_OUTPUT_DIR  # noqa: E402
from library.io.path_utils import ROOT as LIB_ROOT  # noqa: E402


def _guard_cli_module() -> None:
    """Ensure the core ``library.cli.commands.get_data`` module imports cleanly."""

    ensure_project_root(__file__)

    try:
        import_module("library.cli.commands.get_data")
    except ModuleNotFoundError as exc:
        if exc.name != "library":
            raise
        ensure_project_root(__file__)
        import_module("library.cli.commands.get_data")
    except SyntaxError as exc:
        location = f"{exc.filename}:{exc.lineno}"
        raise SystemExit(f"merge conflict detected in {location}") from exc


_guard_cli_module()


@dataclass(frozen=True)
class Stage:
    name: str
    script: str


@dataclass(frozen=True)
class ForwardArgs:
    """Tokenised command line forwarded to pipeline scripts."""

    tokens: tuple[str, ...]
    extras_start: int
    extra_len: int
    output_dir: Path

    @property
    def extras_end(self) -> int:
        """Return the exclusive end offset for forwarded extras."""

        return self.extras_start + self.extra_len

    def as_list(self) -> list[str]:
        """Return a mutable copy of the token sequence."""

        return list(self.tokens)

    def with_default_subcommand(
        self, default_command: str, *, choices: Collection[str]
    ) -> list[str]:
        """Ensure a recognised sub-command is present in the extras slice."""

        tokens = self.as_list()
        extras = tokens[self.extras_start : self.extras_end]
        if any(token in choices for token in extras):
            return tokens
        tokens.insert(self.extras_start, default_command)
        return tokens


STAGES: tuple[Stage, ...] = (
    Stage("document", "get_document_data.py"),
    Stage("target", "get_target_data.py"),
    Stage("assay", "get_assay_data.py"),
    Stage("testitem", "get_testitem_data.py"),
    Stage("activity", "get_activity_data.py"),
)

STAGE_DISPLAY_NAMES: dict[str, str] = {
    "document": "Document",
    "target": "Target",
    "assay": "Assay",
    "testitem": "Test item",
    "activity": "Activity",
}

STAGE_SEQUENCE_LABEL = " → ".join(
    STAGE_DISPLAY_NAMES[stage.name] for stage in STAGES
)

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = LIB_ROOT
DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_INPUT_DIR = Path("data") / "input"
DEFAULT_OUTPUT_DIR = Path("data") / "output"
OUTPUT_DIR = LIB_OUTPUT_DIR
LOGS_DIR = DATA_DIR / "logs"
_PUBCHEM_ENV_VAR = "CHEMBL_DA_PUBCHEM_ENABLE"
_BASE_PATH_ENV_VAR = "CHEMBL_DA_BASE_PATH"

_STAGE_EXPECTED_OUTPUTS: dict[str, int] = {
    "document": 3,
    "target": 3,
    "assay": 3,
    "testitem": 3,
    "activity": 3,
}

_EXPECTED_CSV_COUNT = sum(
    _STAGE_EXPECTED_OUTPUTS[stage.name] for stage in STAGES
)

_CLEANUP_FILE_PATTERNS: tuple[str, ...] = (
    "*.lock",
    "*.pkl.lock",
    "*.tmp",
    "*.tmp.*",
    "*_intermediate.*",
    "*_debug.log",
    "cache_*.pkl",
)

_CLEANUP_DIRECTORY_NAMES: frozenset[str] = frozenset({
    "raw",
    "tmp",
    "interim",
    "cache",
})


def _has_explicit_pubchem_flag(tokens: Sequence[str]) -> bool:
    """Return ``True`` when PubChem CLI overrides are already present."""

    for token in tokens:
        if token in {"--pubchem-enable", "--no-pubchem-enable"}:
            return True
        if token.startswith("--pubchem-enable="):
            return True
    return False


def _extract_option_value(args: Sequence[str], option: str) -> str | None:
    """Return the value for ``option`` (supports ``--opt value`` and ``--opt=value``)."""

    prefixed = f"{option}="
    for idx, token in enumerate(args):
        if token == option:
            return args[idx + 1] if idx + 1 < len(args) else ""
        if token.startswith(prefixed):
            return token[len(prefixed) :]
    return None


def _normalise_cli_path(value: str) -> Path:
    """Return a :class:`Path` handling Windows-style separators on POSIX hosts."""

    stripped = value.strip()
    if os.name != "nt" and "\\" in stripped:
        windows_path = PureWindowsPath(stripped)
        return Path(windows_path.as_posix())
    return Path(stripped)


def _resolve_forward_output_dir(tokens: Sequence[str]) -> Path:
    """Return the absolute output directory encoded within ``tokens``."""

    raw_output = _extract_option_value(tokens, "--output-dir")
    if not raw_output:
        raw_output = os.fspath(DEFAULT_OUTPUT_DIR)
    output_candidate = _normalise_cli_path(raw_output)

    if output_candidate.is_absolute():
        return output_candidate.resolve()

    raw_base = _extract_option_value(tokens, "--base-path")
    if raw_base:
        base_candidate = _normalise_cli_path(raw_base)
        if not base_candidate.is_absolute():
            base_candidate = (PROJECT_ROOT / base_candidate).resolve()
        else:
            base_candidate = base_candidate.resolve()
    else:
        base_candidate = PROJECT_ROOT

    return (base_candidate / output_candidate).resolve()


def _has_flag(tokens: Sequence[str], option: str) -> bool:
    """Return ``True`` when *option* is present in *tokens*."""

    prefix = f"{option}="
    for token in tokens:
        if token == option or token.startswith(prefix):
            return True
    return False


def _normalize_env_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    return None


def _resolve_config_location(args: Sequence[str]) -> tuple[Path, Path | None]:
    config_value = _extract_option_value(args, "--config")
    if config_value:
        config_path = _normalise_cli_path(config_value)
        if not config_path.is_absolute():
            config_path = (PROJECT_ROOT / config_path).resolve()
        else:
            config_path = config_path.resolve()
    else:
        config_path = DEFAULT_CONFIG_PATH
    base_value = _extract_option_value(args, "--base-path")
    if base_value:
        base_candidate = _normalise_cli_path(base_value)
        if not base_candidate.is_absolute():
            base_path = (PROJECT_ROOT / base_candidate).resolve()
        else:
            base_path = base_candidate.resolve()
    else:
        base_path = None
    return config_path, base_path


def _pubchem_enabled_from_config(args: Sequence[str]) -> bool | None:
    config_path, base_path = _resolve_config_location(args)
    try:
        config = load_config(config_path, base_path=base_path)
    except Exception as exc:  # pragma: no cover - defensive logging
        logging.debug(
            "Не удалось загрузить конфигурацию %s: %s",
            config_path,
            exc,
        )
        return None

    sources = getattr(config, "sources", None)
    if sources is None:
        return None
    pubchem_cfg = getattr(sources, "pubchem", None)
    if pubchem_cfg is None:
        return None
    return getattr(pubchem_cfg, "enable", None)


def _ensure_pubchem_env(args: Sequence[str], env: dict[str, str]) -> None:
    """Ensure PubChem enrichment is enabled for the test item subprocess."""

    env_state = _normalize_env_bool(env.get(_PUBCHEM_ENV_VAR))
    config_enabled = _pubchem_enabled_from_config(args)

    if env_state is False:
        return

    if config_enabled is False:
        return

    if env_state is True:
        logging.info("[PUBCHEM] Enrichment enabled for testitem table")
        return

    if config_enabled is True:
        env[_PUBCHEM_ENV_VAR] = "true"
        logging.info("[PUBCHEM] Enrichment enabled for testitem table")
        return

    env[_PUBCHEM_ENV_VAR] = "true"
    logging.info("[PUBCHEM] Enrichment enabled for testitem table")


def _ensure_base_path_env(args: Sequence[str], env: dict[str, str]) -> None:
    """Populate ``CHEMBL_DA_BASE_PATH`` for subprocesses when unset."""

    current_value = env.get(_BASE_PATH_ENV_VAR)
    if current_value:
        return

    base_path_value = _extract_option_value(args, "--base-path")
    if not base_path_value:
        env[_BASE_PATH_ENV_VAR] = os.fspath(PROJECT_ROOT)
        return

    candidate = _normalise_cli_path(base_path_value)
    if not candidate.is_absolute():
        candidate = (PROJECT_ROOT / candidate).resolve()
    else:
        candidate = candidate.resolve()

    env[_BASE_PATH_ENV_VAR] = os.fspath(candidate)


def _resolve_forward_base_path(
    args: argparse.Namespace,
    forwarded_extras: Sequence[str],
    *,
    fallback: Path | str | None = None,
) -> Path:
    """Return the base path propagated to pipeline subprocesses."""

    tokens: list[str] = []
    if args.config is not None:
        tokens.extend(["--config", str(args.config)])
    tokens.extend(forwarded_extras)

    config_path, cli_base_path = _resolve_config_location(tokens)
    if cli_base_path is not None:
        return cli_base_path

    if fallback is None:
        default_base_path = PROJECT_ROOT
    else:
        fallback_path = _normalise_cli_path(os.fspath(fallback))
        if not fallback_path.is_absolute():
            default_base_path = (PROJECT_ROOT / fallback_path).resolve()
        else:
            default_base_path = fallback_path.resolve()

    try:
        config = load_config(config_path, base_path=None)
    except Exception as exc:  # pragma: no cover - defensive logging
        logging.debug("Не удалось вычислить базовый путь из %s: %s", config_path, exc)
    else:
        local_cfg = getattr(config, "local", None)
        io_cfg = getattr(local_cfg, "io", None) if local_cfg is not None else None
        output_dir = getattr(io_cfg, "output_dir", None) if io_cfg is not None else None
        if output_dir is not None:
            output_path = _normalise_cli_path(os.fspath(output_dir))
            if not output_path.is_absolute():
                output_path = (PROJECT_ROOT / output_path).resolve()
            else:
                output_path = output_path.resolve()
            return output_path.parent

    return default_base_path


def parse_args(
    argv: Sequence[str] | None = None,
) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(
        description="Оркестратор последовательного запуска скриптов сбора данных ChEMBL.",
    )
    parser.add_argument(
        "--skip",
        nargs="+",
        choices=[stage.name for stage in STAGES],
        default=(),
        help="Пропустить перечисленные этапы (например: --skip assay activity)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Ограничение числа обрабатываемых записей, прокидывается во все подпроцессы.",
    )
    parser.add_argument(
        "--log-level",
        dest="log_level",
        type=_parse_log_level,
        default=None,
        help="Уровень логирования (для оркестратора и подпроцессов).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Путь до пользовательского конфигурационного файла.",
    )

    args, unknown = parser.parse_known_args(argv)

    extras = list(unknown)

    def _looks_like_positional_limit(token: str, index: int) -> bool:
        if not token:
            return False
        if token.startswith("-"):
            return False
        if index > 0 and extras[index - 1].startswith("-"):
            return False
        return token.isdigit()

    positional_index: int | None = None
    for idx, candidate in enumerate(extras):
        if _looks_like_positional_limit(candidate, idx):
            positional_index = idx
            break

    if positional_index is not None:
        value_token = extras[positional_index]
        if len(extras) > 1:
            parser.error(
                "Неизвестный позиционный аргумент "
                f"'{value_token}'. Используйте '--limit {value_token}'."
            )
        if args.limit is not None:
            parser.error(
                "Ограничение уже указано через '--limit'. "
                "Удалите позиционный аргумент или используйте только флаг."
            )
        args.limit = int(value_token, 10)
        extras.pop(positional_index)

    return args, extras


def _parse_log_level(value: str) -> str:
    normalized = value.upper()
    if normalized not in logging._nameToLevel:
        valid = ", ".join(sorted(logging._nameToLevel))
        raise argparse.ArgumentTypeError(
            f"Недопустимый уровень логирования '{value}'. Доступные значения: {valid}"
        )
    return normalized


def configure_logging(level_name: str | None) -> Path:
    level = logging.INFO
    if level_name is not None:
        level = logging._nameToLevel.get(level_name, logging.INFO)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    log_file = LOGS_DIR / f"get_data_{timestamp}.log"

    handlers: list[logging.Handler] = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, encoding="utf-8"),
    ]

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )
    logging.getLogger(__name__).debug("Логирование настроено на уровень %s", level)
    return log_file


TARGET_SUBCOMMANDS: tuple[str, ...] = ("uniprot", "chembl", "iuphar", "all")
DOCUMENT_SUBCOMMANDS: tuple[str, ...] = ("chembl", "pubmed", "all")

_DEFAULT_TARGET_SUBCOMMAND = "all"
_DEFAULT_DOCUMENT_SUBCOMMAND = "all"


def _coerce_forward_args(
    forward_args: ForwardArgs | Sequence[str],
) -> ForwardArgs:
    if isinstance(forward_args, ForwardArgs):
        return forward_args
    tokens = tuple(forward_args)
    return ForwardArgs(
        tokens,
        extras_start=len(tokens),
        extra_len=0,
        output_dir=_resolve_forward_output_dir(tokens),
    )


def run_stage(stage: Stage, forward_args: ForwardArgs | Sequence[str]) -> float:
    script_path = SCRIPTS_DIR / stage.script
    if not script_path.exists():
        logging.error("❌ Скрипт %s не найден по пути %s", stage.script, script_path)
        sys.exit(1)

    forward = _coerce_forward_args(forward_args)

    if not isinstance(forward_args, ForwardArgs):
        tokens = tuple(forward_args)
        forward_args = ForwardArgs(
            tokens,
            extras_start=0,
            extra_len=len(tokens),
            output_dir=_resolve_forward_output_dir(tokens),
        )

    if stage.name == "target":
        stage_args = forward.with_default_subcommand(
            _DEFAULT_TARGET_SUBCOMMAND, choices=TARGET_SUBCOMMANDS
        )
    elif stage.name == "document":
        stage_args = forward.with_default_subcommand(
            _DEFAULT_DOCUMENT_SUBCOMMAND, choices=DOCUMENT_SUBCOMMANDS
        )
    else:
        stage_args = forward.as_list()

    if stage.name == "testitem":
        # Ensure PubChem enrichment mirrors direct CLI invocation of the
        # get_testitem_data.py script.
        if not _has_explicit_pubchem_flag(stage_args):
            stage_args.append("--pubchem-enable")

    env = os.environ.copy()
    _ensure_base_path_env(stage_args, env)
    if stage.name == "testitem":
        _ensure_pubchem_env(stage_args, env)

    command = [sys.executable, str(script_path), *stage_args]
    quoted_command = shlex.join(command)
    logging.info("▶ Запуск %s", stage.script)
    logging.info("   Команда: %s", quoted_command)
    logging.info("   Рабочая директория: %s", os.fspath(PROJECT_ROOT))

    start = datetime.now()
    result = subprocess.run(command, check=False, env=env, cwd=PROJECT_ROOT)
    duration = (datetime.now() - start).total_seconds()

    if result.returncode != 0:
        logging.error(
            "❌ Ошибка при выполнении %s (exit code %s)",
            stage.script,
            result.returncode,
        )
        sys.exit(result.returncode)

    logging.info("✅ %s завершён успешно за %.1f сек.", stage.script, duration)
    return duration


def build_forward_args(args: argparse.Namespace, extra: Sequence[str]) -> ForwardArgs:
    forward: list[str] = []
    forwarded_extras = list(extra)
    if args.limit is not None:
        forwarded_extras = ["--limit", str(args.limit), *forwarded_extras]
    if args.log_level is not None:
        forward.extend(["--log-level", args.log_level])
    if args.config is not None:
        forward.extend(["--config", str(args.config)])
    extras_start = len(forward)
    forward.extend(forwarded_extras)
    extra_len = len(forwarded_extras)

    def _has_option(option: str) -> bool:
        prefixed = f"{option}="
        return any(
            token == option or token.startswith(prefixed)
            for token in forward
        )

    # if not _has_option("--base-path"):
    #    base_path = _resolve_forward_base_path(args, forwarded_extras)
    #    forward.extend(["--base-path", str(base_path)])
    if not _has_option("--input-dir"):
        forward.extend(["--input-dir", os.fspath(DEFAULT_INPUT_DIR)])
    if not _has_option("--output-dir"):
        forward.extend(["--output-dir", os.fspath(DEFAULT_OUTPUT_DIR)])

    tokens = tuple(forward)
    output_dir = _resolve_forward_output_dir(tokens)
    return ForwardArgs(tokens, extras_start, extra_len, output_dir)


def log_summary(durations: list[tuple[str, float]]) -> None:
    if not durations:
        logging.info("Все этапы были пропущены по флагу --skip.")
        return

    logging.info("⏱️ Сводка по длительности этапов:")
    for name, value in durations:
        logging.info(" • %s: %.1f сек.", name, value)


def count_output_files(output_dir: Path) -> int:
    resolved = output_dir.resolve()
    if not resolved.exists():
        return 0
    return sum(1 for path in resolved.glob("*.csv"))


def _resolve_output_directory(
    args: argparse.Namespace, forward_args: ForwardArgs
) -> Path:
    """Return the absolute output directory referenced by the orchestrator."""

    # ``ForwardArgs`` already carries a normalised output directory derived from
    # the forwarded command line tokens.  Rely on this value directly to avoid
    # drifting away from the paths that the pipeline stages actually use.
    return forward_args.output_dir.resolve()


def _is_cleanup_directory(path: Path) -> bool:
    """Return ``True`` when ``path`` should be removed during cleanup."""

    name = path.name.lower()
    return name in _CLEANUP_DIRECTORY_NAMES


def _relative_to_output(path: Path, output_dir: Path) -> str:
    try:
        return os.fspath(path.relative_to(output_dir))
    except ValueError:
        return os.fspath(path)


def _remove_file(path: Path, *, output_dir: Path) -> bool:
    try:
        path.unlink()
    except FileNotFoundError:
        return False
    except OSError as exc:  # pragma: no cover - filesystem variance
        logging.warning("[CLEANUP] Не удалось удалить %s: %s", path, exc)
        return False
    logging.info("[CLEANUP] Удалён файл %s", _relative_to_output(path, output_dir))
    return True


def _remove_directory(path: Path, *, output_dir: Path) -> bool:
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        return False
    except OSError as exc:  # pragma: no cover - filesystem variance
        logging.warning("[CLEANUP] Не удалось удалить каталог %s: %s", path, exc)
        return False
    logging.info("[CLEANUP] Удалён каталог %s", _relative_to_output(path, output_dir))
    return True


def cleanup_intermediate_files(output_dir: Path) -> int:
    """Remove temporary and diagnostic artefacts from ``output_dir``."""

    resolved = output_dir.resolve()
    if not resolved.exists():
        logging.debug(
            "[CLEANUP] Пропуск очистки: каталог %s не найден",
            resolved,
        )
        return 0

    removed = 0
    seen_files: set[Path] = set()
    for pattern in _CLEANUP_FILE_PATTERNS:
        for candidate in sorted(resolved.rglob(pattern)):
            candidate = candidate.resolve()
            if candidate in seen_files or not candidate.is_file():
                continue
            if _remove_file(candidate, output_dir=resolved):
                removed += 1
            seen_files.add(candidate)

    directories: set[Path] = set()
    for candidate in resolved.rglob("*"):
        if candidate.is_dir() and _is_cleanup_directory(candidate):
            directories.add(candidate.resolve())

    # Remove deepest directories first to avoid orphaned parents lingering.
    for directory in sorted(directories, key=lambda item: len(item.parts), reverse=True):
        if _remove_directory(directory, output_dir=resolved):
            removed += 1

    return removed


def _should_run_cleanup(forward_args: ForwardArgs) -> bool:
    """Return ``True`` when the orchestrator should prune intermediate artefacts."""

    tokens = forward_args.tokens
    debug_requested = _has_flag(tokens, "--debug")
    keep_requested = _has_flag(tokens, "--keep-intermediate")
    legacy_requested = _has_flag(tokens, "--emit-legacy-artifacts")
    return not (debug_requested or keep_requested or legacy_requested)


def main(argv: Sequence[str] | None = None) -> int:
    args, unknown = parse_args(argv)
    log_file = configure_logging(args.log_level)
    logging.info("Логи сохраняются в %s", log_file)
    logging.info("Последовательность этапов: %s", STAGE_SEQUENCE_LABEL)

    forward_args = build_forward_args(args, unknown)
    resolved_output_dir = _resolve_output_directory(args, forward_args)

    durations: list[tuple[str, float]] = []
    skipped_stages = set(args.skip)

    for stage in STAGES:
        if stage.name in skipped_stages:
            logging.info("⏭ Пропуск этапа %s по флагу --skip", stage.name)
            continue
        duration = run_stage(stage, forward_args)
        durations.append((stage.name, duration))

    log_summary(durations)

    csv_count = count_output_files(resolved_output_dir)
    logging.info(
        "🎉 Все выбранные этапы завершены. Найдено %d CSV-файлов в %s.",
        csv_count,
        resolved_output_dir,
    )

    if not skipped_stages and csv_count != _EXPECTED_CSV_COUNT:
        logging.warning(
            "Ожидалось получить %d CSV-файлов (по 3 на этап: %s), фактически найдено %d.",
            _EXPECTED_CSV_COUNT,
            STAGE_SEQUENCE_LABEL,
            csv_count,
        )

    if _should_run_cleanup(forward_args):
        removed = cleanup_intermediate_files(resolved_output_dir)
        logging.info("[CLEANUP] Завершено: удалено %d артефакт(ов)", removed)
    else:
        logging.info(
            "[CLEANUP] Пропущено: сохранение промежуточных файлов запрошено пользователем"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
