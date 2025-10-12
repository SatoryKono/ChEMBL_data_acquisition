from __future__ import annotations

import argparse
import logging
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import import_module
from pathlib import Path
from typing import Collection, Sequence

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


def _guard_cli_module() -> None:
    """Ensure the core ``library.cli.commands.get_data`` module imports cleanly."""

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


from library.config import DEFAULT_CONFIG_PATH, load_config
from library.config.env import _default_base_path as _config_default_base_path


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
    Stage("testitem", "get_testitem_data.py"),
    Stage("target", "get_target_data.py"),
    Stage("document", "get_document_data.py"),
    Stage("assay", "get_assay_data.py"),
    Stage("activity", "get_activity_data.py"),
)

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_INPUT_DIR = "input"
DEFAULT_OUTPUT_DIR = "output"
OUTPUT_DIR = DATA_DIR / DEFAULT_OUTPUT_DIR
LOGS_DIR = PROJECT_ROOT / "logs"
_PUBCHEM_ENV_VAR = "CHEMBL_DA_PUBCHEM_ENABLE"
_BASE_PATH_ENV_VAR = "CHEMBL_DA_BASE_PATH"


def _extract_option_value(args: Sequence[str], option: str) -> str | None:
    """Return the value for ``option`` (supports ``--opt value`` and ``--opt=value``)."""

    prefixed = f"{option}="
    for idx, token in enumerate(args):
        if token == option:
            return args[idx + 1] if idx + 1 < len(args) else ""
        if token.startswith(prefixed):
            return token[len(prefixed) :]
    return None


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
        config_path = Path(config_value)
    else:
        config_path = DEFAULT_CONFIG_PATH
    base_value = _extract_option_value(args, "--base-path")
    base_path = Path(base_value).expanduser() if base_value else None
    return Path(config_path).expanduser(), base_path


def _pubchem_enabled_from_config(args: Sequence[str]) -> bool | None:
    config_path, base_path = _resolve_config_location(args)
    try:
        config = load_config(config_path, base_path=base_path)
    except Exception as exc:  # pragma: no cover - defensive logging
        logging.debug(
            "Не удалось загрузить конфигурацию %s: %s", config_path, exc,
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

    if env_state is True or (config_enabled is True and env_state is not False):
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
        return

    candidate = Path(base_path_value).expanduser()
    if not candidate.is_absolute():
        candidate = (Path.cwd() / candidate).resolve()
    else:
        candidate = candidate.resolve()

    env[_BASE_PATH_ENV_VAR] = str(candidate)


def _resolve_forward_base_path(
    args: argparse.Namespace, forwarded_extras: Sequence[str]
) -> Path:
    """Return the base path propagated to pipeline subprocesses."""

    tokens: list[str] = []
    if args.config is not None:
        tokens.extend(["--config", str(args.config)])
    tokens.extend(forwarded_extras)

    config_path, cli_base_path = _resolve_config_location(tokens)
    if cli_base_path is not None:
        return cli_base_path

    try:
        config = load_config(config_path, base_path=None)
    except Exception as exc:  # pragma: no cover - defensive logging
        logging.debug(
            "Не удалось вычислить базовый путь из %s: %s", config_path, exc
        )
    else:
        local_cfg = getattr(config, "local", None)
        io_cfg = getattr(local_cfg, "io", None) if local_cfg is not None else None
        output_dir = getattr(io_cfg, "output_dir", None) if io_cfg is not None else None
        if output_dir is not None:
            output_path = Path(output_dir).expanduser()
            if output_path.is_absolute():
                return output_path.resolve().parent

    return _config_default_base_path()


def parse_args(argv: Sequence[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
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


def _coerce_forward_args(
    forward_args: ForwardArgs | Sequence[str],
) -> ForwardArgs:
    if isinstance(forward_args, ForwardArgs):
        return forward_args
    return ForwardArgs(tuple(forward_args), extras_start=len(forward_args), extra_len=0)


def run_stage(stage: Stage, forward_args: ForwardArgs | Sequence[str]) -> float:
    script_path = SCRIPTS_DIR / stage.script
    if not script_path.exists():
        logging.error("❌ Скрипт %s не найден по пути %s", stage.script, script_path)
        sys.exit(1)

    forward = _coerce_forward_args(forward_args)

    if not isinstance(forward_args, ForwardArgs):
        forward_args = ForwardArgs(tuple(forward_args), 0, len(forward_args))

    if stage.name == "target":
        stage_args = forward.with_default_subcommand(
            "all", choices=TARGET_SUBCOMMANDS
        )
    else:
        stage_args = forward.as_list()

    if stage.name == "testitem":
        # Ensure PubChem enrichment mirrors direct CLI invocation of the
        # get_testitem_data.py script.
        if "--pubchem-enable" not in stage_args:
            stage_args.append("--pubchem-enable")

    env = os.environ.copy()
    _ensure_base_path_env(stage_args, env)
    if stage.name == "testitem":
        _ensure_pubchem_env(stage_args, env)

    command = [sys.executable, str(script_path), *stage_args]
    quoted_command = shlex.join(command)
    logging.info("▶ Запуск %s", stage.script)
    logging.info("   Команда: %s", quoted_command)
    logging.info("   Рабочая директория: %s", os.getcwd())

    start = datetime.now()
    result = subprocess.run(command, check=False, env=env)
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
        return option in forward

    if not _has_option("--base-path"):
        base_path = _resolve_forward_base_path(args, forwarded_extras)
        forward.extend(["--base-path", str(base_path)])
    if not _has_option("--input-dir"):
        forward.extend(["--input-dir", DEFAULT_INPUT_DIR])
    if not _has_option("--output-dir"):
        forward.extend(["--output-dir", DEFAULT_OUTPUT_DIR])
    return ForwardArgs(tuple(forward), extras_start, extra_len)


def log_summary(durations: list[tuple[str, float]]) -> None:
    if not durations:
        logging.info("Все этапы были пропущены по флагу --skip.")
        return

    logging.info("⏱️ Сводка по длительности этапов:")
    for name, value in durations:
        logging.info(" • %s: %.1f сек.", name, value)


def count_output_files() -> int:
    if not OUTPUT_DIR.exists():
        return 0
    return sum(1 for path in OUTPUT_DIR.glob("*.csv"))


def main(argv: Sequence[str] | None = None) -> int:
    args, unknown = parse_args(argv)
    log_file = configure_logging(args.log_level)
    logging.info("Логи сохраняются в %s", log_file)

    forward_args = build_forward_args(args, unknown)

    durations: list[tuple[str, float]] = []
    for stage in STAGES:
        if stage.name in args.skip:
            logging.info("⏭ Пропуск этапа %s по флагу --skip", stage.name)
            continue
        duration = run_stage(stage, forward_args)
        durations.append((stage.name, duration))

    log_summary(durations)

    csv_count = count_output_files()
    logging.info(
        "🎉 Все выбранные этапы завершены. Найдено %d CSV-файлов в %s.",
        csv_count,
        OUTPUT_DIR,
    )

    if csv_count != 15:
        logging.warning(
            "Ожидалось получить 15 CSV-файлов, фактически найдено %d.",
            csv_count,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
