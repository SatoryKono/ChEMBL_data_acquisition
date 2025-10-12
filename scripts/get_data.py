from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import import_module
from pathlib import Path
from typing import Iterable, Sequence

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


@dataclass(frozen=True)
class Stage:
    name: str
    script: str


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
LOGS_DIR = DATA_DIR / "logs"


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


def run_stage(stage: Stage, extra_args: Iterable[str]) -> float:
    script_path = SCRIPTS_DIR / stage.script
    if not script_path.exists():
        logging.error("❌ Скрипт %s не найден по пути %s", stage.script, script_path)
        sys.exit(1)

    start = datetime.now()
    logging.info("▶ Запуск %s...", stage.script)
    command = [sys.executable, str(script_path), *extra_args]
    result = subprocess.run(command, check=False)
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


def build_forward_args(args: argparse.Namespace, extra: Sequence[str]) -> list[str]:
    forward: list[str] = []
    if args.limit is not None:
        forward.extend(["--limit", str(args.limit)])
    if args.log_level is not None:
        forward.extend(["--log-level", args.log_level])
    if args.config is not None:
        forward.extend(["--config", str(args.config)])
    forward.extend(extra)

    def _has_option(option: str) -> bool:
        return option in forward

    if not _has_option("--base-path"):
        forward.extend(["--base-path", str(DATA_DIR)])
    if not _has_option("--input-dir"):
        forward.extend(["--input-dir", DEFAULT_INPUT_DIR])
    if not _has_option("--output-dir"):
        forward.extend(["--output-dir", DEFAULT_OUTPUT_DIR])
    return forward


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
