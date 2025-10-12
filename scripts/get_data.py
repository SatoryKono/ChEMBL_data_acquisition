from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence


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
OUTPUT_DIR = PROJECT_ROOT / "data" / "output"
LOGS_DIR = PROJECT_ROOT / "logs"


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
    return args, unknown


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
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
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
