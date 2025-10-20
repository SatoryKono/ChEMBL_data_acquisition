"""CLI entry-point for the reproducible target pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

# ruff: noqa: E402 - bootstrap adjusts import order for direct execution
if __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a package module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli
bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module

from library.config import load_config
from library.cli.commands.get_target_data import run_all
from library.utils.logging import get_logger


def _cleanup_intermediate_files(output_dir: Path, final_output: Path) -> None:
    """Удалить промежуточные файлы, оставив только основные результаты."""
    logger = get_logger(__name__)
    
    # Паттерны файлов для удаления
    patterns_to_remove = [
        "*_chembl*.csv",
        "*_chembl*.yaml", 
        "*_chembl*.lock",
        "*_uniprot*.csv",
        "*_uniprot*.yaml",
        "*_uniprot*.lock", 
        "*_iuphar*.csv",
        "*_iuphar*.yaml",
        "*_iuphar*.lock",
        "*_raw*.csv",
        "*_raw*.yaml",
        "*_raw*.lock",
        "output.targets_*.yaml",
        "output.targets_*.lock",
    ]
    
    # Файлы для сохранения
    keep_files = {
        final_output.name,  # Основная таблица
        final_output.with_suffix('.meta.yaml').name,  # Метаданные основной таблицы
        final_output.with_suffix('.meta.yaml.lock').name,  # Lock файл метаданных
        final_output.stem + "_quality_report_table.csv",  # Отчет о качестве
        final_output.stem + "_quality_report_table.csv.meta.yaml",  # Метаданные отчета о качестве
        final_output.stem + "_quality_report_table.csv.meta.yaml.lock",  # Lock файл отчета о качестве
        final_output.stem + "_data_correlation_report_table.csv",  # Отчет о корреляциях
        final_output.stem + "_data_correlation_report_table.csv.meta.yaml",  # Метаданные отчета о корреляциях
        final_output.stem + "_data_correlation_report_table.csv.meta.yaml.lock",  # Lock файл отчета о корреляциях
    }
    
    removed_count = 0
    for pattern in patterns_to_remove:
        for file_path in output_dir.glob(pattern):
            if file_path.name not in keep_files:
                try:
                    file_path.unlink()
                    removed_count += 1
                    logger.debug("removed_intermediate_file", path=str(file_path))
                except OSError as exc:
                    logger.warning("failed_to_remove_file", path=str(file_path), error=str(exc))
    
    logger.info("cleanup_completed", removed_files=removed_count, output_dir=str(output_dir))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Fetch and enrich ChEMBL target data")
    parser.add_argument(
        "command",
        nargs="?",
        default="all",
        help="Execution mode (default: all)",
    )
    parser.add_argument("--limit", type=int, default=1000, help="Number of targets to fetch")
    parser.add_argument("--offset", type=int, default=0, help="Number of targets to skip")
    parser.add_argument(
        "--date-tag",
        type=str,
        default="20250119",
        help="Execution date token in YYYYMMDD format",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory for generated artefacts",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/input"),
        help="Directory containing cached inputs",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        help="Logging level",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to configuration file",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Input CSV file containing target identifiers",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the target pipeline using the full pipeline."""
    args = parse_args(argv)
    logger = get_logger(__name__)

    if args.command != "all":
        logger.error("target_unknown_command", command=args.command)
        return 2

    try:
        config = load_config(args.config) if args.config is not None else load_config(None)
        logger.info("config_loaded", config_type=type(config).__name__)
    except Exception as exc:
        logger.error("target_config_error", error=str(exc))
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.input_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Create namespace with required parameters for run_all()
        namespace = argparse.Namespace(
            input_csv=args.input if args.input else Path("data/input/target.csv"),
            final_out=args.output_dir / f"output.target_{args.date_tag}.csv",
            output_csv=args.output_dir / f"output.target_{args.date_tag}.csv",
            limit=args.limit,
            offset=args.offset,
            command="all",
            chembl_out=None,  # Не сохранять промежуточный ChEMBL файл
            uniprot_out=None,  # Не сохранять промежуточный UniProt файл
            iuphar_out=None,   # Не сохранять промежуточный IUPHAR файл
            raw_out=None,      # Не сохранять raw файл
            raw_format="csv",
            no_reindex_raw=False,
            id_cols=["target_chembl_id"],
            date_tag=args.date_tag,
            disable_gtop=False,
            emit_legacy_artifacts=False,  # Не создавать legacy артефакты
        )
        
        logger.info("target_using_full_pipeline", 
                   input_file=str(namespace.input_csv), 
                   output_file=str(namespace.final_out))
        
        # Run the complete target pipeline
        exit_code = run_all(config, namespace)
        
        if exit_code == 0:
            logger.info("target_pipeline_success", output_file=str(namespace.final_out))
            
            # Очистка промежуточных файлов - оставляем только основные файлы
            _cleanup_intermediate_files(args.output_dir, namespace.final_out)
        else:
            logger.error("target_pipeline_failed", exit_code=exit_code)
        
        return exit_code
    except Exception as exc:
        logger.exception("target_pipeline_failed", error=str(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())