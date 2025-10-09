#!/usr/bin/env python3
"""Utility for converting Markdown documents into DOCX via Pandoc."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path


def check_pandoc() -> None:
    """Ensure that ``pandoc`` is available in ``PATH``."""

    if shutil.which("pandoc"):
        return
    msg = "pandoc не найден в PATH. Установи pandoc и повтори."
    raise RuntimeError(msg)


def build_pandoc_cmd(
    inputs: Sequence[Path],
    output: Path,
    reference_doc: Path | None,
    metadata_file: Path | None,
    number_sections: bool,
    toc: bool,
    toc_depth: int,
    resource_path: list[Path],
    filters: list[str],
    variables: list[str],
    standalone: bool,
) -> list[str]:
    cmd = ["pandoc"]
    # входные файлы (порядок важен)
    for inp in inputs:
        cmd.append(str(inp))

    cmd += ["-o", str(output), "--from=markdown", "--to=docx"]

    if standalone:
        cmd.append("--standalone")

    if reference_doc:
        cmd += ["--reference-doc", str(reference_doc)]

    if metadata_file:
        cmd += ["--metadata-file", str(metadata_file)]

    if number_sections:
        cmd.append("--number-sections")

    if toc:
        cmd.append("--toc")
        cmd += ["--toc-depth", str(toc_depth)]

    # Пути к ресурсам (изображения и т.п.)
    if resource_path:
        # pandoc поддерживает множественные пути через разделитель платформы
        sep = ";" if sys.platform.startswith("win") else ":"
        cmd += ["--resource-path", sep.join(str(p) for p in resource_path)]

    # Фильтры (напр., pandoc-crossref, mermaid-filter, citeproc)
    for f in filters:
        cmd += ["--filter", f]

    # Переменные шаблона (key=value)
    for var in variables:
        cmd += ["-V", var]

    return cmd


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="convert_md_to_docx",
        description="Конвертация Markdown → DOCX через pandoc с поддержкой шаблонов и метаданных.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Входные .md файлы (один или несколько, склеиваются по порядку).",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Путь к выходному .docx"
    )
    parser.add_argument(
        "--refdoc",
        "--reference-doc",
        dest="reference_doc",
        type=Path,
        help="DOCX-шаблон стилей (reference.docx).",
    )
    parser.add_argument(
        "--meta",
        "--metadata-file",
        dest="metadata_file",
        type=Path,
        help="YAML-файл с метаданными (title, author, date, ...).",
    )
    parser.add_argument(
        "--number-sections", action="store_true", help="Нумеровать разделы."
    )
    parser.add_argument(
        "--toc", action="store_true", help="Добавить таблицу содержимого."
    )
    parser.add_argument(
        "--toc-depth", type=int, default=3, help="Глубина TOC (по умолчанию 3)."
    )
    parser.add_argument(
        "--resource-path",
        nargs="*",
        type=Path,
        default=[],
        help="Список директорий для ресурсов (картинки и т.п.).",
    )
    parser.add_argument(
        "--filter",
        dest="filters",
        nargs="*",
        default=[],
        help="Список pandoc-фильтров (например: pandoc-crossref, pandoc-mermaid, citeproc).",
    )
    parser.add_argument(
        "-V",
        "--variable",
        dest="variables",
        nargs="*",
        default=[],
        help="Переменные шаблона (key=value), прокидываются в pandoc (-V).",
    )
    parser.add_argument(
        "--no-standalone",
        action="store_true",
        help="Не добавлять --standalone (по умолчанию standalone включён).",
    )

    args = parser.parse_args(argv)

    try:
        check_pandoc()
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 127
    for p in args.inputs:
        if not p.exists():
            print(f"ERROR: входной файл не найден: {p}", file=sys.stderr)
            return 2

    if args.reference_doc and not args.reference_doc.exists():
        print(f"ERROR: reference DOCX не найден: {args.reference_doc}", file=sys.stderr)
        return 2

    if args.metadata_file and not args.metadata_file.exists():
        print(f"ERROR: metadata YAML не найден: {args.metadata_file}", file=sys.stderr)
        return 2

    cmd = build_pandoc_cmd(
        inputs=args.inputs,
        output=args.output,
        reference_doc=args.reference_doc,
        metadata_file=args.metadata_file,
        number_sections=args.number_sections,
        toc=args.toc,
        toc_depth=args.toc_depth,
        resource_path=args.resource_path,
        filters=args.filters,
        variables=args.variables,
        standalone=not args.no_standalone,
    )

    print("Running:", " ".join(str(x) for x in cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: pandoc завершился с ошибкой (exit {exc.returncode}).",
            file=sys.stderr,
        )
        return int(exc.returncode)
    except FileNotFoundError as exc:
        print(f"ERROR: не удалось запустить pandoc: {exc}", file=sys.stderr)
        return 127

    # немножко гуманности
    print(f"OK: {args.output} создан.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
