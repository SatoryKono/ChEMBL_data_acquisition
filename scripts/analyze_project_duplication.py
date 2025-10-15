"""Audit duplication between an installed project tree and repository sources."""
from __future__ import annotations

import argparse
import csv
import hashlib
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable


MAX_SIMILARITY_SAMPLE_BYTES = 200_000
SIMILARITY_THRESHOLD = 0.90


@dataclass
class FileRecord:
    path: Path
    rel_path: str
    size: int
    sha256: str
    sample_text: str | None


def iter_files(root: Path) -> Iterable[Path]:
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if "__pycache__" in path.parts:
            continue
        yield path


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def maybe_read_text(path: Path) -> str | None:
    try:
        size = path.stat().st_size
    except OSError:
        return None
    if size > MAX_SIMILARITY_SAMPLE_BYTES:
        return None
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None


def build_records(base: Path, *, rel_base: Path | None = None) -> list[FileRecord]:
    records: list[FileRecord] = []
    rel_base = rel_base or base
    for path in iter_files(base):
        rel = path.relative_to(rel_base).as_posix()
        size = path.stat().st_size
        sha256 = compute_sha256(path)
        sample_text = maybe_read_text(path)
        records.append(FileRecord(path=path, rel_path=rel, size=size, sha256=sha256, sample_text=sample_text))
    return records


def find_best_match(project: FileRecord, source_by_rel: dict[str, FileRecord], source_by_hash: dict[str, list[FileRecord]]):
    match_type = "unique"
    similarity = 0.0
    source_record: FileRecord | None = None

    rel_candidate = source_by_rel.get(project.rel_path)
    if rel_candidate is not None:
        if rel_candidate.sha256 == project.sha256:
            return "exact_path_content", 1.0, rel_candidate
        # Only attempt similarity for reasonably small text files
        if project.sample_text and rel_candidate.sample_text:
            similarity = SequenceMatcher(None, project.sample_text, rel_candidate.sample_text).ratio()
            if similarity >= SIMILARITY_THRESHOLD:
                return "path_similar", similarity, rel_candidate
        match_type = "path_mismatch"
        similarity = similarity or 0.0
        source_record = rel_candidate

    if project.sha256 in source_by_hash:
        candidates = source_by_hash[project.sha256]
        # Exclude candidate already considered to avoid reusing same file object
        for candidate in candidates:
            if source_record and candidate.path == source_record.path:
                continue
            return "content_match", 1.0, candidate

    if match_type != "unique":
        return match_type, similarity, source_record
    return "unique", 0.0, None


def write_csv(
    destination: Path,
    mapping: list[tuple[FileRecord, str, float, FileRecord | None]],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "project_path",
                "project_size_bytes",
                "project_sha256",
                "match_type",
                "similarity",
                "source_path",
                "source_size_bytes",
                "source_sha256",
            ]
        )
        for project, match_type, similarity, source in mapping:
            writer.writerow(
                [
                    project.rel_path,
                    project.size,
                    project.sha256,
                    match_type,
                    f"{similarity:.4f}",
                    source.rel_path if source else "",
                    source.size if source else "",
                    source.sha256 if source else "",
                ]
            )


def summarise(mapping: list[tuple[FileRecord, str, float, FileRecord | None]]):
    total_files = len(mapping)
    summary: dict[str, int] = {}
    total_sizes: dict[str, int] = {}
    dir_counts: dict[str, int] = {}
    dir_sizes: dict[str, int] = {}
    for project, match_type, _, _ in mapping:
        summary[match_type] = summary.get(match_type, 0) + 1
        total_sizes[match_type] = total_sizes.get(match_type, 0) + project.size
        top_level = project.rel_path.split("/", 1)[0]
        dir_counts[top_level] = dir_counts.get(top_level, 0) + 1
        dir_sizes[top_level] = dir_sizes.get(top_level, 0) + project.size
    return {
        "total_files": total_files,
        "match_type_counts": summary,
        "match_type_sizes": total_sizes,
        "dir_counts": dir_counts,
        "dir_sizes": dir_sizes,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit reconstructed project tree against repository sources.")
    parser.add_argument("--project", type=Path, default=Path("project"), help="Path to the reconstructed project directory")
    parser.add_argument(
        "--sources",
        type=Path,
        nargs="*",
        default=[Path("library"), Path("scripts"), Path("config")],
        help="Repository source directories to compare against",
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd(), help="Repository root for relative paths")
    parser.add_argument("--csv", type=Path, default=Path("reports/project_dup_map.csv"), help="Destination CSV path")
    parser.add_argument("--summary", type=Path, default=Path("reports/project_dup_summary.json"), help="Summary JSON output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project.resolve()
    repo_root = args.repo_root.resolve()

    project_records = build_records(project_root)

    source_records: list[FileRecord] = []
    for src in args.sources:
        root = (repo_root / src).resolve()
        if not root.exists():
            continue
        source_records.extend(build_records(root, rel_base=repo_root))

    source_by_rel = {record.rel_path: record for record in source_records}
    source_by_hash: dict[str, list[FileRecord]] = {}
    for record in source_records:
        source_by_hash.setdefault(record.sha256, []).append(record)

    mapping: list[tuple[FileRecord, str, float, FileRecord | None]] = []
    for record in sorted(project_records, key=lambda r: r.rel_path):
        match_type, similarity, source = find_best_match(record, source_by_rel, source_by_hash)
        mapping.append((record, match_type, similarity, source))

    write_csv(args.csv, mapping)

    # Write summary JSON for downstream documentation
    import json

    summary_payload = summarise(mapping)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
