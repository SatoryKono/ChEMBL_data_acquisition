"""Batch processing pipeline for review classification."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

import pandas as pd

from .normalization import normalize_mesh, normalize_publication_types
from .scoring import Record, score_record


def _make_uid(ids: Mapping[str, str]) -> str:
    """Combine available identifiers into a single uid."""
    for key in ("doi", "pmid", "pmcid", "openalexid"):
        if ids.get(key):
            return ids[key]
    # fallback to joined ids
    return "|".join(f"{k}:{v}" for k, v in ids.items())


def process_records(records: Sequence[Mapping[str, object]], output_dir: Path, logger: logging.Logger) -> None:
    """Process an iterable of raw records and persist artifacts.

    Parameters
    ----------
    records:
        Sequence of raw record dictionaries.
    output_dir:
        Base directory where artifacts are written.
    logger:
        Logger for decision trail output.
    """
    raw_dir = output_dir / "raw"
    norm_dir = output_dir / "normalized"
    feat_dir = output_dir / "features"
    score_dir = output_dir / "scores"
    decision_dir = output_dir / "decisions"
    log_dir = output_dir / "logs"

    for d in [raw_dir, norm_dir, feat_dir, score_dir, decision_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    raw_path = raw_dir / "records.ndjson"
    with raw_path.open("w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")

    pt_rows: List[Dict[str, object]] = []
    mesh_rows: List[Dict[str, object]] = []
    score_rows: List[Dict[str, object]] = []
    label_rows: List[Dict[str, object]] = []

    for rec in records:
        ids = rec.get("ids", {})
        uid = _make_uid(ids)
        raw_pt = rec.get("publication_types", {})
        raw_mesh = rec.get("mesh_terms", {})
        raw_desc = {k: v.get("descriptors", []) for k, v in raw_mesh.items()}
        raw_qual = {k: v.get("qualifiers", []) for k, v in raw_mesh.items()}

        norm_pt = normalize_publication_types(raw_pt)
        norm_desc, norm_qual = normalize_mesh(raw_desc, raw_qual)

        for src, pts in norm_pt.items():
            for pt in pts:
                pt_rows.append({"uid": uid, "source": src, "pub_type": pt})
        for src, descs in norm_desc.items():
            for dsc in descs:
                mesh_rows.append({"uid": uid, "kind": "descriptor", "value": dsc})
        for src, quals in norm_qual.items():
            for q in quals:
                mesh_rows.append({"uid": uid, "kind": "qualifier", "value": q})

        record_obj = Record(ids=ids, pts=norm_pt, mesh_desc=norm_desc, mesh_qual=norm_qual)
        result = score_record(record_obj, logger)

        score_rows.append(
            {
                "uid": uid,
                "score_review": result.score_review,
                "score_non_review": result.score_non_review,
            }
        )
        label_rows.append(
            {
                "uid": uid,
                "label": result.label,
                "reason": result.reason,
            }
        )

    pd.DataFrame(pt_rows).to_parquet(norm_dir / "publication_type.parquet")
    pd.DataFrame(mesh_rows).to_parquet(feat_dir / "mesh.parquet")
    pd.DataFrame(score_rows).to_parquet(score_dir / "scores.parquet")
    pd.DataFrame(label_rows).to_parquet(decision_dir / "labels.parquet")

    run_meta = {
        "rules": "v1",
        "records": len(records),
    }
    (output_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2), encoding="utf-8")
