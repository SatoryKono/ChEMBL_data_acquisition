from __future__ import annotations

"""High-level pipeline orchestration."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import pandas as pd

from . import io, mesh_probs, normalize, rules


@dataclass
class Config:
    separators: str = "|;,/"
    delta: float = 0.5
    k_min: int = 3
    unknown_mode: bool = False
    prefer_pubmed_epsilon: float = 0.0


def process_dataframe(df: pd.DataFrame, mesh_map: Dict[str, float], cfg: Config) -> pd.DataFrame:
    """Process ``df`` and return classified dataframe."""
    df = normalize.normalize_publication_types(df, cfg.separators)
    df = normalize.normalize_mesh_fields(df, cfg.separators)
    df = rules.vote_2of3(df)
    needs_mesh = (df["review_votes"] == 1) | df["review_votes"].isna()
    mesh_scores = rules.compute_mesh_scores(
        df[needs_mesh], mesh_map, cfg.prefer_pubmed_epsilon
    )
    for col in mesh_scores.columns:
        df.loc[needs_mesh, col] = mesh_scores[col]
    decisions = rules.assign_labels(
        df, delta=cfg.delta, k_min=cfg.k_min, unknown_mode=cfg.unknown_mode
    )
    for col in decisions.columns:
        df[col] = decisions[col]
    return df


def run_pipeline(
    input_all: Path,
    input_mesh: Path,
    output: Path,
    cfg: Config,
    chunksize: int | None = None,
    output_format: str = "csv",
) -> None:
    """Execute the full classification pipeline."""
    mesh_map = mesh_probs.load_mesh_probabilities(input_mesh)
    reader = io.read_documents(input_all, chunksize=chunksize)
    if chunksize:
        frames = [process_dataframe(chunk, mesh_map, cfg) for chunk in reader]
        result = pd.concat(frames, ignore_index=True)
    else:
        result = process_dataframe(reader, mesh_map, cfg)  # type: ignore[arg-type]
    io.write_output(result, output, fmt=output_format)
