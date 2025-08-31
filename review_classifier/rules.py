from __future__ import annotations

"""Classification rules for review vs non-review."""

import json
from typing import Dict

import pandas as pd

from .normalize import collect_mesh_terms, is_review


def vote_2of3(df: pd.DataFrame) -> pd.DataFrame:
    """Compute review flags and vote counts."""
    df["flag.review.pubmed"] = df["norm.PubMed.PublicationType"].apply(lambda xs: is_review(set(xs)))
    df["flag.review.openalex"] = df["norm.OpenAlex.PublicationTypes"].apply(lambda xs: is_review(set(xs)))
    df["flag.review.scholar"] = df["norm.scholar.PublicationTypes"].apply(lambda xs: is_review(set(xs)))
    df["review_votes"] = (
        df[["flag.review.pubmed", "flag.review.openalex", "flag.review.scholar"]]
        .astype(int)
        .sum(axis=1)
    )
    return df


def compute_mesh_scores(
    df: pd.DataFrame, mesh_map: Dict[str, float], prefer_pubmed_epsilon: float = 0.0
) -> pd.DataFrame:
    """Compute MeSH-based experimental/review scores for rows in ``df``."""

    def _calc(row: pd.Series) -> pd.Series:
        terms = [t for t in collect_mesh_terms(row) if t in mesh_map]
        score_exp = float(sum(mesh_map[t] for t in terms))
        score_review = float(sum(1 - mesh_map[t] for t in terms))
        if row.review_votes == 1 and row["flag.review.pubmed"]:
            score_review += prefer_pubmed_epsilon
        delta = score_exp - score_review
        term_probs = [{"term": t, "p": mesh_map[t]} for t in terms]
        top_terms = sorted(term_probs, key=lambda d: abs(d["p"] - 0.5), reverse=True)[:5]
        top_str = ", ".join(f"{d['term']}:{d['p']:.4f}" for d in top_terms)
        return pd.Series(
            {
                "mesh.terms.used": terms,
                "mesh.term_probs": json.dumps(term_probs, ensure_ascii=False),
                "mesh.top_terms": top_str,
                "score.exp": score_exp,
                "score.review": score_review,
                "score.delta": delta,
                "mesh.k_terms": len(terms),
            }
        )

    return df.apply(_calc, axis=1)


def assign_labels(
    df: pd.DataFrame, *, delta: float, k_min: int, unknown_mode: bool
) -> pd.DataFrame:
    """Assign final labels based on votes and MeSH scores."""

    def _decide(row: pd.Series) -> pd.Series:
        note = ""
        evidence = ""
        if row.review_votes >= 2:
            label = "review"
            rule = "2of3"
            sources = [
                src
                for src, flag in [
                    ("pubmed", row["flag.review.pubmed"]),
                    ("openalex", row["flag.review.openalex"]),
                    ("scholar", row["flag.review.scholar"]),
                ]
                if flag
            ]
            evidence = f"2of3: {', '.join(sources)}"
        elif row.review_votes == 0:
            label = "non-review"
            rule = "2of3"
            evidence = "2of3: "
        else:
            # MeSH refinement
            rule = "mesh_prob_refine"
            delta_score = row.get("score.delta", 0.0)
            k_terms = row.get("mesh.k_terms", 0)
            if k_terms == 0:
                note = "no_mesh_signal"
                if unknown_mode:
                    label = "unknown"
                else:
                    label = "non-review"
            else:
                evidence = (
                    f"mesh_prob_refine: top=[{row.get('mesh.top_terms', '')}];"
                    f" delta={delta_score:+.2f}"
                )
                if delta_score >= delta:
                    label = "non-review"
                elif -delta_score >= delta:
                    label = "review"
                else:
                    if unknown_mode and (abs(delta_score) < delta or k_terms < k_min):
                        label = "unknown"
                        note = "mesh_ambiguous" if k_terms >= k_min else "low_terms"
                    else:
                        label = "non-review"
                        note = "mesh_ambiguous" if k_terms >= k_min else "low_terms"
        return pd.Series(
            {
                "label": label,
                "decision_rule": rule,
                "decision_note": note,
                "evidence": evidence,
            }
        )

    return df.apply(_decide, axis=1)
