from __future__ import annotations

from typing import Any, Iterable, Mapping, MutableMapping

from bioetl.clients.base.normalizers import INormalizer

_VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")


def _safe_get(mapping: Mapping[str, Any], *keys: str) -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _extract_sequence(payload: Mapping[str, Any]) -> str | None:
    raw_value = _safe_get(payload, "sequence", "value")
    if not isinstance(raw_value, str):
        return None

    sequence = raw_value.strip().upper()
    if not sequence:
        return None

    invalid = [aa for aa in sequence if aa not in _VALID_AMINO_ACIDS]
    if invalid:
        filtered = "".join(aa for aa in sequence if aa in _VALID_AMINO_ACIDS)
        return filtered or None
    return sequence


def _extract_annotations(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    features = payload.get("features")
    if not isinstance(features, Iterable):
        return annotations

    for feature in features:
        if not isinstance(feature, Mapping):
            continue

        feature_type = str(feature.get("type") or "").lower()
        if not feature_type:
            continue

        location = feature.get("location") if isinstance(feature, Mapping) else None
        start = _safe_get(location or {}, "start", "value") or _safe_get(
            location or {}, "position", "value"
        )
        end = _safe_get(location or {}, "end", "value") or start

        try:
            start_int = int(start) if start is not None else None
        except (TypeError, ValueError):
            start_int = None
        try:
            end_int = int(end) if end is not None else None
        except (TypeError, ValueError):
            end_int = None

        description = feature.get("description") if isinstance(feature, Mapping) else None
        evidence = feature.get("evidences") if isinstance(feature, Mapping) else None

        mapped_type = _map_feature_type(feature_type)

        annotation: dict[str, Any] = {
            "type": mapped_type,
            "start": start_int,
            "end": end_int,
        }
        if description:
            annotation["description"] = description
        if evidence:
            annotation["evidence"] = evidence

        annotations.append(annotation)

    return annotations


def _map_feature_type(feature_type: str) -> str:
    mapping = {
        "domain": "domain",
        "transmembrane": "transmembrane",
        "intramembrane": "transmembrane",
        "topological domain": "domain",
        "chain": "domain",
        "region": "domain",
        "active site": "active_site",
        "mod_res": "post-translational_modification",
        "lipidation": "post-translational_modification",
        "glycosylation": "post-translational_modification",
        "carbohyd": "post-translational_modification",
        "disulfide bond": "post-translational_modification",
        "ptm": "post-translational_modification",
    }
    return mapping.get(feature_type, feature_type.replace(" ", "_"))


def _extract_genes(payload: Mapping[str, Any]) -> list[str]:
    genes: list[str] = []
    raw_genes = payload.get("genes")
    if not isinstance(raw_genes, Iterable):
        return genes

    for gene in raw_genes:
        if not isinstance(gene, Mapping):
            continue
        gene_name = _safe_get(gene, "geneName", "value")
        if isinstance(gene_name, str):
            genes.append(gene_name)
        synonyms = _safe_get(gene, "synonyms")
        if isinstance(synonyms, Iterable):
            for synonym in synonyms:
                synonym_val = _safe_get(synonym, "value")
                if isinstance(synonym_val, str):
                    genes.append(synonym_val)

    return genes


def _extract_title(payload: Mapping[str, Any]) -> str | None:
    recommended = _safe_get(payload, "proteinDescription", "recommendedName", "fullName", "value")
    if isinstance(recommended, str) and recommended:
        return recommended

    alt_names = _safe_get(payload, "proteinDescription", "alternativeNames")
    if isinstance(alt_names, Iterable):
        for alt in alt_names:
            value = _safe_get(alt, "fullName", "value")
            if isinstance(value, str) and value:
                return value
    return None


def _extract_organism(payload: Mapping[str, Any]) -> str | None:
    organism = _safe_get(payload, "organism", "scientificName")
    if isinstance(organism, str) and organism:
        return organism
    return None


def _extract_length(payload: Mapping[str, Any], sequence: str | None) -> int | None:
    raw_length = _safe_get(payload, "sequence", "length")
    if isinstance(raw_length, int):
        return raw_length
    if sequence:
        return len(sequence)
    return None


class UniprotNormalizerImpl(INormalizer):
    """Нормализатор записей UniProt."""

    def normalize(self, record: Mapping[str, Any]) -> MutableMapping[str, Any]:
        sequence = _extract_sequence(record)

        normalized: MutableMapping[str, Any] = {
            "id": record.get("primaryAccession") or record.get("accession"),
            "title": _extract_title(record),
            "genes": _extract_genes(record),
            "organism": _extract_organism(record),
            "sequence": sequence,
            "length": _extract_length(record, sequence),
            "annotations": _extract_annotations(record),
        }
        return normalized

    def normalize_batch(
        self, records: Iterable[Mapping[str, Any]]
    ) -> Iterable[MutableMapping[str, Any]]:
        for record in records:
            yield self.normalize(record)
