from __future__ import annotations

from typing import Any, Iterable, Mapping, Sequence

from bioetl.clients.base.normalizers import INormalizer


class PubChemNormalizerImpl(INormalizer):
    """Нормализатор данных PubChem."""

    def __init__(self, *, synonym_keys: Sequence[str] | None = None) -> None:
        self._synonym_keys = tuple(synonym_keys or ("synonyms", "Synonym", "Synonyms"))

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        cid = self._extract_cid(record)
        props = self._extract_props(record)

        name = (
            record.get("name")
            or record.get("IUPACName")
            or record.get("Name")
            or props.get("IUPAC Name")
            or props.get("Name")
        )
        smiles = (
            record.get("smiles")
            or record.get("canonical_smiles")
            or record.get("CanonicalSMILES")
            or props.get("Canonical SMILES")
            or props.get("SMILES")
        )
        molecular_formula = (
            record.get("molecular_formula")
            or record.get("MolecularFormula")
            or props.get("Molecular Formula")
        )
        molecular_weight = (
            record.get("molecular_weight")
            or record.get("MolecularWeight")
            or props.get("Molecular Weight")
        )

        synonyms = self._extract_synonyms(record)

        return {
            "id": cid,
            "name": name,
            "smiles": smiles,
            "molecular_formula": molecular_formula,
            "molecular_weight": molecular_weight,
            "synonyms": synonyms,
        }

    def _extract_props(self, record: Mapping[str, Any]) -> dict[str, Any]:
        props: dict[str, Any] = {}
        for prop in record.get("props", []):
            urn = prop.get("urn", {}) if isinstance(prop, Mapping) else {}
            label = urn.get("label") if isinstance(urn, Mapping) else None
            name = urn.get("name") if isinstance(urn, Mapping) else None
            value = prop.get("value", {}) if isinstance(prop, Mapping) else {}
            if not isinstance(value, Mapping):
                continue
            content = value.get("sval") or value.get("fval") or value.get("ival")
            if label:
                props[label] = content
            if name:
                props[name] = props.get(name) or content
            if label and name:
                compound_key = f"{label} {name}"
                props[compound_key] = props.get(compound_key) or content
        return props

    def _extract_cid(self, record: Mapping[str, Any]) -> Any:
        if "cid" in record:
            return record.get("cid")
        if "CID" in record:
            return record.get("CID")

        identifier = record.get("id")
        if isinstance(identifier, Mapping):
            inner = identifier.get("id")
            if isinstance(inner, Mapping):
                return inner.get("cid") or inner.get("CID")
            if isinstance(identifier.get("cid"), (str, int)):
                return identifier.get("cid")
        return identifier

    def _extract_synonyms(self, record: Mapping[str, Any]) -> list[str]:
        for key in self._synonym_keys:
            if key in record and isinstance(record[key], Iterable) and not isinstance(record[key], (str, bytes)):
                return [str(item) for item in record[key]]
        return []
