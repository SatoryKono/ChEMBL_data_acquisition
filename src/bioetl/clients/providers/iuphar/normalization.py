from __future__ import annotations

from typing import Any, Iterable, Iterator, Mapping

from bioetl.clients.base.normalizers import INormalizer


class ProviderNormalizer(INormalizer):
    """Нормализатор для ответов IUPHAR (контракт NormalizerProtocol).

    Имплементация размещена в ``providers/iuphar``; фабрика по умолчанию
    находится в ``clients/target/factories.py``.
    """

    def normalize(self, record: Mapping[str, Any]) -> Mapping[str, Any]:
        target_id = record.get("targetId") or record.get("target_id") or record.get("id")
        family_id = record.get("familyId") or record.get("family_id")
        target_name = record.get("name") or record.get("targetName") or record.get("target_name")
        family_name = record.get("familyName") or record.get("family_name")

        normalized = {
            "iuphar_target_id": str(target_id) if target_id is not None else None,
            "iuphar_family_id": str(family_id) if family_id is not None else None,
            "iuphar_type": record.get("type") or record.get("targetType") or record.get("target_type"),
            "iuphar_class": record.get("class") or record.get("targetClass") or record.get("target_class"),
            "iuphar_subclass": record.get("subclass") or record.get("targetSubclass") or record.get("target_subclass"),
            "iuphar_chain": record.get("chain"),
            "iuphar_name": target_name,
        }

        full_id_tokens = [token for token in (target_id, family_id) if token]
        full_name_tokens = [token for token in (target_name, family_name) if token]
        normalized["iuphar_full_id_path"] = "#".join(map(str, full_id_tokens)) if full_id_tokens else None
        normalized["iuphar_full_name_path"] = "#".join(map(str, full_name_tokens)) if full_name_tokens else None

        return normalized

    def normalize_batch(
        self, records: Iterable[Mapping[str, Any]]
    ) -> Iterator[Mapping[str, Any]]:
        for record in records:
            yield self.normalize(record)
