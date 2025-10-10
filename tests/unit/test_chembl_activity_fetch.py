from urllib.parse import parse_qs, urlparse

from library.config import ApiCfg
from library.pipelines.assay.chembl_assay import (
    ACTIVITY_QUERY_FIELDS,
    get_activities,
)


class _RecordingClient:
    def __init__(self) -> None:
        self.urls: list[str] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float
    ) -> dict[str, object]:
        self.urls.append(url)
        return {"activities": [{"activity_id": "ACT1"}], "page_meta": {}}


def test_get_activities__requests_required_fields(monkeypatch) -> None:
    cfg = ApiCfg()
    client = _RecordingClient()

    frame = get_activities(["ACT1"], cfg=cfg, client=client, chunk_size=1)

    assert not frame.empty
    assert "standard_lower_value" in ACTIVITY_QUERY_FIELDS
    assert "standard_upper_value" in ACTIVITY_QUERY_FIELDS

    requested_url = client.urls[0]
    query = parse_qs(urlparse(requested_url).query)
    assert "fields" in query
    requested_fields = query["fields"][0].split(",")

    assert "standard_lower_value" in requested_fields
    assert "standard_upper_value" in requested_fields
    assert {"activity_id", "standard_lower_value", "standard_upper_value"}.issubset(
        set(frame.columns)
    )
