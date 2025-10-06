from __future__ import annotations

from typing import Iterator, Tuple, cast

import pandas as pd
import pytest
import requests

from library.common.log import logger
from library.pipelines.target import chembl_target as ct


@pytest.mark.unit
def test_iter_target_batches_with_retry__shrinks_on_timeout(
    monkeypatch: pytest.MonkeyPatch, cfg, caplog: pytest.LogCaptureFixture
) -> None:
    calls: list[Tuple[Tuple[str, ...], int]] = []

    def fake_iter_target_batches(
        ids: Iterator[str],
        *,
        cfg,
        client,
        mapping_cfg,
        chunk_size: int = 0,
        timeout: float | None = None,
    ):
        ids_list = tuple(ids)
        calls.append((ids_list, chunk_size))
        if chunk_size > 1:
            raise requests.ReadTimeout("simulated timeout")
        payloads = [{"target_chembl_id": ident} for ident in ids_list]
        raw = pd.DataFrame(payloads)
        parsed = pd.DataFrame(payloads)
        yield payloads, raw, parsed

    monkeypatch.setattr(ct, "iter_target_batches", fake_iter_target_batches)

    retry_cfg = cfg.target.chembl.batch_retry
    retry_cfg.enable = True
    retry_cfg.shrink_factor = 0.5
    retry_cfg.min_size = 1

    caplog.set_level("WARNING")

    ids = ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    batches = list(
        ct.iter_target_batches_with_retry(
            ids,
            cfg=cfg.api,
            client=cast(object, None),
            mapping_cfg=cfg.uniprot_mapping,
            chunk_size=2,
            timeout=cfg.target.chembl.timeout,
            retry_cfg=retry_cfg,
            log=logger,
        )
    )

    assert len(batches) == 3
    collected = {
        parsed.iloc[0]["target_chembl_id"]
        for _, _, parsed in batches
        if not parsed.empty
    }
    assert collected == set(ids)

    assert calls[0] == (("CHEMBL1", "CHEMBL2"), 2)
    assert all(size == 1 for _, size in calls[1:])
    assert "chembl_chunk_retry" in caplog.text


@pytest.mark.unit
def test_iter_target_batches_with_retry__disabled_propagates(
    monkeypatch: pytest.MonkeyPatch, cfg
) -> None:
    def failing_iter_target_batches(*args, **kwargs):
        raise requests.ReadTimeout("hard timeout")

    monkeypatch.setattr(ct, "iter_target_batches", failing_iter_target_batches)

    retry_cfg = cfg.target.chembl.batch_retry
    retry_cfg.enable = False

    with pytest.raises(requests.ReadTimeout):
        list(
            ct.iter_target_batches_with_retry(
                ["CHEMBL1"],
                cfg=cfg.api,
                client=cast(object, None),
                mapping_cfg=cfg.uniprot_mapping,
                chunk_size=2,
                timeout=cfg.target.chembl.timeout,
                retry_cfg=retry_cfg,
                log=logger,
            )
        )
