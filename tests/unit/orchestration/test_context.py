"""Tests for :mod:`library.orchestration.context`."""

from __future__ import annotations

import pytest

from library.orchestration import ETLContext


@pytest.mark.unit
def test_etl_context_reuses_single_client(stub_etl_context):
    context, clients = stub_etl_context

    with context as active:
        first = active.chembl_client
        second = active.chembl_client

    assert first is second
    assert clients
    assert clients[0].close_calls == 1


@pytest.mark.unit
def test_etl_context_recreates_client_after_close(stub_etl_context):
    context, clients = stub_etl_context

    with context as active:
        first = active.chembl_client

    with context as active:
        second = active.chembl_client

    assert first is not second
    assert len(clients) == 2
    assert all(client.close_calls == 1 for client in clients)


@pytest.mark.unit
def test_etl_context_global_limiter(cfg):
    cfg.rate.global_rps = 5
    cfg.rate.global_burst = 7

    with ETLContext(cfg) as context:
        limiter = context.global_limiter
        assert limiter is not None
        assert limiter.rps == cfg.rate.global_rps
        assert limiter.burst == cfg.rate.global_burst
        assert context.global_limiter is limiter

    with ETLContext(cfg) as reopened:
        assert reopened.global_limiter is limiter


def test_etl_context_custom_cleanup(stub_etl_context):
    context, clients = stub_etl_context

    cleanup_calls: list[str] = []

    with context as active:
        active.register_cleanup(lambda: cleanup_calls.append("cleanup"))
        _ = active.chembl_client

    assert cleanup_calls == ["cleanup"]
    assert clients and clients[0].close_calls == 1


@pytest.mark.unit
def test_etl_context_register_cleanup_after_close(stub_etl_context):
    context, _ = stub_etl_context

    with context:
        pass

    with pytest.raises(RuntimeError):
        context.register_cleanup(lambda: None)
