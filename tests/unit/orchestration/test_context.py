"""Tests for :mod:`library.orchestration.context`."""

from __future__ import annotations

import pytest


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
