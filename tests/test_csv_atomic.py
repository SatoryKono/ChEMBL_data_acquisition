from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import library.common.csv_utils as csv_utils
from library.common.csv_utils import write_csv_chunks_deterministic


def test_write_csv_chunks_failure_preserves_original(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "out.csv"
    path.write_text("sentinel\n", encoding="utf-8-sig")

    df = pd.DataFrame({"id": [1], "value": [2]})

    original_open_atomic = csv_utils.open_atomic

    @contextmanager
    def failing_open_atomic(*args: Any, **kwargs: Any) -> Any:
        with original_open_atomic(*args, **kwargs) as handle:

            class FaultyFile:
                def __init__(self, wrapped: Any) -> None:
                    self._wrapped = wrapped
                    self._writes = 0

                def write(self, data: Any) -> Any:
                    self._writes += 1
                    if self._writes < 2:
                        return self._wrapped.write(data)
                    raise RuntimeError("boom")

                def __getattr__(self, name: str) -> Any:
                    return getattr(self._wrapped, name)

            yield FaultyFile(handle)

    monkeypatch.setattr(csv_utils, "open_atomic", failing_open_atomic)

    with pytest.raises(RuntimeError, match="boom"):
        write_csv_chunks_deterministic(
            [df],
            path,
            key_cols=["id"],
        )

    assert path.read_text(encoding="utf-8-sig") == "sentinel\n"
