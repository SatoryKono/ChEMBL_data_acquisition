"""Shared utility helpers for data acquisition pipelines."""

from __future__ import annotations

from importlib import import_module
from typing import Any

__all__ = [
  "ChunkFailureTracker",
  "Logger",
  "LoggerConfig",
  "RateLimiter",
  "SidecarErrors",
  "Stats",
  "compute_backoff_delay",
  "configure_logger",
  "file_sha256",
  "get_global_limiter",
  "get_limiter",
  "log_duration",
  "logger",
  "record_quality_failure",
  "sha256_file",
  "sleep",
  "write_csv_chunks_deterministic",
  "write_csv_deterministic",
  "write_meta_yaml",
]

_MAPPING: dict[str, tuple[str, str]] = {
  "ChunkFailureTracker": ("library.common.fetch_retry", "ChunkFailureTracker"),
  "Logger": ("library.common.logging_setup", "Logger"),
  "LoggerConfig": ("library.common.logging_setup", "LoggerConfig"),
  "RateLimiter": ("library.common.rate_limiter", "RateLimiter"),
  "SidecarErrors": ("library.common.sidecar", "SidecarErrors"),
  "Stats": ("library.common.metadata", "Stats"),
  "compute_backoff_delay": ("library.common.fetch_retry", "compute_backoff_delay"),
  "configure_logger": ("library.common.logging_setup", "configure_logger"),
  "file_sha256": ("library.common.metadata", "file_sha256"),
  "get_global_limiter": ("library.common.rate_limiter", "get_global_limiter"),
  "get_limiter": ("library.common.rate_limiter", "get_limiter"),
  "log_duration": ("library.common.timing", "log_duration"),
  "logger": ("library.common.log", "logger"),
  "record_quality_failure": ("library.common.metadata", "record_quality_failure"),
  "sha256_file": ("library.common.csv_utils", "sha256_file"),
  "sleep": ("library.common.rate_limiter", "sleep"),
  "write_csv_chunks_deterministic": (
    "library.common.csv_utils",
    "write_csv_chunks_deterministic",
  ),
  "write_csv_deterministic": ("library.common.csv_utils", "write_csv_deterministic"),
  "write_meta_yaml": ("library.common.metadata", "write_meta_yaml"),
}


def __getattr__(name: str) -> Any:  # pragma: no cover - thin delegation
  try:
    module_name, attribute = _MAPPING[name]
  except KeyError as exc:  # pragma: no cover - defensive programming
    raise AttributeError(f"module 'library.common' has no attribute {name!r}") from exc
  module = import_module(module_name)
  value = getattr(module, attribute)
  globals()[name] = value
  return value
