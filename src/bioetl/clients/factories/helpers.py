from __future__ import annotations

from typing import Any, Callable, Mapping, TypeVar

from bioetl.clients.base.http import BaseHttpClientABC
from bioetl.clients.base.normalizers import INormalizer, IdentityNormalizerImpl

TransportBuilder = Callable[[], BaseHttpClientABC]
NormalizerBuilder = Callable[[], INormalizer]

TTransport = TypeVar("TTransport", bound=BaseHttpClientABC)
TNormalizer = TypeVar("TNormalizer", bound=INormalizer)


def _build_component(
    value: Any,
    expected_type: type[TTransport] | type[TNormalizer],
    *,
    default_factory: Callable[[], Any] | None = None,
    error_message: str,
) -> TTransport | TNormalizer:
    if value is None:
        if default_factory is not None:
            default_value = default_factory()
            if isinstance(default_value, expected_type):
                return default_value
        raise ValueError(error_message)

    if not isinstance(value, type) and isinstance(value, expected_type):
        return value

    if callable(value):
        built = value()
        if isinstance(built, expected_type):
            return built

    raise ValueError(error_message)


def build_transport(config: Mapping[str, Any]) -> BaseHttpClientABC:
    """Construct an HTTP transport from config.

    Expected config keys:
    - ``transport``: either a ``BaseHttpClientABC`` instance or a callable returning it.
    """

    return _build_component(
        config.get("transport"),
        BaseHttpClientABC,
        error_message="Transport is not configured or has invalid type",
    )


def build_normalizer(config: Mapping[str, Any]) -> INormalizer:
    """Construct a normalizer from config.

    Expected config keys:
    - ``normalizer``: optional ``INormalizer`` instance or callable returning it.
      If omitted, ``IdentityNormalizerImpl`` is used.
    """

    return _build_component(
        config.get("normalizer"),
        INormalizer,
        default_factory=IdentityNormalizerImpl,
        error_message="Normalizer has invalid type",
    )
